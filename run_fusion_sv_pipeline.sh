#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Fusion / SV Detection Pipeline (5-callers)
# Author: Nihar Garg
#
# Overview:
#   Targeted paired-end FASTQ -> QC -> trimming -> alignment -> dedup ->
#   optional BQSR -> Manta / GRIDSS / DELLY / SvABA -> AnnotSV ->
#   panel-gene filtering -> breakpoint evidence summary
#
# Usage:
#   bash run_fusion_sv_pipeline.sh config/config.example.sh
#
# Notes:
#   - Tested for hg19 / GRCh37-style resources.
#   - This workflow is intended for research / assay development use.
#   - Clinical deployment requires local validation and SOP compliance.
###############################################################################

usage() {
  cat <<EOF
Usage:
  $0 <config.sh>

Example:
  $0 config/config.example.sh
EOF
}

[[ $# -eq 1 ]] || { usage; exit 1; }

CONFIG="$1"
[[ -f "$CONFIG" ]] || { echo "[FATAL] Config not found: $CONFIG" >&2; exit 1; }

# shellcheck source=/dev/null
source "$CONFIG"

###############################################################################
# PRE-RUN CHECKS
###############################################################################
if ! command -v conda >/dev/null 2>&1; then
  echo "[FATAL] conda not found in PATH" >&2
  exit 1
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$SVPANEL_ENV"

###############################################################################
# TOOL DEFAULTS (override in config if needed)
###############################################################################
FASTP="${FASTP:-fastp}"
FASTQC="${FASTQC:-fastqc}"
BWA="${BWA:-bwa}"
SAMTOOLS="${SAMTOOLS:-samtools}"
BCFTOOLS="${BCFTOOLS:-bcftools}"
BEDTOOLS="${BEDTOOLS:-bedtools}"
GATK="${GATK:-gatk}"
GRIDSS="${GRIDSS:-gridss}"
DELLY="${DELLY:-delly}"
SVABA="${SVABA:-svaba}"
BGZIP="${BGZIP:-bgzip}"
TABIX="${TABIX:-tabix}"

###############################################################################
# HELPERS
###############################################################################
log(){ echo "[$(date '+%F %T')] $*"; }
need(){ [[ -s "$1" ]] || { echo "[FATAL] Missing: $1" >&2; exit 1; }; }

mean_depth(){
  local bam="$1" region="$2"
  "$SAMTOOLS" depth -a -r "$region" "$bam" | awk '{s+=$3;n++} END{ if(n) printf("%.3f", s/n); else printf("0"); }'
}

sa_bridge_count(){
  local bam="$1" regionA="$2" partner_chr="$3" partner_pos="$4" window="$5"
  "$SAMTOOLS" view "$bam" "$regionA" \
    | awk -v chr="$partner_chr" -v pos="$partner_pos" -v w="$window" '
      BEGIN{FS="\t"; c=0}
      {
        for(i=12;i<=NF;i++){
          if($i ~ /^SA:Z:/){
            sa=$i; sub(/^SA:Z:/,"",sa)
            n=split(sa,arr,";")
            for(j=1;j<=n;j++){
              if(arr[j]=="") continue
              m=split(arr[j],f,",")
              if(m>=2 && f[1]==chr){
                p=f[2]+0
                if(p>=pos-w && p<=pos+w){ c++; break }
              }
            }
          }
        }
      }
      END{print c}
    '
}

ensure_rg(){
  local inbam="$1"
  local outbam="$2"
  local sample="$3"

  if "$SAMTOOLS" view -H "$inbam" | grep -q '^@RG'; then
    echo "$inbam"
    return
  fi

  log "[RG] No @RG in BAM. Adding read groups..."
  "$GATK" AddOrReplaceReadGroups \
    -I "$inbam" -O "$outbam" \
    -RGID "$sample" -RGSM "$sample" \
    -RGLB "lib1" -RGPL "ILLUMINA" -RGPU "unit1" \
    --CREATE_INDEX true
  echo "$outbam"
}

ensure_bai(){
  local bam="$1"
  if [[ -s "${bam}.bai" ]]; then
    return
  elif [[ -s "${bam%.bam}.bai" ]]; then
    ln -sf "${bam%.bam}.bai" "${bam}.bai"
  else
    log "[INDEX] Creating BAM index: $bam"
    "$SAMTOOLS" index "$bam"
  fi
}

run_fastqc(){
  local outdir="$1"; shift
  mkdir -p "$outdir"
  conda run -n "$SVPANEL_ENV" bash -lc "
    set -euo pipefail
    $FASTQC -t ${FASTQC_THREADS} -o '$outdir' $*
  "
}

run_fastp_pe(){
  local sample="$1" r1="$2" r2="$3" outdir="$4"
  mkdir -p "$outdir"

  local p1="$outdir/${sample}.R1.trimmed.fastq.gz"
  local p2="$outdir/${sample}.R2.trimmed.fastq.gz"
  local html="$outdir/${sample}.fastp.html"
  local json="$outdir/${sample}.fastp.json"

  if [[ -s "$p1" && -s "$p2" ]]; then
    echo -e "$p1\t$p2"
    return
  fi

  echo "[$(date '+%F %T')] [TRIM] fastp PE -> $outdir" >&2

  conda run -n "$SVPANEL_ENV" bash -lc "
    set -euo pipefail
    $FASTP \
      -i '$r1' -I '$r2' \
      -o '$p1' -O '$p2' \
      --thread ${FASTP_THREADS} \
      --detect_adapter_for_pe \
      --qualified_quality_phred 20 \
      --unqualified_percent_limit 40 \
      --n_base_limit 5 \
      --length_required 30 \
      --trim_front1 0 --trim_front2 0 \
      --trim_tail1 0  --trim_tail2 0 \
      --html '$html' --json '$json'
  "

  echo -e "$p1\t$p2"
}

run_annotsv_panel(){
  local invcf="$1"
  local outprefix="$2"
  local panel="$3"
  local outtsv="${outprefix}.AnnotSV.tsv"
  local outpanel="${outprefix}.PANEL.tsv"

  [[ -s "$invcf" ]] || { log "[AnnotSV] missing VCF: $invcf"; return 0; }
  need "$panel"

  if [[ ! -s "$outtsv" ]]; then
    log "[AnnotSV] annotate -> $outtsv"
    conda run -n "$ANNOTSV_ENV" bash -lc "
      set -euo pipefail
      AnnotSV \
        -genomeBuild GRCh37 \
        -annotationsDir '$ANNOTSVDIR' \
        -SVinputFile '$invcf' \
        -SVinputInfo 1 \
        -outputFile '$outprefix'
    "
  fi

  if [[ -s "$outtsv" && ! -s "$outpanel" ]]; then
    log "[AnnotSV] panel filter -> $outpanel"
    awk 'BEGIN{FS=OFS="\t"}
         NR==FNR{
           g=toupper($1); sub(/^[ \t]+|[ \t]+$/,"",g);
           if(g!="" && g!~ /^#/) genes[g]=1;
           next
         }
         FNR==1{print; next}
         {
           row=toupper($0);
           hit=0;
           for(g in genes){
             if(row ~ ("(^|[^A-Z0-9_])" g "([^A-Z0-9_]|$)")){ hit=1; break }
           }
           if(hit) print
         }' "$panel" "$outtsv" > "$outpanel"
  fi
}

###############################################################################
# REQUIRED CONFIG VARIABLES
###############################################################################
: "${READ_DIR:?Missing READ_DIR in config}"
: "${WORKROOT:?Missing WORKROOT in config}"
: "${REF:?Missing REF in config}"
: "${BED:?Missing BED in config}"
: "${FUSION_TSV:?Missing FUSION_TSV in config}"
: "${GENE_PANEL:?Missing GENE_PANEL in config}"
: "${ANNOTSVDIR:?Missing ANNOTSVDIR in config}"
: "${THREADS:?Missing THREADS in config}"
: "${SVPANEL_ENV:?Missing SVPANEL_ENV in config}"
: "${MANTA_ENV:?Missing MANTA_ENV in config}"
: "${ANNOTSV_ENV:?Missing ANNOTSV_ENV in config}"
: "${DO_BQSR:?Missing DO_BQSR in config}"
: "${GRIDSS_BLACKLIST:?Missing GRIDSS_BLACKLIST in config}"

if [[ "$DO_BQSR" -eq 1 ]]; then
  : "${DBSNP:?Missing DBSNP in config while DO_BQSR=1}"
fi

need "$REF"
need "${REF}.fai"
need "$BED"
need "$FUSION_TSV"
need "$GENE_PANEL"
need "$ANNOTSVDIR"
need "$GRIDSS_BLACKLIST"

mkdir -p "$WORKROOT"/{summary,logs}

SUMMARY="$WORKROOT/summary/summary.tsv"
echo -e "sample\tfusion\tlbp\trbp\tfinal_bam\tmeanDepth_L400\tmeanDepth_R400\tSA_LtoR_500\tSA_RtoL_500\tmanta_bnd_hits\tgridss_bnd_hits\tdelly_bnd_hits\tsvaba_bnd_hits\tmanta_rundir" > "$SUMMARY"

###############################################################################
# MAIN LOOP
###############################################################################
tail -n +2 "$FUSION_TSV" | while IFS=$'\t' read -r SAMPLE FUSION LBP RBP; do
  log "======================================================"
  log "SAMPLE=$SAMPLE  FUSION=$FUSION  LBP=$LBP  RBP=$RBP"

  SAMPLEDIR="$WORKROOT/$SAMPLE"
  mkdir -p "$SAMPLEDIR"/{qc/raw,qc/trimmed,trim,bam,tmp,sv/manta,sv/gridss,sv/delly,sv/svaba,annotsv,logs}

  R1=$(ls -1 "$READ_DIR"/"${SAMPLE}"*_R1*.fastq.gz 2>/dev/null | head -n 1 || true)
  R2=$(ls -1 "$READ_DIR"/"${SAMPLE}"*_R2*.fastq.gz 2>/dev/null | head -n 1 || true)
  if [[ -z "${R1}" || -z "${R2}" ]]; then
    log "[FATAL] FASTQs not found for $SAMPLE in $READ_DIR"
    exit 1
  fi
  log "[OK] R1=$R1"
  log "[OK] R2=$R2"

  run_fastqc "$SAMPLEDIR/qc/raw" "$R1" "$R2"

  read -r TR1 TR2 < <(run_fastp_pe "$SAMPLE" "$R1" "$R2" "$SAMPLEDIR/trim")
  log "[OK] TR1=$TR1"
  log "[OK] TR2=$TR2"

  run_fastqc "$SAMPLEDIR/qc/trimmed" "$TR1" "$TR2"

  SORTED="$SAMPLEDIR/bam/${SAMPLE}.sorted.bam"
  DEDUP="$SAMPLEDIR/bam/${SAMPLE}.dedup.bam"
  METRICS="$SAMPLEDIR/bam/${SAMPLE}.markdup.metrics.txt"
  DEDUP_RG="$SAMPLEDIR/bam/${SAMPLE}.dedup.RG.bam"
  RECAL="$SAMPLEDIR/bam/${SAMPLE}.recal.table"
  BQSR_BAM="$SAMPLEDIR/bam/${SAMPLE}.bqsr.bam"

  if [[ ! -s "$SORTED" ]]; then
    log "[BWA] Align + sort (trimmed reads, with @RG)..."
    "$BWA" mem -t "$THREADS" \
      -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:lib1\tPL:ILLUMINA\tPU:unit1" \
      "$REF" "$TR1" "$TR2" \
      | "$SAMTOOLS" sort -@ "$THREADS" -o "$SORTED" -
    "$SAMTOOLS" index "$SORTED"
  else
    log "[BWA] Sorted BAM exists; skipping"
  fi

  if [[ ! -s "$DEDUP" ]]; then
    log "[GATK] MarkDuplicates"
    "$GATK" MarkDuplicates \
      -I "$SORTED" -O "$DEDUP" -M "$METRICS" \
      --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT
  else
    log "[GATK] Dedup BAM exists; skipping"
  fi
  ensure_bai "$DEDUP"

  FINAL_IN=$(ensure_rg "$DEDUP" "$DEDUP_RG" "$SAMPLE")
  ensure_bai "$FINAL_IN"
  FINAL_BAM="$FINAL_IN"

  if [[ "$DO_BQSR" -eq 1 ]]; then
    need "$DBSNP"
    log "[BQSR] BaseRecalibrator + ApplyBQSR"
    "$GATK" BaseRecalibrator -R "$REF" -I "$FINAL_IN" --known-sites "$DBSNP" -O "$RECAL"
    "$GATK" ApplyBQSR -R "$REF" -I "$FINAL_IN" --bqsr-recal-file "$RECAL" -O "$BQSR_BAM"
    "$SAMTOOLS" index "$BQSR_BAM"
    FINAL_BAM="$BQSR_BAM"
  fi

  ts=$(date '+%Y%m%d_%H%M%S')
  MANTA_RUNDIR="$SAMPLEDIR/sv/manta/${SAMPLE}_nocallRegions_${ts}"
  mkdir -p "$MANTA_RUNDIR"

  if [[ ! -s "$MANTA_RUNDIR/results/variants/candidateSV.vcf.gz" ]]; then
    log "[MANTA] Running -> $MANTA_RUNDIR"
    conda run -n "$MANTA_ENV" bash -lc "
      set -euo pipefail
      configManta.py --bam '$FINAL_BAM' --referenceFasta '$REF' --runDir '$MANTA_RUNDIR'
      '$MANTA_RUNDIR/runWorkflow.py' -m local -j '$THREADS'
    "
  else
    log "[MANTA] Outputs exist; skipping"
  fi
  MANTA_VCF="$MANTA_RUNDIR/results/variants/diploidSV.vcf.gz"

  GRIDSS_VCF="$SAMPLEDIR/sv/gridss/${SAMPLE}.gridss.vcf.gz"
  GRIDSS_ASM="$SAMPLEDIR/sv/gridss/${SAMPLE}.gridss.assembly.bam"
  GRIDSS_THREADS="${GRIDSS_THREADS:-8}"
  GRIDSS_HEAP="${GRIDSS_HEAP:-10g}"

  if [[ ! -s "$GRIDSS_VCF" ]]; then
    log "[GRIDSS] Running (threads=$GRIDSS_THREADS heap=$GRIDSS_HEAP)"
    "$GRIDSS" \
      -r "$REF" \
      -o "$GRIDSS_VCF" \
      -a "$GRIDSS_ASM" \
      -w "$SAMPLEDIR/sv/gridss/work" \
      -t "$GRIDSS_THREADS" \
      --jvmheap "$GRIDSS_HEAP" \
      -b "$GRIDSS_BLACKLIST" \
      "$FINAL_BAM"
    "$TABIX" -p vcf "$GRIDSS_VCF" || true
  else
    log "[GRIDSS] VCF exists; skipping"
  fi

  DELLY_BCF="$SAMPLEDIR/sv/delly/${SAMPLE}.delly.bcf"
  DELLY_VCF="$SAMPLEDIR/sv/delly/${SAMPLE}.delly.vcf.gz"
  if [[ ! -s "$DELLY_VCF" ]]; then
    log "[DELLY] Running (tumor-only)"
    "$DELLY" call -g "$REF" -o "$DELLY_BCF" "$FINAL_BAM"
    "$BCFTOOLS" view -Oz -o "$DELLY_VCF" "$DELLY_BCF"
    "$TABIX" -p vcf "$DELLY_VCF" || true
  else
    log "[DELLY] VCF exists; skipping"
  fi

  SVABA_PREFIX="$SAMPLEDIR/sv/svaba/${SAMPLE}"
  SVABA_VCF="${SVABA_PREFIX}.svaba.sv.vcf"
  if [[ ! -s "$SVABA_VCF" ]]; then
    log "[SvABA] Running"
    "$SVABA" run -t "$FINAL_BAM" -G "$REF" -a "$SVABA_PREFIX" -p "$THREADS"
  else
    log "[SvABA] VCF exists; skipping"
  fi

  run_annotsv_panel "$MANTA_VCF"  "$SAMPLEDIR/annotsv/${SAMPLE}.manta"  "$GENE_PANEL"
  run_annotsv_panel "$GRIDSS_VCF" "$SAMPLEDIR/annotsv/${SAMPLE}.gridss" "$GENE_PANEL"
  run_annotsv_panel "$DELLY_VCF"  "$SAMPLEDIR/annotsv/${SAMPLE}.delly"  "$GENE_PANEL"
  run_annotsv_panel "$SVABA_VCF"  "$SAMPLEDIR/annotsv/${SAMPLE}.svaba"  "$GENE_PANEL"

  mdL="0"; mdR="0"; saLR="0"; saRL="0"
  manta_hit="0"; grid_hit="0"; delly_hit="0"; svaba_hit="0"

  if [[ "$LBP" != "NA" && "$RBP" != "NA" ]]; then
    LCHR=${LBP%:*}; LPOS=${LBP#*:}
    RCHR=${RBP%:*}; RPOS=${RBP#*:}

    L1=$((LPOS-200)); L2=$((LPOS+200))
    R1P=$((RPOS-200)); R2P=$((RPOS+200))
    mdL=$(mean_depth "$FINAL_BAM" "${LCHR}:${L1}-${L2}")
    mdR=$(mean_depth "$FINAL_BAM" "${RCHR}:${R1P}-${R2P}")

    saLR=$(sa_bridge_count "$FINAL_BAM" "${LCHR}:${L1}-${L2}" "$RCHR" "$RPOS" 500)
    saRL=$(sa_bridge_count "$FINAL_BAM" "${RCHR}:${R1P}-${R2P}" "$LCHR" "$LPOS" 500)

    manta_hit=$("$BCFTOOLS" view -H "$MANTA_VCF" 2>/dev/null \
      | awk -F'\t' -v lc="$LCHR" -v lp="$LPOS" -v rc="$RCHR" -v rp="$RPOS" '
        function near(x,p){return (x>=p-2000 && x<=p+2000)}
        $8~/(^|;)SVTYPE=BND(;|$)/ && ( ($1==lc && near($2,lp)) || ($1==rc && near($2,rp)) || ($5~lc":" && $5~rc":") ) {c++}
        END{print c+0}'
    )

    grid_hit=$("$BCFTOOLS" view -H "$GRIDSS_VCF" 2>/dev/null \
      | awk -F'\t' -v lc="$LCHR" -v lp="$LPOS" -v rc="$RCHR" -v rp="$RPOS" '
        function near(x,p){return (x>=p-2000 && x<=p+2000)}
        $8~/(^|;)SVTYPE=BND(;|$)/ && ( ($1==lc && near($2,lp)) || ($1==rc && near($2,rp)) || ($5~lc":" && $5~rc":") ) {c++}
        END{print c+0}'
    )

    delly_hit=$("$BCFTOOLS" view -H "$DELLY_VCF" 2>/dev/null \
      | awk -F'\t' -v lc="$LCHR" -v lp="$LPOS" -v rc="$RCHR" -v rp="$RPOS" '
        function near(x,p){return (x>=p-2000 && x<=p+2000)}
        $8~/(^|;)SVTYPE=BND(;|$)/ && ( ($1==lc && near($2,lp)) || ($1==rc && near($2,rp)) || ($5~lc":" && $5~rc":") ) {c++}
        END{print c+0}'
    )

    svaba_hit=$("$BCFTOOLS" view -H "$SVABA_VCF" 2>/dev/null \
      | awk -F'\t' -v lc="$LCHR" -v lp="$LPOS" -v rc="$RCHR" -v rp="$RPOS" '
        function near(x,p){return (x>=p-2000 && x<=p+2000)}
        $8~/(^|;)SVTYPE=BND(;|$)/ && ( ($1==lc && near($2,lp)) || ($1==rc && near($2,rp)) || ($5~lc":" && $5~rc":") ) {c++}
        END{print c+0}'
    )
  else
    log "[INFO] LBP/RBP=NA; skipping breakpoint evidence checks for this sample."
  fi

  echo -e "${SAMPLE}\t${FUSION}\t${LBP}\t${RBP}\t${FINAL_BAM}\t${mdL}\t${mdR}\t${saLR}\t${saRL}\t${manta_hit}\t${grid_hit}\t${delly_hit}\t${svaba_hit}\t${MANTA_RUNDIR}" >> "$SUMMARY"
  log "[DONE] Sample finished -> appended to summary"
done

log "======================================================"
log "ALL DONE. Summary table: $SUMMARY"
log "Panel-filtered AnnotSV outputs per sample are in: $WORKROOT/<sample>/annotsv/*.PANEL.tsv"
