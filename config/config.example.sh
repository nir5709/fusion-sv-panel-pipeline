# Example configuration for run_fusion_sv_pipeline.sh
# Copy this file and edit paths for your environment.

###############################################################################
# ENVIRONMENTS
###############################################################################
SVPANEL_ENV="svpanel"
MANTA_ENV="manta_py2"
ANNOTSV_ENV="base"

###############################################################################
# INPUT / OUTPUT PATHS
###############################################################################
READ_DIR="/path/to/input_fastq"
WORKROOT="/path/to/workdir/fusion_sv_pipeline"
REF="/path/to/reference/hg19.fa"
BED="/path/to/reference/design.hg19.bed"
FUSION_TSV="/path/to/fusion_targets.tsv"
GENE_PANEL="/path/to/panel_genes.txt"
ANNOTSVDIR="/path/to/annotsv_db_hg19"

###############################################################################
# OPTIONAL RESOURCES
###############################################################################
DO_BQSR=1
DBSNP="/path/to/dbsnp_138.hg19.vcf.gz"
GRIDSS_BLACKLIST="/path/to/hg19-blacklist.v2.bed"

###############################################################################
# THREADS
###############################################################################
THREADS=16
FASTQC_THREADS=4
FASTP_THREADS=8
GRIDSS_THREADS=8
GRIDSS_HEAP="10g"

###############################################################################
# TOOL BINARIES (override only if not in PATH)
###############################################################################
FASTP="fastp"
FASTQC="fastqc"
BWA="bwa"
SAMTOOLS="samtools"
BCFTOOLS="bcftools"
BEDTOOLS="bedtools"
GATK="gatk"
GRIDSS="gridss"
DELLY="delly"
SVABA="svaba"
BGZIP="bgzip"
TABIX="tabix"
