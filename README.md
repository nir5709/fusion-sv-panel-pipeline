# fusion-sv-panel-pipeline

A targeted fusion and structural variant detection workflow for paired-end sequencing data using multiple SV callers and downstream panel-focused annotation.

## What this pipeline does

This workflow processes paired-end FASTQ files and performs:

1. raw FASTQ quality control with FastQC
2. read trimming with fastp
3. trimmed FASTQ quality control
4. BWA-MEM alignment
5. BAM sorting and duplicate marking
6. optional BQSR with GATK
7. SV detection using:
   - Manta
   - GRIDSS2
   - DELLY
   - SvABA
8. AnnotSV-based annotation
9. panel-gene filtering for must-report events
10. breakpoint evidence summarization around expected fusion loci

## Intended use

This repository is designed for research workflows, assay development, and validated internal analysis pipelines.

It is not a substitute for clinical-grade reporting by itself. Any diagnostic or accredited-lab use should be supported by local validation, QC thresholds, SOPs, and reporting governance.

## Highlights

- multi-caller SV strategy
- expected breakpoint evidence support
- panel-focused AnnotSV filtering
- per-sample structured output directories
- summary table for fusion evidence review
- portable config-driven setup instead of hardcoded machine paths

## Repository structure

```text
fusion-sv-panel-pipeline/
├── run_fusion_sv_pipeline.sh
├── README.md
├── LICENSE
├── .gitignore
├── config/
│   └── config.example.sh
├── templates/
│   ├── fusion_targets.example.tsv
│   └── panel_genes.example.txt
├── docs/
│   └── workflow.md
└── env/
    └── tool_versions.txt
```

## Inputs

### 1. Paired-end FASTQ files
The script expects sample-specific FASTQ files in `READ_DIR`, matched with patterns like:

- `Sample1*_R1*.fastq.gz`
- `Sample1*_R2*.fastq.gz`

### 2. Fusion targets TSV
A tab-delimited file with header:

```text
sample	fusion	lbp	rbp
```

Example:

```text
sample	fusion	lbp	rbp
SAMPLE_001	EML4--ALK	chr2:42522694	chr2:29446394
SAMPLE_002	ETV6--NTRK3	chr12:12022900	chr15:88483900
```

### 3. Panel gene list
One gene symbol per line.

## Outputs

For each sample, the workflow creates:

- raw and trimmed FastQC reports
- trimmed FASTQ files
- sorted, deduplicated, and optional BQSR BAMs
- Manta / GRIDSS / DELLY / SvABA outputs
- AnnotSV annotated TSVs
- panel-filtered TSVs
- per-sample logs

A combined summary file is also created:

```text
WORKROOT/summary/summary.tsv
```

This includes:

- sample
- fusion
- expected breakpoints
- final BAM used
- mean depth around breakpoint windows
- supplementary alignment support
- caller-specific BND hit counts
- Manta run directory path

## Requirements

This pipeline assumes access to the following tools:

- conda
- fastqc
- fastp
- bwa
- samtools
- bcftools
- bedtools
- GATK
- Manta
- GRIDSS2
- DELLY
- SvABA
- AnnotSV
- bgzip
- tabix

## Setup

### 1. Clone the repository

```bash
git clone https://github.com/YOUR_USERNAME/fusion-sv-panel-pipeline.git
cd fusion-sv-panel-pipeline
```

### 2. Copy and edit the config

```bash
cp config/config.example.sh config/my_run_config.sh
nano config/my_run_config.sh
```

Edit all paths and environment names.

### 3. Run the pipeline

```bash
bash run_fusion_sv_pipeline.sh config/my_run_config.sh
```

## Notes on references

This version is written for hg19 / GRCh37-style resources. If you adapt it for GRCh38, update:

- reference FASTA
- index files
- dbSNP resource
- panel BED
- blacklist file
- AnnotSV database
- fusion target coordinates

## Suggested tested-version tracking

See `env/tool_versions.txt` and replace the placeholders with your actual validated versions before publishing or tagging a release.

## Recommended publication checklist

Before making the repository public, update:

- your name and affiliation if desired
- tool versions
- example input files
- LICENSE choice
- validation statement in README
- manuscript or preprint citation if applicable

## Suggested citation section

If you use this repository in a study, cite:

- the original software tools used in the workflow
- your assay validation or benchmark study
- this GitHub repository release tag

## License

This repository currently includes an MIT License template. Replace it if your institute or project requires a different license.
