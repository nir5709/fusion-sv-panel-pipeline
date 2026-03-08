# Fusion SV Panel Pipeline

![GitHub repo size](https://img.shields.io/github/repo-size/nir5709/fusion-sv-panel-pipeline)
![GitHub last commit](https://img.shields.io/github/last-commit/nir5709/fusion-sv-panel-pipeline)
![GitHub license](https://img.shields.io/github/license/nir5709/fusion-sv-panel-pipeline)
![Platform](https://img.shields.io/badge/platform-Linux-blue)
![Genome](https://img.shields.io/badge/genome-hg19%20%2F%20GRCh37-green)

A **targeted fusion and structural variant detection pipeline** for paired-end sequencing data using multiple SV callers and panel-focused annotation.

---

# Overview

This workflow processes paired-end FASTQ files and performs:

1. Raw read quality control using **FastQC**
2. Adapter trimming and quality filtering using **fastp**
3. Quality control of trimmed reads
4. Alignment using **BWA-MEM**
5. BAM sorting and indexing using **SAMtools**
6. Duplicate marking using **GATK MarkDuplicates**
7. Optional **Base Quality Score Recalibration (BQSR)**
8. Structural variant detection using:

* Manta
* GRIDSS2
* DELLY
* SvABA

9. Structural variant annotation using **AnnotSV**
10. Panel gene filtering
11. Breakpoint evidence evaluation
12. Final summary report generation

---

## Pipeline Workflow

```
FASTQ
↓
FastQC (raw reads)
↓
fastp trimming
↓
FastQC (trimmed reads)
↓
BWA-MEM alignment
↓
SAMtools sorting
↓
GATK MarkDuplicates
↓
(Optional) BQSR
↓
SV Detection
  ├── Manta
  ├── GRIDSS2
  ├── DELLY
  └── SvABA
↓
AnnotSV Annotation
↓
Panel Gene Filtering
↓
Breakpoint Evidence Analysis
↓
Summary Report
```

---

# Why Multiple SV Callers?

Structural variant detection algorithms use different signals from sequencing data.

| Tool    | Detection Strategy          |
| ------- | --------------------------- |
| Manta   | paired-end + split-read     |
| GRIDSS2 | assembly-based SV detection |
| DELLY   | paired-end + read-depth     |
| SvABA   | local assembly              |

Using multiple callers improves:

* sensitivity
* breakpoint resolution
* detection of complex rearrangements

Consensus evidence across callers increases confidence in detected fusion events.

---

# Methodological Approach

This workflow integrates multiple structural variant detection algorithms and evaluates breakpoint evidence using several complementary signals:

- split-read alignments
- discordant paired-end reads
- local assembly-based variant detection
- supplementary alignment bridging
- depth-of-coverage near breakpoints

Combining these signals improves sensitivity and breakpoint resolution for fusion detection in targeted sequencing panels.

---

# Key Features

* Multi-caller SV detection strategy
* Breakpoint-level evidence evaluation
* Panel-specific gene filtering using AnnotSV
* Structured per-sample output directories
* Automated summary table across samples
* Config-driven setup (no hardcoded machine paths)

---

# Repository Structure

```
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

---

# Input Requirements

## 1. Paired-end FASTQ files

The pipeline expects paired FASTQ files in the input directory (`READ_DIR`):

```
Sample1*_R1*.fastq.gz
Sample1*_R2*.fastq.gz
```

---

## 2. Fusion Target File

Example format:

```
sample	fusion	lbp	rbp
SAMPLE_001	EML4--ALK	chr2:42522694	chr2:29446394
SAMPLE_002	ETV6--NTRK3	chr12:12022900	chr15:88483900
```

---

## 3. Panel Gene List

One gene symbol per line:

```
ALK
RET
ROS1
NTRK1
NTRK2
NTRK3
```

---

# Output

For each sample the pipeline produces:

### Quality Control

* Raw FastQC reports
* Trimmed FastQC reports

### Processed Reads

* Trimmed FASTQ files
* Sorted BAM
* Deduplicated BAM
* Optional BQSR BAM

### Structural Variant Results

* Manta calls
* GRIDSS2 calls
* DELLY calls
* SvABA calls

### Annotation

* AnnotSV annotated variants
* Panel-filtered variant tables

### Logs

* Per-sample log files

---

# Summary Output

A combined summary table is produced:

```
WORKROOT/summary/summary.tsv
```

Key fields include:

* sample
* fusion
* breakpoint coordinates
* mean depth near breakpoints
* supplementary alignment support
* SV caller evidence

---

# Software Requirements

* Conda
* FastQC
* fastp
* BWA
* SAMtools
* BCFtools
* BEDtools
* GATK
* Manta
* GRIDSS2
* DELLY
* SvABA
* AnnotSV
* bgzip
* tabix

The pipeline was developed and tested using **latest stable versions of all tools**.

---

# Reference Genome

```
Human genome reference: hg19 / GRCh37
```

If adapting to **GRCh38**, update:

* reference genome
* dbSNP resource
* panel BED file
* AnnotSV database
* blacklist files
* fusion breakpoint coordinates

---

# Installation

Clone repository:

```
git clone https://github.com/nir5709/fusion-sv-panel-pipeline.git
cd fusion-sv-panel-pipeline
```

---

# Configuration

Create config file:

```
cp config/config.example.sh config/my_run_config.sh
nano config/my_run_config.sh
```

Edit all required paths.

---

# Running the Pipeline

```
bash run_fusion_sv_pipeline.sh config/my_run_config.sh
```

---

# Use Case

This workflow is suitable for:

* targeted oncology panels
* fusion detection
* structural variant analysis
* research sequencing studies

---

# Disclaimer

This pipeline is intended for **research use and assay development**.

Clinical use requires:

* local validation
* QC thresholds
* regulatory compliance
* standardized reporting workflows.

---

# Author

**Nihar Garg**  
Bioinformatics | Cancer Genomics

---

# Citation

If you use this pipeline in research, please cite:

Nihar Garg.  
Fusion SV Panel Pipeline.  
GitHub repository.  
https://github.com/nir5709/fusion-sv-panel-pipeline

---

# License

MIT License
