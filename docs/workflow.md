# Workflow summary

FASTQ
-> FastQC (raw)
-> fastp trimming
-> FastQC (trimmed)
-> BWA-MEM alignment
-> sort BAM
-> MarkDuplicates
-> optional BQSR
-> Manta
-> GRIDSS2
-> DELLY
-> SvABA
-> AnnotSV
-> panel-gene filter
-> breakpoint evidence summary

## Breakpoint evidence extracted

For non-NA left and right breakpoints:

- mean depth around left breakpoint window
- mean depth around right breakpoint window
- supplementary alignment bridging from left to right
- supplementary alignment bridging from right to left
- BND hits from:
  - Manta
  - GRIDSS2
  - DELLY
  - SvABA

## Summary output columns

- sample
- fusion
- lbp
- rbp
- final_bam
- meanDepth_L400
- meanDepth_R400
- SA_LtoR_500
- SA_RtoL_500
- manta_bnd_hits
- gridss_bnd_hits
- delly_bnd_hits
- svaba_bnd_hits
- manta_rundir
