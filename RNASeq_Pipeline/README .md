# RNA-Seq Analysis Pipeline

This folder documents the command-line workflow used for bulk RNA-Seq
data analysis.

The workflow covers the major stages from raw FASTQ quality assessment
to gene-level quantification.

## Workflow

Raw FASTQ
↓
FastQC
↓
Read trimming with fastp
↓
FastQC after trimming
↓
MultiQC
↓
Reference genome indexing
↓
HISAT2 alignment
↓
SAM to BAM conversion and sorting
↓
featureCounts gene-level quantification
↓
DESeq2 differential expression analysis
↓
Functional interpretation

## Tools Used

| Stage | Tool |
|---|---|
| Initial QC | FastQC |
| Read trimming | fastp |
| QC aggregation | MultiQC |
| Reference indexing | HISAT2 |
| Read alignment | HISAT2 |
| SAM/BAM processing | SAMtools |
| Gene quantification | featureCounts |
| Differential expression | DESeq2 |
| Statistical analysis | R |
| Visualization | ggplot2 |

## Reference Genome

Human GRCh38 reference genome.

## Important Note

The commands documented here represent the workflow used during
hands-on RNA-Seq analysis.

Large raw FASTQ files, reference genomes and BAM files are not uploaded
to GitHub because of file-size and storage considerations.

The purpose of this repository is to demonstrate the analysis workflow,
command-line skills, data processing steps and interpretation of RNA-Seq
results.
