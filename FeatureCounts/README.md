# FeatureCounts Quantification

This folder contains the gene-level read count table generated during the RNA-Seq analysis using FeatureCounts.

## Purpose

FeatureCounts was used to assign aligned RNA-Seq reads to annotated genomic features and generate a count matrix for downstream differential expression analysis.

## Output

The `FeatureCounts_results.xlsx` file contains:

- Gene identifiers
- Chromosome
- Genomic start and end positions
- Strand information
- Feature length
- Read counts for the control samples
- Read counts for the disease samples

## Sample Groups

The count table contains samples from two groups:

- Control
- Disease

The resulting count matrix was used as input for the downstream DESeq2 differential expression analysis.

## File

| File | Description |
|---|---|
| `FeatureCounts_results.xlsx` | Gene-level read count table generated from aligned RNA-Seq reads |
