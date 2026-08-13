# Load packages

library(DESeq2)
library(tidyverse)

# Read count data

counts_data <- read.csv(
  "COUNTS.csv",
  header = TRUE,
  row.names = "Geneid"
)

# Select sample count columns

counts <- counts_data[, c(
  "control_rep1",
  "control_rep2",
  "control_rep3",
  "disease_rep1",
  "disease_rep2",
  "disease_rep3"
)]

# Define experimental conditions

condition <- factor(c(
  "control",
  "control",
  "control",
  "disease",
  "disease",
  "disease"
))

# Create sample metadata

colData <- data.frame(
  row.names = colnames(counts),
  condition = condition
)

# Create DESeq2 dataset

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = colData,
  design = ~ condition
)

# Filter low-count genes

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Run DESeq2

dds <- DESeq(dds)

# Extract differential expression results

res <- results(
  dds,
  contrast = c("condition", "control", "disease")
)

# Order by adjusted p-value

res_ordered <- res[order(res$padj), ]

# Save complete results

write.csv(
  as.data.frame(res_ordered),
  "DESeq2_results.csv"
)

# Identify significant genes

significant_genes <- subset(
  res_ordered,
  padj < 0.05
)

# Save significant genes

write.csv(
  as.data.frame(significant_genes),
  "Significant_Genes.csv"
)

# Display summary

summary(res)
