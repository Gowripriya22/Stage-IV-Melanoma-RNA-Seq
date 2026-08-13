# ============================================================
# DESeq2 Differential Expression Analysis
# Stage IV Melanoma RNA-Seq Project
# ============================================================

# Load required libraries
library(DESeq2)
library(tidyverse)

# ------------------------------------------------------------
# 1. Read count data
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# 2. Define experimental conditions
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# 3. Create DESeq2 dataset
# ------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = colData,
  design = ~ condition
)

# ------------------------------------------------------------
# 4. Filter low-count genes
# ------------------------------------------------------------

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# ------------------------------------------------------------
# 5. Run DESeq2
# ------------------------------------------------------------

dds <- DESeq(dds)

# ------------------------------------------------------------
# 6. Extract differential expression results
# ------------------------------------------------------------

res <- results(
  dds,
  contrast = c("condition", "control", "disease")
)

# Order results by adjusted p-value
res_ordered <- res[order(res$padj), ]

# Save complete DESeq2 results
write.csv(
  as.data.frame(res_ordered),
  "DESeq2_results.csv"
)

# ------------------------------------------------------------
# 7. Identify significant genes
# ------------------------------------------------------------

significant_genes <- subset(
  res_ordered,
  padj < 0.05
)

write.csv(
  as.data.frame(significant_genes),
  "Significant_Genes.csv"
)

# Display DESeq2 summary
summary(res)

# ------------------------------------------------------------
# 8. MA Plot
# ------------------------------------------------------------

pdf("MA_plot.pdf")

plotMA(
  res,
  main = "MA Plot of Differential Expression",
  ylim = c(-5, 5)
)

dev.off()

# ------------------------------------------------------------
# 9. PCA Plot
# ------------------------------------------------------------

rld <- rlog(dds)

pcaData <- plotPCA(
  rld,
  intgroup = "condition",
  returnData = TRUE
)

percentVar <- round(
  100 * attr(pcaData, "percentVar")
)

pca_plot <- ggplot(
  pcaData,
  aes(PC1, PC2, color = condition)
) +
  geom_point(size = 3) +
  xlab(
    paste0(
      "PC1: ",
      percentVar[1],
      "% variance"
    )
  ) +
  ylab(
    paste0(
      "PC2: ",
      percentVar[2],
      "% variance"
    )
  ) +
  ggtitle("PCA Plot of Samples") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  "PCA_plot.png",
  plot = pca_plot,
  width = 6,
  height = 4,
  dpi = 300
)

print(pca_plot)

# ------------------------------------------------------------
# 10. Volcano Plot
# ------------------------------------------------------------

volcano_data <- as.data.frame(res_ordered)

volcano_data$significance <- ifelse(
  !is.na(volcano_data$padj) &
  volcano_data$padj < 0.05 &
  abs(volcano_data$log2FoldChange) > 1,
  "Significant",
  "Not Significant"
)

volcano_plot <- ggplot(
  volcano_data,
  aes(
    x = log2FoldChange,
    y = -log10(pvalue),
    color = significance
  )
) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10(p-value)",
    title = "Volcano Plot"
  ) +
  theme(
    legend.position = "top"
  )

ggsave(
  "Volcano_plot.png",
  plot = volcano_plot,
  width = 7,
  height = 5,
  dpi = 300
)

print(volcano_plot)
