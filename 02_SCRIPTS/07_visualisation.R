library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)

log2_data <- readRDS(
  "01_RAW_DATA/analysis_expression.rds"
)

metadata <- readRDS(
  "01_RAW_DATA/analysis_metadata.rds"
)
limma_results <- read.csv("01_RAW_DATA/limma_results.csv", row.names = 1)


if (!dir.exists("03_OUTPUTS/figures")) {
  dir.create("03_OUTPUTS/figures", recursive = TRUE)
}

limma_results$direction <- case_when(
  limma_results$logFC > 0.5 & limma_results$adj.P.Val < 0.05 ~ "upregulated",
  limma_results$logFC < -0.5 & limma_results$adj.P.Val < 0.05 ~ "downregulated",
  TRUE ~ "not significant"
)

print(table(limma_results$direction))

significant_genes <- limma_results %>%
  filter(abs(logFC) > 0.5 & adj.P.Val < 0.05)

cat("   Total significant:", nrow(significant_genes), "\n")
cat("   Upregulated:", sum(significant_genes$logFC > 0), "\n")
cat("   Downregulated:", sum(significant_genes$logFC < 0), "\n")

volcano_plot <- ggplot(limma_results, aes(x = logFC, y = -log10(adj.P.Val), color = direction)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  scale_color_manual(
    values = c("upregulated" = "red", 
               "downregulated" = "blue", 
               "not significant" = "gray")
  ) +
  theme_classic() +
  labs(
    title = "Volcano Plot: Liver Tumor vs Non-Tumor",
    x = "Log2 Fold Change",
    y = "-Log10(Adjusted P-Value)",
    color = "Expression"
  ) +
  theme(legend.position = "top")

print(volcano_plot)
ggsave("03_OUTPUTS/figures/volcano_plot.png", volcano_plot, width = 8, height = 6, dpi = 300)

top_genes <- limma_results %>%
  filter(abs(logFC) > 0.5 & adj.P.Val < 0.05) %>%
  arrange(adj.P.Val) %>%
  head(50)


heatmap_data <- log2_data[
  rownames(top_genes),
  ,
  drop = FALSE
]

annotation <- data.frame(
  Tissue = metadata$characteristics_ch1,
  row.names = rownames(metadata)
)

annotation_colors <- list(
  Tissue = c(
    "Tumor Tissue" = "red",
    "Non-Tumor Tissue" = "blue"
  )
)

if (nrow(heatmap_data) >= 2 && ncol(heatmap_data) >= 2) {
  
  pheatmap(
    heatmap_data,
    scale = "row",
    annotation_col = annotation,
    annotation_colors = annotation_colors,
    show_rownames = FALSE,
    show_colnames = FALSE,
    main = paste0(
      "Top ", nrow(top_genes),
      " Differentially Expressed Probe Sets"
    ),
    fontsize = 8,
    border_color = NA,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    filename = "03_OUTPUTS/figures/heatmap.png",
    width = 10,
    height = 8
  )
  
} else {
  
  message("Not enough significant probe sets or samples for clustering.")
  
}

 pca <- prcomp(
   t(log2_data),
   center = TRUE,
   scale. = FALSE
 )
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Sample = rownames(pca$x)
)

pca_df$Condition <- metadata[
  pca_df$Sample,
  "characteristics_ch1"
]

var_explained <- round(summary(pca)$importance[2, 1:2] * 100, 1)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Condition), level = 0.95) +
  theme_classic() +
  labs(
    title = "PCA Plot: Liver Tumor vs Non-Tumor",
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)"),
    color = "Tissue"
  ) +
  theme(legend.position = "top") +
  scale_color_manual(values = c("Tumor Tissue" = "red", "Non-Tumor Tissue" = "blue"))

print(pca_plot)
ggsave("03_OUTPUTS/figures/pca_plot.png", pca_plot, width = 8, height = 6, dpi = 300)

ma_plot <- ggplot(limma_results, aes(x = AveExpr, y = logFC, color = direction)) +
  geom_point(alpha = 0.5, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray50") +
  scale_color_manual(
    values = c("upregulated" = "red", 
               "downregulated" = "blue", 
               "not significant" = "gray")
  ) +
  theme_classic() +
  labs(
    title = "MA Plot: Liver Tumor vs Non-Tumor",
    x = "Average Expression (log2)",
    y = "Log2 Fold Change",
    color = "Expression"
  ) +
  theme(legend.position = "top")

print(ma_plot)
ggsave("03_OUTPUTS/figures/ma_plot.png", ma_plot, width = 8, height = 6, dpi = 300)

