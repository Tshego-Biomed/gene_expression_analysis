library(dplyr)

log2_data <- readRDS(
  "01_RAW_DATA/analysis_expression.rds"
)

metadata <- readRDS(
  "01_RAW_DATA/analysis_metadata.rds"
)
limma_results <- read.csv("01_RAW_DATA/limma_results.csv", row.names = 1)

metadata$characteristics_ch1 <- gsub("tissue: ", "", metadata$characteristics_ch1)
metadata$characteristics_ch1 <- gsub("Tissue: ", "", metadata$characteristics_ch1)
metadata$characteristics_ch1 <- gsub("Liver ", "", metadata$characteristics_ch1)
metadata$characteristics_ch1 <- trimws(metadata$characteristics_ch1)

metadata <- metadata[metadata$characteristics_ch1 != "", ]
metadata <- metadata[!is.na(metadata$characteristics_ch1), ]


print(head(metadata$characteristics_ch1, 20))
print(unique(metadata$characteristics_ch1))

metadata$characteristics_ch1[metadata$characteristics_ch1 == "Tumor Tissue"] <- "Tumor Tissue"
metadata$characteristics_ch1[metadata$characteristics_ch1 == "Non-Tumor Tissue"] <- "Non-Tumor Tissue"

metadata <- metadata[metadata$characteristics_ch1 %in% c("Tumor Tissue", "Non-Tumor Tissue"), ]

print(unique(metadata$characteristics_ch1))
print(table(metadata$characteristics_ch1))

if (length(unique(metadata$characteristics_ch1)) < 2) {
  stop("❌ Only one group found! Check metadata column.")
}

common_samples <- intersect(colnames(log2_data), rownames(metadata))
log2_data <- log2_data[, common_samples]
metadata <- metadata[common_samples, ]

cat("Samples after alignment:", ncol(log2_data), "\n")

significant_genes <- limma_results %>%
  filter(abs(logFC) > 0.5 & adj.P.Val < 0.05)

cat("Significant Probe Sets:", nrow(significant_genes), "\n")

summary_stats <- data.frame(
  Metric = c(
    "Total Samples",
    "Tumor Samples",
    "Non-Tumor Samples",
    "Total Probe Sets Tested",
    "Significant Probe Sets (|logFC| > 0.5, adj.P.Val < 0.05)",
    "Upregulated in Tumor",
    "Downregulated in Tumor",
    "Most Upregulated Probe Set",
    "Most Upregulated Probe Set LogFC",
    "Most Downregulated Probe Set",
    "Most Downregulated Probe Set LogFC"
  ),
  Value = c(
    ncol(log2_data),
    sum(metadata$characteristics_ch1 == "Tumor Tissue"),
    sum(metadata$characteristics_ch1 == "Non-Tumor Tissue"),
    nrow(limma_results),
    nrow(significant_genes),
    sum(significant_genes$logFC > 0),
    sum(significant_genes$logFC < 0),
    if(nrow(significant_genes) > 0) rownames(significant_genes)[which.max(significant_genes$logFC)] else "N/A",
    if(nrow(significant_genes) > 0) round(max(significant_genes$logFC), 2) else "N/A",
    if(nrow(significant_genes) > 0) rownames(significant_genes)[which.min(significant_genes$logFC)] else "N/A",
    if(nrow(significant_genes) > 0) round(min(significant_genes$logFC), 2) else "N/A"
  )
)

print(summary_stats)

top_upregulated <- significant_genes %>%
  filter(logFC > 0) %>%
  arrange(desc(logFC)) %>%
  head(10)


print(top_upregulated[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")])


top_downregulated <- significant_genes %>%
  filter(logFC < 0) %>%
  arrange(logFC) %>%
  head(10)

print(top_downregulated[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")])

if(!dir.exists("03_OUTPUTS/tables")) {
  dir.create("03_OUTPUTS/tables", recursive = TRUE)
}

write.csv(summary_stats, "03_OUTPUTS/tables/summary_statistics.csv", row.names = FALSE)

write.csv(top_upregulated, "03_OUTPUTS/tables/top_upregulated_genes.csv", row.names = TRUE)
write.csv(top_downregulated, "03_OUTPUTS/tables/top_downregulated_genes.csv", row.names = TRUE)


cat("\n Results Summary:\n")
cat("   - Total samples:", ncol(log2_data), "\n")
cat("   - Tumor samples:", sum(metadata$characteristics_ch1 == "Tumor Tissue"), "\n")
cat("   - Non-Tumour samples:", sum(metadata$characteristics_ch1 == "Non-Tumor Tissue"), "\n")
cat("   - Total probe sets tested:", nrow(limma_results), "\n")
cat("   - Significant Probe Sets (|logFC| > 0.5, adj.P.Val < 0.05)", nrow(significant_genes), "\n")
cat("   - Upregulated in Tumor:", sum(significant_genes$logFC > 0), "\n")
cat("   - Downregulated in Tumor:", sum(significant_genes$logFC < 0), "\n")
cat("   - Most upregulated:", summary_stats$Value[8], "(", summary_stats$Value[9], ")\n")
cat("   - Most downregulated:", summary_stats$Value[10], "(", summary_stats$Value[11], ")\n") 

