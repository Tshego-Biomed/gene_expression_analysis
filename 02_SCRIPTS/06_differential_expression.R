library(limma) 

library(dplyr)

log2_data <- readRDS("01_RAW_DATA/log2_data.rds") 
metadata <- readRDS("01_RAW_DATA/metadata.rds")

cat("\n RAW values (first 20):\n")
print(head(metadata$characteristics_ch1, 20))

cat("\n ALL UNIQUE VALUES:\n")
print(unique(metadata$characteristics_ch1))


metadata$characteristics_ch1 <- gsub("Tissue: ", "", metadata$characteristics_ch1)
metadata$characteristics_ch1 <- gsub("tissue: ", "", metadata$characteristics_ch1)
metadata$characteristics_ch1 <- gsub("Liver ", "", metadata$characteristics_ch1)
metadata$characteristics_ch1 <- trimws(metadata$characteristics_ch1)


metadata <- metadata[metadata$characteristics_ch1 != "", ]
metadata <- metadata[!is.na(metadata$characteristics_ch1), ]


cat("\n Cleaned values (first 20):\n")
print(head(metadata$characteristics_ch1, 20))

cat("\n ALL UNIQUE VALUES AFTER CLEANING:\n")
print(unique(metadata$characteristics_ch1))


metadata$characteristics_ch1[metadata$characteristics_ch1 == "Tumor Tissue"] <- "Tumor Tissue"
metadata$characteristics_ch1[metadata$characteristics_ch1 == "Non-Tumor Tissue"] <- "Non-Tumor Tissue"


metadata <- metadata[metadata$characteristics_ch1 %in% c("Tumor Tissue", "Non-Tumor Tissue"), ]
common_samples <- colnames(log2_data)[
  colnames(log2_data) %in% rownames(metadata)
]

log2_data <- log2_data[
  ,
  common_samples,
  drop = FALSE
]

metadata <- metadata[
  common_samples,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(log2_data),
    rownames(metadata)
  )
)

print(
  table(metadata$characteristics_ch1)
)

cat("\n FINAL GROUPS:\n")
print(unique(metadata$characteristics_ch1))
cat("\n GROUP COUNTS:\n")
print(table(metadata$characteristics_ch1))

if (length(unique(metadata$characteristics_ch1)) < 2) {
  stop("Only one group found! Check metadata column.")
}


# Use the sample source labels to identify patient pairs.
sample_labels <- trimws(
  as.character(metadata$source_name_ch1)
)

# Check that every label ends with a tissue suffix.
unexpected <- is.na(sample_labels) |
  !grepl("^.+[AB]$", sample_labels)

if (any(unexpected)) {
  print(
    metadata[
      unexpected,
      c("title", "source_name_ch1", "characteristics_ch1.2"),
      drop = FALSE
    ]
  )
  
  stop("Some source labels need inspection.")
}

# Verify that A means tumour and B means non-tumour.
suffix <- substr(
  sample_labels,
  nchar(sample_labels),
  nchar(sample_labels)
)

expected_suffix <- ifelse(
  metadata$characteristics_ch1 == "Tumor Tissue",
  "A",
  "B"
)

if (any(suffix != expected_suffix)) {
  stop("Source-label suffixes do not match the tissue groups.")
}

# Remove only the final A or B; preserve the full prefix.
metadata$Patient_ID <- sub(
  "[AB]$",
  "",
  sample_labels
)
# Count samples for each patient and tissue.
pair_counts <- with(
  metadata,
  table(Patient_ID, characteristics_ch1)
)

# Keep patients with exactly one sample of each tissue.
complete_patients <- rownames(pair_counts)[
  pair_counts[, "Tumor Tissue"] == 1 &
    pair_counts[, "Non-Tumor Tissue"] == 1
]
# Exclude one pair with unresolved individual identifiers.
complete_patients <- setdiff(
  complete_patients,
  "LCS-190"
)
metadata <- metadata[
  metadata$Patient_ID %in% complete_patients,
  ,
  drop = FALSE
]

# Match expression columns to the retained metadata.
log2_data <- log2_data[
  ,
  rownames(metadata),
  drop = FALSE
]

stopifnot(
  identical(colnames(log2_data), rownames(metadata))
)

# Set non-tumour tissue as the reference.
metadata$Tissue <- factor(
  metadata$characteristics_ch1,
  levels = c("Non-Tumor Tissue", "Tumor Tissue"),
  labels = c("NonTumor", "Tumor")
)

metadata$Patient_ID <- factor(metadata$Patient_ID)

cat("Complete patient pairs:", nlevels(metadata$Patient_ID), "\n")
cat("Samples analysed:", ncol(log2_data), "\n")
print(table(metadata$Tissue))

# Account for patient identity and estimate the tissue effect.
design_matrix <- model.matrix(
  ~ Patient_ID + Tissue,
  data = metadata
)

if (qr(design_matrix)$rank < ncol(design_matrix)) {
  stop("The design matrix is not full rank.")
}

fit <- lmFit(log2_data, design_matrix)
fit <- eBayes(fit)

limma_results <- topTable(
  fit,
  coef = "TissueTumor",
  number = Inf,
  adjust.method = "BH"
)

significant_genes <- limma_results %>%
  filter(
    abs(logFC) > 0.5,
    adj.P.Val < 0.05
  )

write.csv(
  limma_results,
  "01_RAW_DATA/limma_results.csv"
)

saveRDS(
  log2_data,
  "01_RAW_DATA/analysis_expression.rds"
)

saveRDS(
  metadata,
  "01_RAW_DATA/analysis_metadata.rds"
)

cat(
  "Significant probe sets:",
  nrow(significant_genes),
  "\n"
)

write.csv(
  metadata[, c(
    "source_name_ch1",
    "characteristics_ch1.2",
    "characteristics_ch1",
    "Patient_ID"
  )],
  "03_OUTPUTS/tables/patient_mapping.csv",
  row.names = TRUE
)
