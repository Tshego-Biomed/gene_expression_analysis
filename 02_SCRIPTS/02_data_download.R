library(GEOquery) 

gse_id <- "GSE14520"

gset <- getGEO(
  gse_id,
  GSEMatrix = TRUE
)

platform_ids <- vapply(
  gset,
  annotation,
  character(1)
)

print(platform_ids)

# Analyse the HT Human Genome U133A platform.
selected_platform <- "GPL3921"

platform_index <- which(
  platform_ids == selected_platform
)

if (length(platform_index) != 1) {
  stop("The selected platform was not uniquely identified.")
}

gse_data <- gset[[platform_index]]

expression_data <- exprs(gse_data)
metadata <- pData(gse_data)

stopifnot(
  identical(
    colnames(expression_data),
    rownames(metadata)
  )
)

cat("Platform:", annotation(gse_data), "\n")
cat("Probe sets:", nrow(expression_data), "\n")
cat("Samples:", ncol(expression_data), "\n")

saveRDS(
  expression_data,
  "01_RAW_DATA/expression_data.rds"
)

saveRDS(
  metadata,
  "01_RAW_DATA/metadata.rds"
)
