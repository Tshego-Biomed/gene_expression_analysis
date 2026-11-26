filtered_data <- readRDS(
  "01_RAW_DATA/filtered_data.rds"
)

log2_data <- filtered_data

saveRDS(
  log2_data,
  "01_RAW_DATA/log2_data.rds"
)
