expression_data <- readRDS(
  "01_RAW_DATA/expression_data.rds"
)

finite_rows <- rowSums(
  !is.finite(expression_data)
) == 0

finite_data <- expression_data[
  finite_rows,
  ,
  drop = FALSE
]

probe_variance <- apply(
  finite_data,
  1,
  var
)

filtered_data <- finite_data[
  is.finite(probe_variance) & probe_variance > 0,
  ,
  drop = FALSE
]

cat(
  "Original probe sets:",
  nrow(expression_data),
  "\n"
)

cat(
  "Retained probe sets:",
  nrow(filtered_data),
  "\n"
)

cat(
  "Removed probe sets:",
  nrow(expression_data) - nrow(filtered_data),
  "\n"
)

saveRDS(
  filtered_data,
  "01_RAW_DATA/filtered_data.rds"
)

write.csv(
  filtered_data,
  "01_RAW_DATA/filtered_data.csv"
)
