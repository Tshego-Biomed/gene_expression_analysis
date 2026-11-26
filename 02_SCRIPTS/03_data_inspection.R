expression_data<- readRDS("01_RAW_DATA/expression_data.rds") 
metadata<- readRDS("01_RAW_DATA/metadata.rds") 

print(dim(expression_data))
print(dim(metadata)) 

print(
  quantile(
    expression_data,
    probs = c(0, 0.25, 0.5, 0.75, 0.99, 1),
    na.rm = TRUE
  )
)

head(expression_data[, 1:5])
head(metadata) 

colnames(metadata) 

hist(
  as.vector(expression_data),
  breaks = 50,
  main = "Downloaded expression distribution",
  xlab = "Log2 expression intensity"
)

metadata_columns <- grep(
  "title|source_name|characteristics",
  colnames(metadata),
  value = TRUE
)

print(metadata_columns)

print(
  head(
    metadata[, metadata_columns, drop = FALSE],
    12
  )
)

print(
  unique(metadata$characteristics_ch1)
)
