library(tidyverse)
library(GEOquery)
library(limma)
library(pheatmap)
library(RColorBrewer)

dir.create("01_RAW_DATA", showWarnings = FALSE)

dir.create(
  "03_OUTPUTS/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "03_OUTPUTS/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

print("Packages loaded and folders prepared.")
