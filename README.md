# Gene Expression Analysis of Hepatocellular Carcinoma

## Project overview
This beginner R portfolio project compares expression in hepatocellular carcinoma tumour tissue and adjacent non-tumour liver tissue using the public GSE14520 dataset from NCBI GEO.
I used the project to practise downloading public data, inspecting sample metadata, cleaning expression matrices, identifying matched samples, fitting a differential-expression model and visualising results.

## Dataset
- Dataset: GSE14520.
- Selected platform: GPL3921, Affymetrix HT Human Genome U133A Array.
- Data type: microarray expression profiling.
- Input: processed, RMA-normalised log2 expression values downloaded through GEOquery.
- Final analysis: 213 matched patients, comprising 213 tumour and 213 adjacent non-tumour samples.
The downloaded expression matrix contains processed values, not raw sequencing counts or raw CEL-file measurements. No additional log transformation is applied.

## Sample selection and pairing
Samples are restricted to tumour and non-tumour tissue from the selected platform. Expression columns and metadata rows are aligned by sample ID.
Patient pairs are identified from the source labels by removing the terminal tissue suffix: A for tumour and B for non-tumour. Full identifier prefixes are preserved. Only identifiers with exactly one sample of each tissue type are retained.
One source-based pair, LCS-190, is excluded conservatively because its individual identifiers are inconsistent: 03-187A for the tumour sample and LCS-190B for the non-tumour sample. This is an unresolved naming discrepancy, not proof that the samples belong to different patients.
For the 213 retained pairs, the individual labels agree within each pair after removing the tissue suffix and surrounding whitespace. The exported patient mapping records the sample assignments.

## Analysis workflow
1. Load packages and create working folders.
2. Download GSE14520 and explicitly select GPL3921.
3. Inspect expression values and sample metadata.
4. Remove probe sets containing non-finite values or having zero variance across the downloaded samples.
5. Preserve the supplied normalised log2 expression values.
6. Select matched samples and fit a paired limma model.
7. Generate volcano, heatmap, PCA and MA plots.
8. Export summary statistics and the top significant upregulated and downregulated probe sets.
The filter is a basic data-quality filter; it does not specifically remove low-expression probe sets.
The differential-expression design is:
~ Patient_ID + Tissue
Patient identity accounts for the matching. Non-tumour tissue is the reference, and the TissueTumor coefficient estimates tumour versus non-tumour expression. limma's empirical Bayes procedure is applied, and p-values are adjusted using the Benjamini–Hochberg method.
Probe sets are classified as significant when both conditions are met:
- Absolute log2 fold change greater than 0.5.
- BH-adjusted p-value below 0.05.
Positive logFC indicates higher expression in tumour tissue; negative logFC indicates lower expression.

## Key findings
Measure	Result
Matched patients	213
Total samples	426
Tumour samples	213
Non-tumour samples	213
Probe sets tested	22,268
Significant probe sets	4,704
Upregulated in tumour	2,976
Downregulated in tumour	1,728


Among significant probe sets, 209220_at has the largest positive logFC, approximately 4.87. Probe set 220491_at has the most negative logFC, approximately -5.45.
These are probe-set identifiers, not gene names. Multiple probe sets may represent the same gene, so 4,704 significant probe sets does not mean 4,704 distinct genes.

## Visualisations
- Volcano plot: log2 fold change against negative log10 BH-adjusted p-value, coloured by the stated significance criteria.
- Heatmap: up to 50 significant probe sets ranked by adjusted p-value. Row scaling shows relative expression patterns within each probe set, not absolute expression differences between probe sets.
- PCA plot: an unsupervised view of the retained samples using centred log2 expression values without feature scaling. The group ellipses are descriptive and are not a statistical test of separation.
- MA plot: log2 fold change against average log2 expression.
PCA is calculated from the expression matrix, not from patient-adjusted model residuals. The heatmap uses features selected by the tissue comparison, so it is not independent validation of that comparison.

## Limitations
This is an exploratory reanalysis for learning. It relies on the submitted expression processing and sample annotations rather than reprocessing raw CEL files.
The model accounts for patient pairing but does not explicitly adjust for additional clinical variables or technical batch effects. Adjacent non-tumour tissue should not be interpreted as healthy donor liver tissue.
The analysis remains at probe-set level without gene annotation or independent validation. Differential expression describes associations in these samples and does not establish causation, clinical usefulness or validated biomarkers.

## Project files
Location	Contents
01_RAW_DATA/	Downloaded expression and metadata, plus generated intermediate files and full limma results
02_SCRIPTS/	Eight numbered R scripts
03_OUTPUTS/figures/	Four PNG figures
03_OUTPUTS/tables/	Summary, top-result and patient-mapping CSV files
Gene Expression Analysis.Rproj	RStudio project
README.md	Project description and instructions
Disclaimer.md	Learning background and acknowledgements
Although named 01_RAW_DATA, that folder also contains processed and intermediate data. The current .gitignore excludes it, generated CSV files and the output-table folder from Git. These files are regenerated by the workflow; the original dataset remains publicly accessible through GEO.

## How to run
Download or clone the project and open Gene Expression Analysis.Rproj in RStudio. The working directory should be the project folder containing 02_SCRIPTS.
Install the required packages once:
install.packages(c("tidyverse", "pheatmap", "RColorBrewer"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("GEOquery", "limma"))
Restart R and ensure that saved workspace objects have not been restored. Run the scripts in order:
``` r
source("02_SCRIPTS/01_load_packages.R")
source("02_SCRIPTS/02_data_download.R")
source("02_SCRIPTS/03_data_inspection.R")
source("02_SCRIPTS/04_data_cleaning.R")
source("02_SCRIPTS/05_normalization.R")
source("02_SCRIPTS/06_differential_expression.R")
source("02_SCRIPTS/07_visualisation.R")
source("02_SCRIPTS/08_statistics.R")
```
Internet access is required for package installation and downloading GEO data. Script 1 creates the required folders, and script 2 saves the downloaded expression matrix and metadata.
If a file cannot be found, check the working directory and whether the preceding script completed:
getwd()
file.exists("01_RAW_DATA/metadata.rds")
file.exists("01_RAW_DATA/log2_data.rds")

## Generated outputs
Figures:
- 03_OUTPUTS/figures/volcano_plot.png
- 03_OUTPUTS/figures/heatmap.png
- 03_OUTPUTS/figures/pca_plot.png
- 03_OUTPUTS/figures/ma_plot.png
Tables:
- 01_RAW_DATA/limma_results.csv
- 03_OUTPUTS/tables/summary_statistics.csv
- 03_OUTPUTS/tables/top_upregulated_genes.csv
- 03_OUTPUTS/tables/top_downregulated_genes.csv
- 03_OUTPUTS/tables/patient_mapping.csv
The top-result filenames retain the word “genes”, but their rows contain probe sets. The heatmap is generated only when enough significant probe sets and samples are available for clustering.

## Tools and skills practised
R, tidyverse (including ggplot2), GEOquery, limma, pheatmap and RColorBrewer.
Skills include public-data retrieval, matrix handling, metadata cleaning, sample alignment, paired modelling, multiple-testing correction, visualisation and documenting limitations.

## Author and acknowledgements
Tshegofatso Leope
BHSc Biomedical Sciences, University of the Witwatersrand
Initially completed on 10 July 2026; refined in October 2026.
The original study's researchers made the data available through NCBI GEO. Online courses, reference materials and AI assistance supported learning and troubleshooting. See Disclaimer.md for further acknowledgements.