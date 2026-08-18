## Module 2: Data Acquisition & Preprocessing

In this module, we bridge the gap between raw sequencing outputs and analytical readiness. Firstly, we will fetch the real clinical cancer patient data from The Cancer Genome Atlas (TCGA). The format it lands in will be somewhat preprocessed, which often times is what you as a wet lab researcher might begin with, or will be made available online by other researchers. We can discuss raw data reads at another stage.  

When dealing with transcriptomic data, such as RNA sequencing (RNA-Seq), you will often begin with a large matrix of raw read counts. This is the output of processing raw FASTQ files (which requires QC steps and is usually done outside of R), so for now we will start with a counts matrix. These counts simply represent the absolute number of sequencing reads that successfully map to a given gene in the genome. However, these raw numbers cannot be compared directly across different samples. Variations in sequencing depth (the total number of reads generated for a specific sample) and overall RNA composition dictate that the data must be mathematically adjusted, a process known as normalisation. Furthermore, thousands of genes are either not expressed or expressed at exceptionally low levels in our tissue of interest. Retaining these lowly expressed genes introduces significant statistical noise and ultimately weakens the power of our subsequent analyses.  

In the `src02_preprocessing.R` script below, we will programmatically query and download a real transcriptomic dataset from The Cancer Genome Atlas (TCGA). Once the raw count matrix is assembled, we will filter out uninformative low-count genes and apply Trimmed Mean of M-values (TMM) normalisation to correct for technical sequencing biases. Finally, we will apply a `voom` transformation. This crucial step converts our raw counts into log2 values while assigning mathematical weights to each observation based on its variance, completely preparing the dataset for robust linear modelling and differential expression testing in Module 3.  

```r
# ==============================================================================
# Script: scripts/src02_preprocessing.R
# Purpose: Fetch real TCGA RNA-Seq data, process counts, filter, and normalise
# ==============================================================================

# We load our project libraries quietly to keep our console clean.
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(tidyr)
  library(edgeR)
  library(limma)
})

# ------------------------------------------------------------------------------
# 1. Configuration & Cohort Definition
# ------------------------------------------------------------------------------

# We define the project code and clean label here so we can easily swap cohort 
# identifiers (e.g., "TCGA-BRCA", "TCGA-LUAD", "TCGA-COAD") without rewriting code.
cancer_project    <- "TCGA-BRCA"
cancer_label      <- "Breast Cancer (BRCA)"
n_samples_per_grp <- 20 # Number of samples per condition to keep download sizes fast

message(sprintf("Initiating data query for project: %s", cancer_project))

# ------------------------------------------------------------------------------
# 2. Query GDC Portal & Download Data
# ------------------------------------------------------------------------------

# We query the Genomic Data Commons (GDC) portal for primary solid tumour 
# and solid tissue normal gene expression quantification files.
query_tcga <- GDCquery(
  project = cancer_project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

# We extract the metadata table from our query to subset a balanced sub-cohort.
results_df <- getResults(query_tcga)

# We select sample barcodes for our specified sample size per group.
tumor_barcodes <- results_df$cases[results_df$sample_type == "Primary Tumor"][1:n_samples_per_grp]
normal_barcodes <- results_df$cases[results_df$sample_type == "Solid Tissue Normal"][1:n_samples_per_grp]

# We build a focused query targeting only our selected sample barcodes.
query_sub <- GDCquery(
  project = cancer_project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = c(tumor_barcodes, normal_barcodes)
)

# We download the file payload directly via the GDC API.
GDCdownload(query_sub, method = "api", files.per.chunk = 10)

# We parse the downloaded files into a SummarizedExperiment container.
tcga_se <- GDCprepare(query_sub)

# ------------------------------------------------------------------------------
# 3. Extract Assay Data & Format Metadata
# ------------------------------------------------------------------------------

# We extract the raw unstranded count matrix from our SummarizedExperiment object.
counts_mat <- assay(tcga_se, "unstranded")

# We extract gene annotation metadata (e.g., symbol, gene type).
gene_info <- as.data.frame(rowData(tcga_se))

# We extract clinical/sample metadata and create a clean experimental condition factor, 
# explicitly setting "Normal" as our baseline level.
sample_metadata <- as.data.frame(colData(tcga_se)) %>%
  transmute(
    barcode = barcode,
    sample_type = sample_type,
    Condition = factor(
      ifelse(definition == "Primary solid Tumor", "Tumour", "Normal"),
      levels = c("Normal", "Tumour")
    )
  )

# ------------------------------------------------------------------------------
# 4. Low-Count Filtering & Voom Normalisation
# ------------------------------------------------------------------------------

# We construct a DGEList object to structure our counts and experimental groups for edgeR.
dge <- DGEList(counts = counts_mat, group = sample_metadata$Condition)

# We apply edgeR's filterByExpr heuristic to automatically discard uninformative, 
# low-count genes based on library sizes and group proportions.
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# We compute Trimmed Mean of M-values (TMM) scale factors to correct for composition bias.
dge <- calcNormFactors(dge)

# We transform our count data into log2 counts-per-million (CPM) using voom, 
# which models the mean-variance relationship to generate observational weights for limma.
v <- voom(dge, plot = FALSE)
log_expr <- v$E

# ------------------------------------------------------------------------------
# 5. Save Processed Artifacts
# ------------------------------------------------------------------------------

# We save our cleaned metadata, gene information, and normalised log2 expression matrix 
# to the data directory for downstream modular analyses.
write.csv(sample_metadata, "data/sample_metadata.csv", row.names = FALSE)
write.csv(gene_info, "data/gene_info.csv", row.names = TRUE)
saveRDS(log_expr, "data/normalised_log_expression.rds")
saveRDS(v, "data/voom_object.rds")

message("Preprocessing complete! Cleaned data saved to data/")
```  

---
