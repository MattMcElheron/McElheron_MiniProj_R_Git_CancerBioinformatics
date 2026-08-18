# Mini-Project: Introduction to R, Git, and Basic Cancer Transcriptomic Analyses  

[![Cohort](https://img.shields.io/badge/Cohort-TCGA-blue.svg)](https://portal.gdc.cancer.gov/)
[![Cohort](https://img.shields.io/badge/Cohort-TCGA%20(cBioPortal)-008080.svg)](https://www.cbioportal.org/)

Welcome to your introductory dry-lab mini-project. This hands-on tutorial is designed to build foundational skills in **computational biology**, **version control**, and **data science using R**, designed for someone with a wet lab background interested in adding some computational skills to their week-to-week as a researcher.  


### Learning Objectives
By completing this project, you will learn how to:
1. Set up a reproducible project structure using **RStudio** and **Git/GitHub**.
2. Query and download real cancer transcriptomic data directly from the GDC portal.
3. Perform data filtering, TMM normalisation, and variance-stabilising transformations (`voom`).
4. Conduct hypothesis testing (`limma`), multiple testing correction (FDR), and dimensionality reduction (**PCA**).
5. Generate publication-ready visualisations (**PCA plots** & **Volcano plots**) using `ggplot2`.

---  

## Software & Setup

Please complete the following installations before starting Module 1:

* [ ] Install **R (v4.3+)**: [cran.r-project.org](https://cran.r-project.org/)
* [ ] Install **RStudio Desktop**: [posit.co/download/rstudio-desktop](https://posit.co/download/rstudio-desktop/)
* [ ] Install **Git**: [git-scm.com](https://git-scm.com/)
* [ ] Create a **GitHub** account: [github.com](https://github.com/)

Installing R & RStudio should be relatively straightforward, although it can be a little slow depending on your machine. RStudio is our graphical interface which runs R, and R is a programming language and software environment designed for statistical computing, data analysis, and graphical visualisation. It is really commonly used in Bioinformatics, similar to Python. In my opinion, R is better for Viz and standard statistics, not as good for bigger scale machine learning projects, but overall much more forgiving for someone first experiencing coding.  


GitHub is a cloud-based platform that hosts software development projects and allows us as researchers or code developers to store, manage, track, and share code using Git, a distributed version control system. Think of it like a private or public Sharepoint/Onedrive where we can share our code documents with lab/team members or the wider public. In R, we often install and use "Packages" of code - GitHub is where one might host a package, or something as simple as a single chunk of code that we are working on. "Git" specifically is a program we use to send/pull our files to/from GitHub.  

---  

### Recommended Preparation: Software Carpentry
If you are completely new to programming, I strongly recommend completing the free, self-paced **Software Carpentry: Programming with R** course before or alongside this project. It provides a beginner-friendly overview of basic R syntax, data structures, and best practices using a clinical example:

* **Course Link:** [Software Carpentry – Programming with R](https://swcarpentry.github.io/r-novice-inflammation/)

---  

## Module 1: Project Setup & Version Control

In this module, we will lay the foundation for our dry-lab environment. Transitioning from the wet lab to computational biology can feel daunting, but the core principles of good science remain exactly the same. Just as we would never throw all our reagents, buffers, and samples into a single unlabelled box in the fridge, we should never throw all our data, code, and results into a single messy folder on our computer.  

By organising our project into standard subdirectories (folders) right from the outset, we ensure our work remains structured and logical. We want our code to work flawlessly on our machine today, next month, and on someone else's computer (a colleague, a reviewer, or a client) a year from now. Similar to maintaining a meticulous physical lab notebook with a very specific layout and clear protocols, we are creating a digital lab notebook that is actually readable for others. Think of recipe books: they usually follow a predictable structure. We do not want to be on the final step of baking a cake only to realise we are missing a vital ingredient because our notes were disorganised. Coding follows this exact same logic.  

Furthermore, we will establish "version control" using Git and GitHub. If you have ever saved manuscript drafts as `thesis_v1`, `thesis_final`, and `thesis_final_FINAL2`, you already understand the need for version control. Imageine your supervisor gave you edits on a paper without tracking revisions? Git acts as a digital time machine for our code, silently tracking every single change, deletion, or addition we make. If we accidentally break our script, we can easily revert to a previous working version. GitHub is simply the cloud-based platform where we store this time machine. It acts like a highly specialised Google Drive, allowing our updates to be backed up, tracked, and made fully reproducible for the wider scientific community.  

### Step 1: Create local directory structure
Open **RStudio**, create a new R Project linked to your cloned GitHub repository, and run the following in your Console:

```r
# Create standard project folders
dir.create("data")
dir.create("scripts")
dir.create("output")
```

Now, when someone else opens our script, they will know how our data/input/output should be structured.

### Task 1: Write a Setup Script  

Create a new file in `scripts/src01_setup.R`. This file will install and load any libraries we need across CRAN and Bioconductor.

```r
# ==============================================================================
# Script: scripts/src01_setup.R
# Purpose: Configure dependencies and environment for TCGA transcriptomic workflow
# ==============================================================================

# We first verify if BiocManager is available on our machine; if not, we install it 
# so we can seamlessly manage packages hosted on Bioconductor.
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# We list the standard CRAN packages required for data manipulation and visualisation.
cran_packages <- c("dplyr", "tidyr", "ggplot2", "ggrepel")

# We list the specialised Bioconductor packages needed for fetching GDC data, 
# managing expression matrices, and running differential expression.
bioc_packages <- c("TCGAbiolinks", "SummarizedExperiment", "limma", "edgeR")

# We identify which CRAN packages are currently missing from our library.
new_cran <- cran_packages[!(cran_packages %in% installed.packages()[, "Package"])]

# We install missing CRAN packages if needed.
if (length(new_cran) > 0) {
  install.packages(new_cran)
}

# We identify which Bioconductor packages are missing from our library.
new_bioc <- bioc_packages[!(bioc_packages %in% installed.packages()[, "Package"])]

# We install missing Bioconductor packages if needed.
if (length(new_bioc) > 0) {
  BiocManager::install(new_bioc)
}

message("Environment successfully setup! All dependencies are installed.")
```  

### Git Checkpoint:  
Stage your files in RStudio's Git tab, write a commit message ("Setup directory structure and setup script"), and push to GitHub.  

---

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

## Module 3: Statistical Analysis & PCA

In this module, we transition from data preparation to answering fundamental biological questions. Our primary objective is to identify which specific genes are "turned up" (up-regulated) or "turned down" (down-regulated) in tumour tissues compared to healthy tissues. In bioinformatics, this is known as Differential Expression analysis. This is how we ask questions like "what genes are expressed differently between Group A and Group B?", or, "what genes are more expressed as Variable A increases?". Cancer effects, gene mutation effects, age effects, sex effects, treatment effects, and many many more.  

To achieve this, we rely on a highly regarded statistical framework in R called `limma`. Because biological data is inherently noisy and we often have a limited number of clinical samples, `limma` employs an advanced technique called "empirical Bayes". Put simply, the software borrows statistical confidence from the thousands of genes present in the entire dataset to make much more accurate and robust decisions about each individual gene's behaviour.  

Furthermore, when we ask the question "is this gene significantly different?" across 20,000 genes simultaneously, basic probability breaks down. By pure chance alone, standard p-values would yield hundreds of false positives. To combat this "multiple testing" problem, we calculate adjusted p-values using the False Discovery Rate (FDR). This places a strict mathematical penalty on our results, ensuring our final list of genes represents genuine biological discoveries rather than statistical flukes.  

Finally, we will perform Principal Component Analysis (PCA). Our brains cannot easily visualise 20,000 variables (genes) at once. PCA is a dimensionality reduction technique that mathematically compresses all that complex, multi-dimensional gene expression data down into just two or three main axes of variation. By plotting these "Principal Components", we can visually inspect whether our tumour samples naturally cluster away from our normal samples. This serves as a vital quality control step, confirming that the underlying disease biology is driving the main differences in our data. It is also how we inspect at scale other effects across many dimensions simultaneously.  

Create `scripts/src03_analysis.R` to run this analytical workflow:  

```r
# ==============================================================================
# Script: scripts/src03_analysis.R
# Purpose: Perform limma differential expression and Principal Component Analysis
# ==============================================================================

# We load our analytical packages quietly.
suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(tibble)
  library(matrixStats)
})

# We read in our preprocessed data artifacts from disk.
sample_metadata <- read.csv("data/sample_metadata.csv")
gene_info       <- read.csv("data/gene_info.csv", row.names = 1)
log_expr        <- readRDS("data/normalised_log_expression.rds")
v               <- readRDS("data/voom_object.rds")

# Ensure factor ordering is maintained
sample_metadata$Condition <- factor(sample_metadata$Condition, levels = c("Normal", "Tumour"))

# ------------------------------------------------------------------------------
# 1. Differential Expression Analysis (limma)
# ------------------------------------------------------------------------------

# We build a linear model design matrix parameterised by sample condition.
design <- model.matrix(~ Condition, data = sample_metadata)

# We fit gene-wise linear models using our voom-weighted expression object.
fit <- lmFit(v, design)

# We apply empirical Bayes moderation to shrink gene-wise variances toward a global trend, 
# stabilising inference even with modest sample sizes.
fit <- eBayes(fit)

# We extract differential expression statistics across all tested genes and merge 
# official gene symbol annotations.
de_results <- topTable(fit, coef = "ConditionTumour", number = Inf) %>%
  rownames_to_column("gene_id") %>%
  left_join(gene_info, by = "gene_id") %>%
  mutate(
    gene_name = ifelse(is.na(gene_name), gene_id, gene_name),
    Significance = case_when(
      adj.P.Val < 0.05 & logFC > 1.5  ~ "Up-regulated",
      adj.P.Val < 0.05 & logFC < -1.5 ~ "Down-regulated",
      TRUE                            ~ "Not Significant"
    )
  )

# ------------------------------------------------------------------------------
# 2. Dimensionality Reduction (PCA)
# ------------------------------------------------------------------------------

# We isolate the top 500 most variable genes across our dataset to capture major variance.
var_genes <- head(order(rowVars(log_expr), decreasing = TRUE), 500)

# We perform Principal Component Analysis on transposed expression data for these variable genes.
pca_res <- prcomp(t(log_expr[var_genes, ]), scale. = TRUE)

# We compute the proportion of variance explained by each principal component.
var_explained <- summary(pca_res)$importance[2, ] * 100

# We structure our PCA coordinates into a data frame merged with sample condition labels.
pca_df <- as.data.frame(pca_res$x) %>%
  rownames_to_column("barcode") %>%
  inner_join(sample_metadata, by = "barcode")

# We attach variance metrics as an attribute so our visualisation script can format axis labels dynamically.
attr(pca_df, "var_explained") <- var_explained

# ------------------------------------------------------------------------------
# 3. Save Analytical Outputs
# ------------------------------------------------------------------------------

# We save our statistical results and PCA tables to our output directory.
write.csv(de_results, "output/de_results.csv", row.names = FALSE)
write.csv(pca_df, "output/pca_results.csv", row.names = FALSE)
saveRDS(var_explained, "output/variance_explained.rds")

message("Statistical modeling and PCA complete! Results exported to output/")
```

---

## Module 4: Data Visualisation

In this final module, we focus on translating our complex statistical outputs into intuitive, publication-ready visualisations. Staring at a spreadsheet containing 20,000 p-values is rarely illuminating; to effectively communicate our findings, we must create clear visual narratives. 

To achieve this, we will rely heavily on `ggplot2`, the gold-standard graphics package in R. `ggplot2` operates on a philosophy of "layers" - we start with a blank canvas, map our data to the axes, add geometric shapes (like points), and finally customise the themes and colours. We will generate two fundamental bioinformatic plots:

**1. The PCA Scatter Plot:** This gives us a bird's-eye view of our entire dataset. Each point on this plot represents an entire patient sample. If our biological groups are truly distinct, we should see the Tumour dots cluster tightly together and separate cleanly from the Normal dots. It is the ultimate sanity check, confirming that the main source of variation in our experiment is actually the disease state, rather than a technical batch effect.

**2. The Volcano Plot:** This is the classic method for visualising differential expression. The horizontal x-axis displays the *magnitude* of the biological change (Log2 Fold Change—how much expression increased or decreased), while the vertical y-axis displays the *statistical significance* (adjusted p-value). The resulting shape resembles an erupting volcano. The most biologically interesting target genes—those that are both highly significant and drastically altered—are pushed to the top-left and top-right corners. We will also use an excellent package called `ggrepel` to automatically arrange the names of these top genes so their text labels do not overlap.

Create `scripts/src04_visualisation.R` to construct these high-resolution figures:

```r
# ==============================================================================
# Script: scripts/src04_visualisation.R
# Purpose: Generate publication-ready PCA and Volcano plots
# ==============================================================================

# We load our visualisation packages quietly.
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

# We define our target cohort clean text label for figure headings.
cancer_label <- "Breast Cancer (BRCA)"

# We read in our processed analytical tables and metrics.
de_results    <- read.csv("output/de_results.csv")
pca_df        <- read.csv("output/pca_results.csv")
var_explained <- readRDS("output/variance_explained.rds")

# Ensure factor ordering is consistent for plot legends
pca_df$Condition <- factor(pca_df$Condition, levels = c("Normal", "Tumour"))

# ------------------------------------------------------------------------------
# 1. Principal Component Analysis (PCA) Plot
# ------------------------------------------------------------------------------

# We construct our PCA scatter plot to inspect global sample clustering.
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3.5, alpha = 0.85) +
  scale_color_manual(values = c("Normal" = "#2b5c8f", "Tumour" = "#d95f02")) +
  theme_minimal(base_size = 12) +
  labs(
    title = sprintf("PCA: %s (Tumour vs Normal)", cancer_label),
    x = sprintf("PC1 (%.1f%% Variance)", var_explained[1]),
    y = sprintf("PC2 (%.1f%% Variance)", var_explained[2]),
    color = "Group"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top"
  )

# We display the PCA plot in our session and save a high-resolution copy to disk.
print(p_pca)
ggsave("output/pca_plot.png", p_pca, width = 6, height = 5, dpi = 300)

# ------------------------------------------------------------------------------
# 2. Volcano Plot Visualisation
# ------------------------------------------------------------------------------

# We filter and select the top 10 most statistically significant genes for label positioning.
top_genes <- de_results %>%
  filter(Significance != "Not Significant") %>%
  arrange(adj.P.Val) %>%
  head(10)

# We construct our volcano plot mapping log2 fold-changes against -log10 p-values.
p_volcano <- ggplot(de_results, aes(x = logFC, y = -log10(P.Value), color = Significance)) +
  geom_point(alpha = 0.5, size = 1.8) +
  scale_color_manual(values = c(
    "Down-regulated" = "#3182bd", 
    "Not Significant" = "#969696", 
    "Up-regulated" = "#de2d26"
  )) +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  geom_text_repel(
    data = top_genes, 
    aes(label = gene_name), 
    size = 3.5, 
    max.overlaps = 15, 
    show.legend = FALSE
  ) +
  theme_classic(base_size = 15) +
  labs(
    title = sprintf("Differential Expression: %s", cancer_label),
    x = expression("Log"[2] ~ "Fold Change"),
    y = expression("-Log"[10] ~ "p-value")
  ) +
  guides(colour = guide_legend(override.aes = list(size = 6))) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top"
  )

# We display the volcano plot in our session and save a high-resolution copy to disk.
print(p_volcano)
ggsave("output/volcano_plot.png", p_volcano, width = 7, height = 6, dpi = 300)

message("Visualisations complete! PNG figures saved to output/")
```

---

## Module 5: Targeted Gene Set Analysis & Visualisation  

In this module, we will move away from looking at the entire transcriptome and focus on a specific subset of genes. Often in biological research, you have a predefined list of genes of interest, perhaps drawn from the literature, a specific biological pathway, or a therapeutic target panel. Many papers will investigate global differences between two datasets (and adjust heavily for false discovery, meaning stringent p value cutoffs) - however, we may have a more specific, hypothesis-driven research question - what do type-I interferon genes look like in older BRCA cancer patients compared to younger BRCA cancer patients? Are senescence genes associated with inflammatory gene responses in paired tumour and non-tumour samples? And many many more that we can investigate.    

We will learn how to programmatically extract a specific gene set from our large expression matrix. For this example, we will investigate a panel of immune checkpoint markers and cellular senescence (SASP) genes, which play critical roles in the tumour microenvironment and the ageing immune system. Once extracted, we will visualise the expression differences between Tumour and Normal samples using two highly effective methods:

1. **Grouped Boxplots (`ggplot2`)**: To directly compare the statistical distribution of expression values for each individual gene.
2. **Clustered Heatmap (`pheatmap`)**: To visualise the expression pattern of the entire gene set across all samples simultaneously, clustering patients with similar expression profiles together.

Create `scripts/src05_geneset_analysis.R` to run this targeted analysis:

```r
# ==============================================================================
# Script: scripts/src05_geneset_analysis.R
# Purpose: Extract a target gene set, generate grouped boxplots and a heatmap
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Setup & Data Loading
# ------------------------------------------------------------------------------

# Ensure pheatmap is installed for our heatmap visualisation
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
})

# We read in our preprocessed data artifacts from Module 2
sample_metadata <- read.csv("data/sample_metadata.csv")
gene_info       <- read.csv("data/gene_info.csv", row.names = 1)
log_expr        <- readRDS("data/normalised_log_expression.rds")

# We ensure our condition factor is correctly ordered for plotting
sample_metadata$Condition <- factor(sample_metadata$Condition, levels = c("Normal", "Tumour"))

# ------------------------------------------------------------------------------
# 2. Gene Set Extraction
# ------------------------------------------------------------------------------

# We define a targeted list of immune checkpoint and cellular senescence genes.
target_genes <- c("CD274", "CTLA4", "PDCD1", "LAG3", "CDKN1A", "CDKN2A", "IL6", "CXCL8")

# We filter our gene metadata to find the Ensembl IDs that match our target symbols
target_info <- gene_info %>% 
  filter(gene_name %in% target_genes)

# We subset our large log-expression matrix using these specific Ensembl IDs
mat_subset <- log_expr[rownames(log_expr) %in% target_info$gene_id, ]

# For intuitive visualisation, we replace the Ensembl IDs with human-readable gene symbols
rownames(mat_subset) <- target_info$gene_name[match(rownames(mat_subset), target_info$gene_id)]

# ------------------------------------------------------------------------------
# 3. Grouped Boxplots (ggplot2)
# ------------------------------------------------------------------------------

# To plot with ggplot2, we must reshape our matrix from "wide" to "long" format.
long_df <- as.data.frame(mat_subset) %>%
  tibble::rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene, names_to = "barcode", values_to = "Expression") %>%
  left_join(sample_metadata, by = "barcode")

p_box <- ggplot(long_df, aes(x = Gene, y = Expression, fill = Condition)) +
  geom_boxplot(alpha = 0.8, outlier.size = 1.2) +
  scale_fill_manual(values = c("Normal" = "#2b5c8f", "Tumour" = "#d95f02")) +
  theme_classic(base_size = 14) +
  labs(
    title = "Expression of Immune & Senescence Markers",
    x = "Gene Symbol",
    y = expression("Log"[2] ~ "Expression (CPM)")
  ) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

print(p_box)
ggsave("output/geneset_boxplot.png", p_box, width = 8, height = 5, dpi = 300)

# ------------------------------------------------------------------------------
# 4. Clustered Heatmap (pheatmap)
# ------------------------------------------------------------------------------

# We prepare an annotation data frame so pheatmap can colour-code our samples by Condition.
anno_df <- data.frame(Condition = sample_metadata$Condition)
rownames(anno_df) <- sample_metadata$barcode

# We define a custom colour palette for our condition annotations to keep styling consistent.
anno_colors <- list(
  Condition = c("Normal" = "#2b5c8f", "Tumour" = "#d95f02")
)

# We generate the clustered heatmap, scaling the rows (genes) so relative differences pop out.
# pheatmap automatically saves the file if we provide a filename argument.
pheatmap(
  mat_subset,
  annotation_col = anno_df,
  annotation_colors = anno_colors,
  scale = "row",
  show_colnames = FALSE,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  main = "Targeted Gene Expression Profile",
  filename = "output/geneset_heatmap.png",
  width = 7, 
  height = 5
)

message("Gene set visualisations complete! PNG figures saved to output/")

```

---

## Final Checklist

To complete this mini-project, ensure your public GitHub repository contains:

- [ ] Clear folder structure (`data/`, `scripts/`, `output/`).
- [ ] Four well-commented R scripts (`src01_setup.R` through `src04_visualisation.R`).
- [ ] Saved CSV outputs and data objects in `output/` and `data/`.
- [ ] Rendered `.png` figures linked directly inside your main `README.md`.

---  

## Task 1: Independent Assignment (Tumour vs. Tumour Comparison)

Now that we have successfully built a complete bioinformatic pipeline from raw data acquisition to publication-ready visualisations, it is time to test your new computational skills. 

Up until this point, we have compared **Tumour vs. Normal** tissue within a single cancer type (Breast Cancer, TCGA-BRCA). For your independent assignment, your objective is to investigate the transcriptomic differences between **two different types of cancer**. You will query a new TCGA cohort, isolate the tumour samples from both cohorts, and perform a differential expression analysis comparing **BRCA Tumours vs. [Your Chosen Cancer] Tumours**. More specifically, you will then investigate if there are differences in hypoxia gene signalling between the cancer types. 

This is a very common exploratory analysis in cancer biology, often used to identify tissue-specific oncogenic drivers or potential pan-cancer therapeutic targets.

### The Task

1. **Select a New Cohort:** Choose a new TCGA project to compare against our existing BRCA dataset. Good options include Lung Adenocarcinoma (`TCGA-LUAD`), Colon Adenocarcinoma (`TCGA-COAD`), or Prostate Adenocarcinoma (`TCGA-PRAD`).
2. **Data Acquisition:** Write a script to query and download RNA-Seq data for your new cohort. 
3. **Data Wrangling:** Filter both your BRCA dataset and your new dataset to retain **only** the "Primary Tumor" samples (drop the "Solid Tissue Normal" samples).
4. **Merge Datasets:** Combine the two tumour count matrices into a single matrix, and combine their metadata into a single table.
5. **Analysis & Visualisation:** Run the combined dataset through your `edgeR`/`limma` pipeline. Generate a PCA plot to see if the two cancer types cluster independently, and a Volcano plot highlighting the top differentially expressed genes between the two malignancies.
6. **Targeted Gene Analysis:** Pull out the transcription levels of genes from the hypoxia hallmark gene set, then compare levels between your groups.

---

### Tips and Hints for Success

Merging two independent datasets can introduce a few computational hurdles. Here are some tips to help you navigate the tricky parts of this assignment:

* **Hint 1: Consistent Gene Lists (The `intersect` function)**
  When you download two different TCGA projects, the exact list of genes (the rownames of your count matrices) might not perfectly match due to slight differences in genome annotations or zero-count filtering. Before binding the matrices together using `cbind()`, you must ensure they have the exact same genes in the exact same order. Use the `intersect()` function to find the overlapping genes, and subset both matrices before merging:  

```r
  common_genes <- intersect(rownames(brca_counts), rownames(luad_counts))
  brca_counts_sub <- brca_counts[common_genes, ]
  luad_counts_sub <- luad_counts[common_genes, ]
  combined_counts <- cbind(brca_counts_sub, luad_counts_sub)
```

* **Hint 2: Updating the Metadata Factor**
In Module 3, our experimental Condition was "Normal" vs "Tumour". For this assignment, you will need to create a new column in your merged metadata table (e.g., Cancer_Type) to act as your new condition. Make sure it is properly formatted as a factor before building your model.matrix():


```r
merged_metadata$Cancer_Type <- factor(merged_metadata$Cancer_Type, levels = c("BRCA", "LUAD"))

```  

* **Hint 3: Memory Management**

Combining two TCGA cohorts can result in very large objects that consume a lot of RAM. Stick to subsetting around 20–30 tumour samples per cancer type to keep your code running smoothly on a standard laptop.  


* **Hint 4: A Note on Batch Effects**  

In real-world bioinformatics, comparing datasets generated from different distinct projects can introduce technical "batch effects" (e.g., different sequencing machines, different technicians, different extraction kits). While we are skipping formal batch correction for this introductory exercise, it is an important biological caveat.    


---



