# Mini-Project: Introduction to R, Git, and Basic Cancer Transcriptomic Analyses  

[![Cohort](https://img.shields.io/badge/Cohort-TCGA-blue.svg)](https://portal.gdc.cancer.gov/)
[![Cohort](https://img.shields.io/badge/Cohort-TCGA%20(cBioPortal)-008080.svg)](https://www.cbioportal.org/)
[![Cohort](https://img.shields.io/badge/Cohort-TCGA--BRCA-darkblue.svg)](https://portal.gdc.cancer.gov/projects/TCGA-BRCA)


Welcome to your introductory dry-lab mini-project. This hands-on tutorial is designed to build foundational skills in **computational biology**, **version control**, and **data science using R**, designed for someone with a wet lab background interested in adding some computational skills to their week-to-week as a researcher.  


### Learning Objectives
By completing this project, you will learn how to:
1. Set up a reproducible project structure using **RStudio** and **Git/GitHub**.
2. Perform basic data wrangling, quality control, and log-transformation.
3. Conduct hypothesis testing ($t$-tests), multiple testing correction (FDR), and dimensionality reduction (**PCA**).
4. Generate publication-ready visualisations (**PCA plots** & **Volcano plots**) using `ggplot2`.

---  

## Software & Setup

Please complete the following installations before starting Module 1:

* [ ] Install **R (v4.3+)**: [cran.r-project.org](https://cran.r-project.org/)
* [ ] Install **RStudio Desktop**: [posit.co/download/rstudio-desktop](https://posit.co/download/rstudio-desktop/)
* [ ] Install **Git**: [git-scm.com](https://git-scm.com/)
* [ ] Create a **GitHub** account: [github.com](https://github.com/)

Installing R & RStudio should be relatively straightforward, although it can be a little slow depending on your machine.  RStudio is our graphical interface which runs R, and R is a programming language and software environment designed for statistical computing, data analysis, and graphical visualisation. It is really commonly used in Bioinformatics, similar to Python. In my opinion, R is better for Viz and standards statistics, not as good for bigger scale machine learning projects, but overall much more forgiving for someone first experiencing coding.  


GitHub is a cloud-based platform that hosts software development projects and allows us as researchers or code developers to store, manage, track, and share code using Git, a distributed version control system. Think of it like a private or public Sharepoint/Onedrive where we can share our code documents with lab/team members or the wider public. In R, we often install and use "Packages" of code - GitHub is where one might host a package, or something as simple as a single chunk of code that we are working on. "Git" specifically is a program we use to send/pull our files to/from GitHub.  

---  

### Recommended Preparation: Software Carpentry
If you are completely new to programming, I strongly recommend completing the free, self-paced **Software Carpentry: Programming with R** course before or alongside this project. It provides a beginner-friendly overview of basic R syntax, data structures, and best practices using a clinical example:

* **Course Link:** [Software Carpentry – Programming with R](https://swcarpentry.github.io/r-novice-inflammation/)

---  

## Module 1: Project Setup & Version Control

In this module, we will set up a structured computational environment and establish version control using Git and GitHub. Organizing your code, data, and outputs into standard subdirectories from the outset ensures your project remains structured, while linking your RProject to GitHub ensures all updates are tracked, backed up, and fully reproducible. We want our code to work time and time again on our machine, AND to be useable in the exact same way on someone else's computer (a colleague, a reviewer, a client). Similar to lab work, we want a strict lab notebook with a very specific layout and set of instructions, and for it to actually be readable for others! Think of recipe books - they usually follow the same structure - we do not want to be on the final step of baking our cake to realize we are missing a vital ingredient, coding follows the same logic.

### Step 1: Create local directory structure
Open **RStudio**, create a new R Project linked to your cloned GitHub repository, and run the following in your Console:

```r
# Create standard project folders
dir.create("data")
dir.create("scripts")
dir.create("output")

```

Now, when someone else opens our script, they will know how our data/input/output should be structured.

# Task 1: Write a Script  

Create a new file in `scripts/srcsrc01_setup.R`. This file will install and load any libraries we need.


### Script: src01_setup.R
## Purpose: Check and install required packages

When we start using RStudio, we will need to install packages, using `install.packages("package_name")`. We then switch on the package, using `library(package_name)`. Note the use of quotations for install but not for library. The next time we open up RStudio, we do not need to reinstall the package, but we do need to switch it on, with the `library()` function.  

The below code gives R the list of packages we need and installs them if they are missing. 

```r
required_packages <- c("dplyr", "ggplot2", "tidyverse", "ggrepel")

new_packages <- required_packages[!(required_packages %in% installed_packages)]
if(length(new_packages)) install.packages(new_packages)

message("Packages installed. Environment setup complete.")

```  


### Git Checkpoint:  
Stage your files in RStudio's Git tab, write a commit message ("Setup directory structure and setup script"), and push to GitHub.  


---

## Module 2: Preprocessing & Quality Control

Create `scripts/src02_preprocessing.R` to download, filter, and normalize real breast cancer expression data from TCGA:

```r
# Script: src02_preprocessing.R
# Purpose: Fetch real TCGA-BRCA transcriptomic data, extract metadata, and log2-transform counts

# 1. Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) BiocManager::install("TCGAbiolinks")
if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) BiocManager::install("SummarizedExperiment")

library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)

message("Downloading TCGA-BRCA data from GDC portal...")

# 2. Query GDC for TCGA-BRCA RNA-seq gene expression (STAR - Counts)
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Download a subset of samples (20 Primary Tumour, 20 Solid Tissue Normal) to keep runtime fast
samples_normal <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
# Download data files locally
GDCdownload(query, files.per.chunk = 10)
brca_se <- GDCprepare(query)

# 3. Extract sample metadata & define conditions
coldata <- as.data.frame(colData(brca_se)) %>%
  select(barcode, sample_type) %>%
  filter(sample_type %in% c("Primary Tumor", "Solid Tissue Normal")) %>%
  mutate(condition = factor(ifelse(sample_type == "Primary Tumor", "Tumour", "Normal"), 
                            levels = c("Normal", "Tumour")))

# Select a balanced subset: 20 Normal and 20 Tumour samples
set.seed(42)
selected_samples <- coldata %>%
  group_by(condition) %>%
  slice_sample(n = 20) %>%
  ungroup()

# 4. Filter count matrix to selected samples & top variable genes
raw_counts <- assay(brca_se, "unstranded")[, selected_samples$barcode]

# Clean gene symbols from row data
gene_info <- as.data.frame(rowData(brca_se))
rownames(raw_counts) <- gene_info$gene_name

# Filter out unmapped/empty gene symbols and low-count genes
raw_counts <- raw_counts[!is.na(rownames(raw_counts)) & rownames(raw_counts) != "", ]
raw_counts <- raw_counts[rowSums(raw_counts) > 50, ]

# Select top 1,000 most variable genes for efficient processing
gene_vars <- apply(raw_counts, 1, var)
top_genes <- head(order(gene_vars, decreasing = TRUE), 1000)
counts_subset <- raw_counts[top_genes, ]

# 5. Quality Control & Normalisation
stopifnot(sum(is.na(counts_subset)) == 0) # Confirm zero missing values

# Log2 transformation (adding pseudo-count +1)
log_counts <- log2(counts_subset + 1)

# Format sample IDs for clean downstream handling
metadata_clean <- selected_samples %>%
  transmute(
    sample_id = barcode,
    condition = condition
  )

colnames(log_counts) <- metadata_clean$sample_id

# 6. Save outputs to disk
write.csv(metadata_clean, "data/sample_metadata.csv", row.names = FALSE)
write.csv(log_counts, "data/processed_counts.csv", row.names = TRUE)

message("Preprocessing complete! Real TCGA-BRCA dataset saved to data/")
```  

---

# Module 3: Statistical Analysis & PCA

```markdown
## Module 3: Statistical Testing & Dimensionality Reduction

Create `scripts/src03_analysis.R` to calculate differential expression and run Principal Component Analysis (PCA):

```r
# Script: src03_analysis.R
# Purpose: t-tests, FDR correction, and PCA

library(tidyverse)

metadata <- read.csv("data/sample_metadata.csv")
log_counts <- read.csv("data/processed_counts.csv", row.names = 1)

# 1. Differential Expression Testing
de_results <- apply(log_counts, 1, function(gene_vals) {
  norm <- gene_vals[metadata$condition == "Normal"]
  tum  <- gene_vals[metadata$condition == "Tumour"]
  
  l2fc <- mean(tum) - mean(norm)
  pval <- t.test(tum, norm)$p.value
  
  return(c(log2FoldChange = l2fc, pvalue = pval))
})

# Reformat and apply Benjamini-Hochberg (FDR) correction
de_df <- as.data.frame(t(de_results)) %>%
  rownames_to_column(var = "gene") %>%
  mutate(
    padj = p.adjust(pvalue, method = "BH"),
    significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Significant", "Not Significant")
  )

# 2. Principal Component Analysis
pca_res <- prcomp(t(log_counts), scale. = TRUE)

pca_df <- as.data.frame(pca_res$x) %>%
  rownames_to_column(var = "sample_id") %>%
  left_join(metadata, by = "sample_id")

# Save tables
write.csv(de_df, "output/de_results.csv", row.names = FALSE)
write.csv(pca_df, "output/pca_results.csv", row.names = FALSE)

message("Analysis complete. Output saved to output/")

```

---

# Module 4: Data Visualisation

```markdown
## Module 4: Data Visualisation (`ggplot2`)

Create `scripts/src04_visualisation.R` to construct high-resolution plots:

```r
# Script: src04_visualisation.R
# Purpose: Generate PCA plot and Volcano plot

library(tidyverse)

de_df  <- read.csv("output/de_results.csv")
pca_df <- read.csv("output/pca_results.csv")

# 1. PCA Plot
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("Normal" = "#1f77b4", "Tumour" = "#d62728")) +
  theme_bw() +
  labs(title = "PCA: Normal vs Tumour", x = "PC1", y = "PC2", color = "Group")

ggsave("output/pca_plot.png", p_pca, width = 6, height = 5, dpi = 300)

# 2. Volcano Plot
p_volcano <- ggplot(de_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "#d62728")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 p-value")

ggsave("output/volcano_plot.png", p_volcano, width = 6, height = 5, dpi = 300)

message("Visualisations created successfully.")

```

---

# README Final Deliverable Checklist

```markdown
## Final Checklist

To complete this mini-project, ensure your public GitHub repository contains:

- [ ] Clear folder structure (`data/`, `scripts/`, `output/`).
- [ ] Four well-commented R scripts (`src01_setup.R` through `src04_visualisation.R`).
- [ ] Saved CSV outputs in `output/`.
- [ ] Rendered `.png` figures linked directly inside your main `README.md`.
