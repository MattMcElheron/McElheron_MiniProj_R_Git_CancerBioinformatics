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

Create a new file in `scripts/src01_setup.R`. This file will install and load any libraries we need.


### Script: 01_setup.R
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

# Module 2: Data Preprocessing & QC

```markdown
## Module 2: Preprocessing & Quality Control

Create `scripts/02_preprocessing.R` to simulate and clean a transcriptomic expression matrix:

```r
# Script: 02_preprocessing.R
# Purpose: Data simulation, quality control, and log2 transformation

library(tidyverse)
set.seed(42)

# 1. Simulate metadata (20 Normal vs 20 Tumour samples)
metadata <- data.frame(
  sample_id = paste0("Sample_", 1:40),
  condition = factor(rep(c("Normal", "Tumour"), each = 20), levels = c("Normal", "Tumour"))
)

# 2. Simulate count matrix (500 genes x 40 samples)
counts <- matrix(rpois(500 * 40, lambda = 50), nrow = 500, ncol = 40)
rownames(counts) <- paste0("Gene_", 1:500)
colnames(counts) <- metadata$sample_id

# Introduce artificial upregulation in Tumour for Genes 1-50
counts[1:50, metadata$condition == "Tumour"] <- counts[1:50, metadata$condition == "Tumour"] * 2.5

# 3. Quality Control Checks
stopifnot(sum(is.na(counts)) == 0) # Check no missing values

# Log2 transformation (adding pseudo-count +1)
log_counts <- log2(counts + 1)

# 4. Save outputs
write.csv(metadata, "data/sample_metadata.csv", row.names = FALSE)
write.csv(log_counts, "data/processed_counts.csv", row.names = TRUE)

message("Preprocessing complete. Output saved to data/")
```
