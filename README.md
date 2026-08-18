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

---

## Module 1: Project Setup & Version Control

### Step 1: Create local directory structure
Open **RStudio**, create a new R Project linked to your cloned GitHub repository, and run the following in your Console:

```r
# Create standard project folders
dir.create("data")
dir.create("scripts")
dir.create("output")

```

# Script: 01_setup.R
# Purpose: Check and install required packages

required_packages <- c("tidyverse", "ggrepel")

new_packages <- required_packages[!(required_packages %in% installed_packages)]
if(length(new_packages)) install.packages(new_packages)

message("Environment setup complete.")








