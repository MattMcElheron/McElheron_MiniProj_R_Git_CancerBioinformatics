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
