
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
