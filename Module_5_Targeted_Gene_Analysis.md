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
