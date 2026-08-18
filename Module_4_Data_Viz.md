
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
