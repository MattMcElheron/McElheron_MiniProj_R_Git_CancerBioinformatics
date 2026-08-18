## Task 1: Independent Assignment (Tumour vs. Tumour Comparison)

Now that we have successfully built a complete bioinformatic pipeline from raw data acquisition to publication-ready visualisations, it is time to test your new computational skills. 

Up until this point, we have compared **Tumour vs. Normal** tissue within a single cancer type (Breast Cancer, TCGA-BRCA). For your independent assignment, your objective is to investigate the transcriptomic differences between **two different types of cancer**. You will query a new TCGA cohort, isolate the tumour samples from both cohorts, and perform a differential expression analysis comparing **BRCA Tumours vs. [Your Chosen Cancer] Tumours**. More specifically, you will then investigate if there are differences in hypoxia gene signalling between the cancer types. 

This is a very common exploratory analysis in cancer biology, often used to identify tissue-specific oncogenic drivers or potential pan-cancer therapeutic targets.  

### The Task

1. **Select a New Cohort:** Choose a new TCGA project to compare against our existing BRCA dataset. Good options include Lung Adenocarcinoma (`TCGA-LUAD`), Colon Adenocarcinoma (`TCGA-COAD`), or Prostate Adenocarcinoma (`TCGA-PRAD`).
2. **Data Acquisition:** Write a script to query and download RNA-Seq data for your new cohort. Sample number here doesn't matter too much, the script is what we are interested in here.
3. **Data Wrangling:** Filter both your BRCA dataset and your new dataset to retain **only** the "Primary Tumor" samples (drop the "Solid Tissue Normal" samples).
4. **Merge Datasets:** Combine the two tumour count matrices into a single matrix, and combine their metadata into a single table.
5. **Analysis & Visualisation:** Run the combined dataset through your `edgeR`/`limma` pipeline. Generate a PCA plot to see if the two cancer types cluster independently, and a Volcano plot highlighting the top differentially expressed genes between the two malignancies.
6. **Targeted Gene Analysis:** Pull out the transcription levels of genes from the hypoxia hallmark gene set, then compare levels between your groups.  

---  

### Preparing a Bioinformatics Report  

When writing up your findings for this assignment, I recommend adopting the standard format of a scientific research article. Your report should guide the reader logically from raw data retrieval to biological interpretation. Below is a suggested structure for your submission:

1. **Introduction:** Briefly introduce the two selected cancer types and outline the biological rationale for comparing their transcriptomes. Introduce hypoxia signalling as a critical microenvironmental factor influencing tumour progression, therapy resistance, and metabolic reprogramming.
2. **Methods:** Clearly document your computational workflow. Specify the TCGA cohorts queried, the filtering and normalisation steps applied (such as TMM normalisation and voom transformation), and the statistical models used for differential expression testing (limma).
3. **Results:** Present your findings with high-resolution visualisations. Make sure and create caption your figure with legends.
4. **Discussion:** Interpret your results in a broader biological context, including potential limitaitons to your work and future steps.  

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

---

## Assignment Final Checklist

To complete this mini-project assignment, you will need to complete and document the following (preferably on your GitHub):

- [ ] **Choose a New Cohort:** Select an additional TCGA project to compare against our existing BRCA dataset (e.g., TCGA-LUAD, TCGA-COAD, or TCGA-PRAD).
- [ ] **Data Acquisition Script:** Write or adapt an R script to query and download RNA-Seq quantification files for your chosen cohort via TCGAbiolinks.
- [ ] **Filter Normal Samples:** Filter both your BRCA dataset and your new dataset to retain only "Primary Tumor" samples, dropping all "Solid Tissue Normal" samples.
- [ ] **Merge Datasets:** Combine the two tumour count matrices using intersect() to match genes correctly, and merge their sample metadata tables.
- [ ] **Differential Expression & PCA:** Run your combined data through the edgeR/limma pipeline, generate a PCA plot to inspect whether the two cancer types cluster independently, construct a Volcano plot, and investigate your gene set of interest.
- [ ] **Final GitHub Deliverable:** Ensure your public repository contains the updated script(s), saved analytical output tables, rendered .png figures, and a brief note in your README.md discussion what you have found.  

---  





