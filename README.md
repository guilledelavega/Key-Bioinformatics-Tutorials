# Key Bioinformatics Tutorials  
*Bioinformatics Data and Tutorial Guides*

---
## Bulk RNA-seq

Bulk RNA-seq analysis involves quantifying gene expression from pooled cell populations. Common tools include:

- **edgeR (R)**: Uses a negative binomial model to detect differentially expressed genes, especially suited for small sample sizes and overdispersed count data.  
  [edgeR User's Guide](https://bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)

- **DESeq2 (R)**: Applies shrinkage estimation for dispersion and fold changes to improve stability and interpretability of differential expression results.  
  [DESeq2 Vignette](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

- **limma-voom (R)**: Transforms count data to log-counts per million with associated precision weights, allowing linear modeling with limma's framework. Suitable for larger sample sizes.  
  [limma-voom Tutorial](https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html)

---

## Single-cell RNA-seq and ATAC-seq Technologies

- **scRNAseq Seurat (R)**: Comprehensive toolkit for single-cell RNA-seq data analysis, including clustering, visualization, and integration.  
  [Seurat Documentation](https://satijalab.org/seurat/)

- **SCVI-tools (Python)**: Probabilistic models for single-cell RNA-seq data analysis enabling batch correction, imputation, and clustering.  
  [SCVI-tools Documentation](https://docs.scvi-tools.org/en/1.0.0/index.html)

- **Signac (R)**: Toolkit for analyzing single-cell chromatin data, such as ATAC-seq, integrated with Seurat workflows.  
  [Signac PBMC Vignette](https://stuartlab.org/signac/articles/pbmc_vignette.html)

---
  
## Mosaic Integration Tools

Mosaic integration in bioinformatics refers to computational methods that integrate heterogeneous single-cell multi-omics datasets—such as transcriptomics, epigenomics, proteomics, and spatial data—often containing missing modalities. This approach enables comprehensive analysis across diverse data types, addressing challenges like modality scalability and batch effects, and is crucial for advancing multi-omics research.

- **MOFA (R)**: Multi-Omics Factor Analysis framework that identifies shared and dataset-specific factors across multi-modal data to facilitate integrative analysis.  
  [MOFA Documentation](https://biofam.github.io/MOFA2/)

- **LIGER (R)**: Uses integrative non-negative matrix factorization to jointly analyze multiple single-cell datasets, allowing identification of shared and dataset-specific features.  
  [LIGER GitHub](https://github.com/welch-lab/liger)

- **STABMAP (R)**: Statistical framework designed to integrate single-cell datasets by accounting for batch effects and other sources of variability.  
  [STABMAP GitHub](https://github.com/MarioniLab/StabMap)

- **MIDAS (Python)**: Toolkit focused on multi-omics integration for single-cell data, enabling discovery of cellular heterogeneity across modalities.  
  [MIDAS GitHub](https://github.com/labomics/midas)


