<table>
  <tr>
    <td style="vertical-align: middle; padding-right: 18px;">
      <img src="Logo.jpeg" alt="Logo" width="100"/>
    </td>
    <td style="vertical-align: middle;">
      <h1 style="margin-bottom: 0;">Key Bioinformatics Tutorials</h1>
      <em>Bioinformatics Data and Tutorial Guides</em>
    </td>
  </tr>
</table>

---

# Key Bioinformatics Tutorials  
*Bioinformatics Data and Tutorial Guides*

---
## Bulk RNA-seq tools

Bulk RNA-seq analysis involves quantifying gene expression from pooled cell populations. Common tools include:

- **edgeR (R)**: Uses a negative binomial model to detect differentially expressed genes, especially suited for small sample sizes and overdispersed count data.  
  [edgeR User's Guide](https://bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)

- **DESeq2 (R)**: Applies shrinkage estimation for dispersion and fold changes to improve stability and interpretability of differential expression results.  
  [DESeq2 Vignette](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

- **limma-voom (R)**: Transforms count data to log-counts per million with associated precision weights, allowing linear modeling with limma's framework. Suitable for larger sample sizes.  
  [limma-voom Tutorial](https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html)

---

## Single-cell RNA-seq and ATAC-seq tools

- **scRNAseq Seurat (R)**: Comprehensive toolkit for single-cell RNA-seq data analysis, including clustering, visualization, and integration.  
  [Seurat Documentation](https://satijalab.org/seurat/)

- **SCVI-tools (Python)**: Probabilistic models for single-cell RNA-seq data analysis enabling batch correction, imputation, and clustering.  
  [SCVI-tools Documentation](https://docs.scvi-tools.org/en/1.0.0/index.html)

- **Signac (R)**: Toolkit for analyzing single-cell chromatin data, such as ATAC-seq, integrated with Seurat workflows.  
  [Signac PBMC Vignette](https://stuartlab.org/signac/articles/pbmc_vignette.html)
  
- **DoubletFinder (R)**: Tool for detecting doublets (artificial merged cells) in single-cell RNA-seq data to improve data quality.  
  [DoubletFinder GitHub](https://github.com/chris-mcginnis-ucsf/DoubletFinder)

- **Seurat Integration (R)**: Workflow for integrating multiple single-cell datasets to correct batch effects and combine data across experiments.  
  [Seurat Integration Tutorial](https://satijalab.org/seurat/articles/seurat5_integration)

---
  
## Mosaic Integration tools

Mosaic integration in bioinformatics refers to computational methods that integrate heterogeneous single-cell multi-omics datasets—such as transcriptomics, epigenomics, proteomics, and spatial data—often containing missing modalities. This approach enables comprehensive analysis across diverse data types, addressing challenges like modality scalability and batch effects, and is crucial for advancing multi-omics research.

- **MOFA (R)**: Multi-Omics Factor Analysis framework that identifies shared and dataset-specific factors across multi-modal data to facilitate integrative analysis.  
  [MOFA Documentation](https://biofam.github.io/MOFA2/)

- **LIGER (R)**: Uses integrative non-negative matrix factorization to jointly analyze multiple single-cell datasets, allowing identification of shared and dataset-specific features.  
  [LIGER GitHub](https://github.com/welch-lab/liger)

- **STABMAP (R)**: Statistical framework designed to integrate single-cell datasets by accounting for batch effects and other sources of variability.  
  [STABMAP GitHub](https://github.com/MarioniLab/StabMap)

- **MIDAS (Python)**: Toolkit focused on multi-omics integration for single-cell data, enabling discovery of cellular heterogeneity across modalities.  
  [MIDAS GitHub](https://github.com/labomics/midas)

---
  ## Spatial Transcriptomics Tools

### Image-based Data
- **Seurat Spatial (Image-based)**: Comprehensive workflow for analyzing image-based spatial transcriptomics data using Seurat.  
  [Seurat Spatial Vignette](https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2)

### VISIUM HD
- **Seurat VISIUM HD**: Analysis pipeline tailored for high-definition VISIUM spatial transcriptomics data.  
  [Seurat VISIUM HD Vignette](https://satijalab.org/seurat/articles/visiumhd_analysis_vignette)

### Sequencing-based Data
- **Seurat Spatial (Sequencing-based)**: Tutorial for processing and analyzing sequencing-based spatial transcriptomics data.  
  [Seurat Spatial Sequencing Vignette](https://satijalab.org/seurat/articles/spatial_vignette)

### Other Tools
- **lisaClust**: Bioconductor package for clustering and analyzing spatial transcriptomics data.  
  [lisaClust Bioconductor](https://www.bioconductor.org/packages/release/bioc/vignettes/lisaClust/inst/doc/lisaClust.html) | [lisaClust GitHub](https://github.com/SydneyBioX/lisaClust)

- **Banksy**: An R package for spatial transcriptomics analysis focusing on spatial autocorrelation and gene expression patterns.  
  [Banksy Documentation](https://prabhakarlab.github.io/Banksy/)
