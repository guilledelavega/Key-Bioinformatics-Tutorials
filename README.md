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

- **NGS Tutorials**: Several NGS tutorials for absolute begginers.  
  [NGS Tutorials](https://ngs101.com/)

# Bulk RNA-seq tools

Bulk RNA-seq analysis involves quantifying gene expression from pooled cell populations. Common tools include:

- **edgeR (R)**: Uses a negative binomial model to detect differentially expressed genes, especially suited for small sample sizes and overdispersed count data.  
  [edgeR User's Guide](https://bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)

- **DESeq2 (R)**: Applies shrinkage estimation for dispersion and fold changes to improve stability and interpretability of differential expression results.  
  [DESeq2 Vignette](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

- **limma-voom (R)**: Transforms count data to log-counts per million with associated precision weights, allowing linear modeling with limma's framework. Suitable for larger sample sizes.  
  [limma-voom Tutorial](https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html)

- **PUREE (Python)**: Compact and fast machine-learning method that estimates tumour purity (cancer-cell fraction) directly from bulk RNA-seq expression profiles.  
  [PUREE GitHub](https://github.com/skandlab/PUREE) 

- **decoupleR (R)**: Framework that bundles multiple statistical methods to infer transcription-factor and pathway activities from omics data.
  [decoupleR GitHub](https://github.com/saezlab/decoupleR) 

- **Gene Set Testing Tutorial (R)**: Hands-on guide for covering gene-set enrichment analysis (GSEA) on RNA-seq data.  
  [GSEA Tutorial](https://bioinformatics-core-shared-training.github.io/RNAseq_September_2019/html/06_Gene_set_testing.html)

- **ComBat-seq (R)**: Performs batch-effect correction on RNA-seq *count* data with a negative-binomial regression.  
  [ComBat-seq tutorial](https://rnabio.org/module-03-expression/0003/06/02/Batch-Correction/)

- **RUVSeq – Remove Unwanted Variation (R)**: Suite of methods that estimate hidden factors causing technical variation and incorporate them into downstream DE models for cleaner differential-expression results.  
  [RUVSeq vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.html)

- **WGCNA (R)**: Builds weighted gene-co-expression networks to cluster genes into modules and relate those modules to phenotypes, enabling discovery of regulatory programs in bulk RNA-seq datasets.  
  [WGCNA tutorial 1](https://bioinformaticsworkbook.org/dataAnalysis/RNA-Seq/RNA-SeqIntro/wgcna.html#gsc.tab=0) · [WGCNA tutorial 2](https://fuzzyatelin.github.io/bioanth-stats/module-F21-Group1/module-F21-Group1.html)

- **tximport (R)**: Imports transcript-level estimates from pseudo - genome aligners(e.g., from kallisto/salmon) and summarizes them to gene-level count/abundance matrices—optionally length-scaled—so they can be used directly in DESeq2/edgeR/limma while retaining proper variance information.
  [tximport tutorial](https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html)

- **variancePartition (R)**: General framework for understanding drivers of variation in gene expression in experiments with complex designs.
  [variancePartition Tutorial](https://www.bioconductor.org/packages/release/bioc/vignettes/variancePartition/inst/doc/variancePartition.html)

- **rrvo (R)**:  It aims at simplifying the redundance of GO sets by grouping similar terms based on their semantic similarity (parental function group).  
  [rrvgo · GitHub](https://ssayols.github.io/rrvgo/articles/rrvgo.html)

---
  
# Whole Exome Sequencing
- **maftools (R)**: Toolkit to analyze single-sample WES analysis, importing MAF/VCF, performing QC, computing metrics like TMB, and producing publication-ready plots (oncoplots, lollipop, summaries).  
  [maftools Bioconductor](https://bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html)

- **MesKit (R)**: Toolkit for multi-region/longitudinal tumor WES that integrates multiple biopsies from the same patient to infer clonal architecture, branch mutations, reconstruct phylogenetic trees, explore metastatic routes, and summarize mutational signatures.  
  [MesKit Bioconductor](https://www.bioconductor.org/packages/devel/bioc/vignettes/MesKit/inst/doc/MesKit.html#6_Measurement_of_intra-tumor_heterogeneity)

---
# Single-cell RNA-seq and ATAC-seq tools

- **scRNAseq Seurat (R)**: Comprehensive toolkit for single-cell RNA-seq data analysis, including clustering, visualization, and integration.  
  [Seurat Documentation](https://satijalab.org/seurat/)

- **Scanpy (Python)**: Comprehensive framework for the analysis of single-cell gene expression data.  
  [Scanpy Tutorial](https://scanpy.readthedocs.io/en/1.10.x/tutorials/spatial/basic-analysis.html)

- **OSCA (R)**: Orchestrating Single-Cell Analysis with Bioconductor book, which teaches users some common workflows for the analysis of single-cell RNA-seq data (scRNA-seq).  
[OSCA Bioconductor](https://bioconductor.org/books/release/OSCA/)

- **SCVI-tools (Python)**: Probabilistic models for single-cell RNA-seq data analysis enabling batch correction, imputation, and clustering.  
  [SCVI-tools Documentation](https://docs.scvi-tools.org/en/1.0.0/index.html)

- **Signac (R)**: Toolkit for analyzing single-cell chromatin data, such as ATAC-seq, integrated with Seurat workflows.  
  [Signac PBMC Vignette](https://stuartlab.org/signac/articles/pbmc_vignette.html)
  
- **DoubletFinder (R)**: Tool for detecting doublets (artificial merged cells) in single-cell RNA-seq data to improve data quality.  
  [DoubletFinder GitHub](https://github.com/chris-mcginnis-ucsf/DoubletFinder)

- **Seurat Integration (R)**: Workflow for integrating multiple single-cell datasets to correct batch effects and combine data across experiments.  
  [Seurat Integration Tutorial](https://satijalab.org/seurat/articles/seurat5_integration)

- **decoupleR (R)**: Framework that bundles multiple statistical methods to infer transcription-factor and pathway activities from omics data.
  [decoupleR GitHub](https://github.com/saezlab/decoupleR)

- **SComatic (Python)**: Functionalities to detect somatic single-nucleotide mutations in single-cell RNA-seq and single-cell ATAC-seq.
  [SComatic GitHub](https://github.com/cortes-ciriano-lab/SComatic)

- **InferCNV (R/Python)**: Tool to identify somatic large-scale chromosomal copy number alterations, such as gains or deletions of entire chromosomes or large segments of chromosomes.  
  [InferCNV GitHub](https://github.com/broadinstitute/inferCNV/wiki)

- **sccomp (R)**: Package for differential composition analysis in single-cell and compositional count data. It uses robust statistical models to identify changes in cell-type proportions between conditions, accounting for variability and overdispersion.  
  [sccomp GitHub](https://github.com/MangiolaLaboratory/sccomp) | [sccomp Bioconductor](https://www.bioconductor.org/packages/devel/bioc/vignettes/sccomp/inst/doc/introduction.html) | [sccomp explanation](https://www.sc-best-practices.org/conditions/compositional.html#why-cell-type-count-data-is-compositional)

- **Liana (R)**: Computational framework designed to analyze and infer cell-cell communication (CCC) by studying interactions between ligands and receptors.  
  [Liana Tutorial](https://saezlab.github.io/liana/)

- **CellChat (R)**: Computational framework designed for inferring, analyzing, and visualizing cell-cell communication.   
  [CellChat GitHub](https://github.com/jinworks/CellChat)

- **Nichenet (R)**:It predicts how ligands from sender cells influence gene expression in receiver cells by linking cell–cell communication signals to target gene regulation.  
  [Nichenet GitHub](https://github.com/saeyslab/nichenetr?tab=readme-ov-file)

---
# Spatial Transcriptomics Tools
 - **Squidpy (Python)** : Tutorials showcasing core Squidpy functionalities by applying them to a diverse set of different spatial datasets.  
  [Squidpy Documentation](https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/index.html)

- **OSTA (R)**: Online book Orchestrating Spatial Transcriptomics Analysis with Bioconductor.  
  [OSTA Tutorial](https://bioconductor.org/books/release/OSTA/)

### Image-based Data
- **Seurat Spatial (Image-based) (R)**: Comprehensive workflow for analyzing image-based spatial transcriptomics data using Seurat.  
  [Seurat Spatial Vignette](https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2)

### Sequencing-based Data
- **Seurat Spatial (Sequencing-based) (R)**: Tutorial for processing and analyzing sequencing-based spatial transcriptomics data.  
  [Seurat Spatial Sequencing Vignette](https://satijalab.org/seurat/articles/spatial_vignette)
  
- **Scanpy (Python)**: Comprehensive framework for the analysis of single-cell gene expression data.  
  [Scanpy Visium Tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/spatial/basic-analysis.html#visualization-in-spatial-coordinates)

### VISIUM HD
- **Seurat VISIUM HD (R)**: Analysis pipeline tailored for high-definition VISIUM spatial transcriptomics data.  
  [Seurat VISIUM HD Vignette](https://satijalab.org/seurat/articles/visiumhd_analysis_vignette)

### Other Tools
- **spacexr (R)**: Spatial cell type deconvolution using RCTD (Robust Cell Type Deconvolution). RCTD estimates the proportion of different cell types at each spatial location by combining spatial data with reference gene expression profiles of cell types.  
  [spacexr GitHub](https://github.com/dmcable/spacexr) | [spacexr Tutorial Bioconductor](https://www.bioconductor.org/packages/release/bioc/vignettes/spacexr/inst/doc/rctd-tutorial.html)

- **lisaClust (R)**: Bioconductor package for clustering and analyzing spatial transcriptomics data.  
  [lisaClust Bioconductor](https://www.bioconductor.org/packages/release/bioc/vignettes/lisaClust/inst/doc/lisaClust.html) | [lisaClust GitHub](https://github.com/SydneyBioX/lisaClust)

- **Banksy (R)**: Package for spatial transcriptomics analysis focusing on spatial autocorrelation and gene expression patterns.  
  [Banksy Documentation](https://prabhakarlab.github.io/Banksy/)

- **SPLIT (R)**: Spatial Purification of Layered Intracellular Transcripts in Xenium Data, a novel method that integrates scRNA-seq with RCTD deconvolution to enhance signal purity.  
[SPLIT Github](https://github.com/bdsc-tds/SPLIT?tab=readme-ov-file)

- **spicyR (R)**:  Provides a framework for performing inference on changes in spatial relationships between pairs of cell types for cell-resolution spatial omics technologiesell imaging data.  
  [spicyR Bioconductor](https://www.bioconductor.org/packages/release/bioc/vignettes/spicyR/inst/doc/spicyR.html)

- **sccomp (R)**: Package for differential composition analysis in single-cell and compositional count data. It uses robust statistical models to identify changes in cell-type proportions between conditions, accounting for variability and overdispersion.  
  [sccomp GitHub](https://github.com/MangiolaLaboratory/sccomp) | [sccomp Bioconductor](https://www.bioconductor.org/packages/devel/bioc/vignettes/sccomp/inst/doc/introduction.html) | [sccomp explanation](https://www.sc-best-practices.org/conditions/compositional.html#why-cell-type-count-data-is-compositional)

- **statial (R)**: Toolkit for spatially resolved transcriptomics data analysis, focusing on spatial neighborhood analysis and cell-cell interaction statistics.  
  [statial Bioconductor](https://www.bioconductor.org/packages/release/bioc/vignettes/Statial/inst/doc/Statial.html#region-analysis-using-lisaclust)

- **CRAWDAD (R)**: Statistical framework that uses cell-type labeled spatial omics data to identify the colocalization or separation of cell types at different length scales.  
  [CRAWDAD Github](https://github.com/JEFworks-Lab/CRAWDAD)

- **Scanpy (Python)**: Comprehensive framework for the analysis of single-cell gene expression data.  
  [Scanpy Tutorial](https://scanpy.readthedocs.io/en/1.10.x/tutorials/spatial/basic-analysis.html)

- **CellChat (R)**: Computational framework designed for inferring, analyzing, and visualizing cell-cell communication.   
  [CellChat GitHub](https://github.com/jinworks/CellChat)

- **SpaceTrooper (R):** Provides quality control for imaging-based spatial transcriptomics data (Xenium, Merfish/Merscope, CosMx), using GLM-based outlier detection, QC metrics on cell morphology and probe counts, and visualization of cell geometries with ggplot2.  
  [SpaceTrooper Bioconductor vignette](https://bioconductor.org/packages/devel/bioc/vignettes/SpaceTrooper/inst/doc/introduction.html)

- **DESpace (R)**: DESpace is a framework for identifying spatially variable genes (SVGs) and differential spatial variable pattern (DSP) genes.  
  [DESpace vignette](https://peicai.github.io/DESpace/articles/SVG.html)

- **SpatialQM (R)**: SpatialQM is a package that supports loading Spatial In-Situ datasets and calculating Quality Control metrics to aid understanding of data quality specifically tailored for In-situ spatial omics data (e.g CosMx, Xenium).  
  [SpatialQM Github](https://github.com/Center-for-Spatial-OMICs/SpatialQM)

---
  
# Mosaic Integration tools

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
  
# Protein Design Tools
- **Generative-protein-design-workshop (Python)**: Hands-on workshop introducing generative models for protein design. Provides tutorials and notebooks to build, train, and evaluate pipelines for creating novel protein sequences.  
  [Generative Protein Design Workshop · GitHub](https://github.com/romainstuder/evosite3d/tree/master/tutorials/generative-protein-design-workshop)

- **RFdiffusion (Python)**: Open source method for generating novel protein backbones using diffusion models, with options for tasks like motif scaffolding, binder design, symmetric oligomers, and macrocyclic peptides.  
  [RFdiffusion · GitHub](https://github.com/RosettaCommons/RFdiffusion)

- **ProteinMPNN (Python)**: Deep learning method that designs amino acid sequences predicted to fold into given protein backbone structures.  
  [ProteinMPNN · GitHub](https://github.com/dauparas/ProteinMPNN)

- **Boltz-2 (Python)**: Energy-based generative model for predicting and designing protein–ligand or protein–protein complexes, capable of estimating binding conformations and relative affinities by sampling from learned Boltzmann distributions.    
  [Boltz2 · GitHub](https://github.com/jwohlwend/boltz)

- **AlphaFold3 (Python)**: Geometry-driven deep learning system that predicts the 3D structures and interactions of proteins, nucleic acids, and small molecules with near-atomic accuracy, but does not estimate binding energies or affinities.  
  [AlpaFold3 · GitHub](https://github.com/google-deepmind/alphafold3)
- **BoltzGen (Python)**: BoltzGen is an all-atom generative model for designing and predicting 3D biomolecular structures and sequences.  
  [Boltzgen · GitHub](https://github.com/HannesStark/boltzgen)
