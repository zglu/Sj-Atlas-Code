# Sj-Atlas-Code
This repository contains the code and analysis pipelines for the study: Z. Lu, X. Wang, S. Li, et al. “Single-Cell and Spatial Transcriptomics Unveil Key Regulators Governing Cell Differentiation for *Schistosoma japonicum* Sexual Development.” Advanced Science (2026): e76329. https://doi.org/10.1002/advs.76329

### 🌐 Interactive Online Resource
The data from this study is hosted as a publicly available, interactive resource. You can explore cell clusters, gene expression profiles, and spatial niches without any coding at: 👉 https://schisto.xyz/sj-atlas

### 📖 Overview
Schistosomiasis remains a major global health challenge, primarily driven by the prolific egg production of sexually mature parasites. This project provides a foundational molecular framework for understanding schistosome biology through:

**Dynamic scRNA-seq**: Covering key stages of sexual maturation and egg production (16, 20, and 26 dpi).

**Spatial Transcriptomics**: Mapping tissue-resolved cellular niches and the male–female pairing interface for the first time.

**Cross-species integration**: Comparative analysis with *S. mansoni* scRNA-seq data (Wendt et al 2020) to identify conserved and species-divergent cellular plasticity.

### 📂 Repository Structure

```text
.
├── 10X_scRNA-seq/
│   ├── 01_preprocessing/           # Data QC, doublet removal, filtering, and normalization
│   ├── 02_integration_clustering/  # RPCA Integration, dimensionality reduction, clustering, and signature scores
│   ├── 03_annotation_markers/      # Cell type annotation and marker gene identification (all, or certain cell types)
│   ├── 04_re-clustering/           # High-resolution re-clustering of stem cells
│   ├── 05_pseudo-bulk_analysis     # Pseudo-bulk differential expression and dynamic genes
│   └── 06_cross-species/           # Integration with S. mansoni data (Wendt et al., 2020)
├── Stereo-seq_spatial/             # Processing, tissue mapping, and niche analysis for spatial data
├── LICENSE                         # MIT License
└── README.md
```


