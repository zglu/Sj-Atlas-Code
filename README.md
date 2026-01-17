# Sj-Atlas-Code
This repository contains the code and analysis pipelines for the study: "A dynamic single-cell and spatial transcriptomic atlas of *Schistosoma japonicum* development and sexual maturation."

### 🌐 Interactive Online Resource
The data from this study is hosted as a publicly available, interactive resource. You can explore cell clusters, gene expression profiles, and spatial niches without any coding at: 👉 https://schisto.xyz/xxx

### 📖 Overview
Schistosomiasis remains a major global health challenge, primarily driven by the prolific egg production of sexually mature parasites. This project provides a foundational molecular framework for understanding schistosome biology through:

**Dynamic scRNA-seq**: Covering key stages of sexual maturation and egg production.

**Spatial Transcriptomics**: Mapping tissue-resolved cellular niches and the male–female pairing interface for the first time.

**Cross-species integration**: Comparative analysis with *S. mansoni* to identify conserved and species-divergent cellular plasticity.

### 📂 Repository Structure

```text
.
├── 10X_scRNA-seq/
│   ├── 01_preprocessing/           # QC, doublet removal, filtering, and normalization
│   ├── 02_integration_clustering/  # Integration, dimensionality reduction, clustering, and signature scores
│   ├── 03_annotation_markers/      # Cell type annotation and marker gene identification
│   ├── 04_sub-clustering/          # High-resolution analysis of specific lineages (e.g., germlines)
│   └── 05_cross-species/           # Integration with S. mansoni data (Wendt et al., 2020)
├── Stereo-seq_spatial/             # Processing, tissue mapping, and niche analysis for spatial data
├── LICENSE                         # MIT License
└── README.md
```

### 📜 Citation
To be added.

*Scripts subjected to update.*
