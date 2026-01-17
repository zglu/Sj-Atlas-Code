# Sj-Atlas-Code
This repository contains the code and analysis pipelines for the study: "A dynamic single-cell and spatial transcriptomic atlas of *Schistosoma japonicum* development and sexual maturation."

🌐 Interactive Online Resource
The data from this study is hosted as a publicly available, interactive resource. You can explore cell clusters, gene expression profiles, and spatial niches without any coding at: 👉 https://schisto.xyz/xxx

📖 Overview
Schistosomiasis remains a major global health challenge, primarily driven by the prolific egg production of sexually mature parasites. This project provides a foundational molecular framework for understanding schistosome biology through:

**Dynamic scRNA-seq**: Covering key stages of sexual maturation and egg production.
**Spatial Transcriptomics**: Mapping tissue-resolved cellular niches and the male–female pairing interface for the first time.
**Cross-species integration**: Comparative analysis with *S. mansoni* to identify conserved and species-divergent cellular plasticity.

📂 Repository Structure

```text
├── 10X_scRNA-seq/
│   ├── 01_preprocessing/  # QC, doublet removal, filtering, normalization, etc
│   ├── 02_integration_clustering/     # Integration, dimensionality reduction, clustering, and signatures
│   ├── 03_annotation_markers/        # Spatial transcriptomics mapping and niche analysis
│   ├── 04_sub-clustering/   # Sub-clustering for specific cell types
│   └── 05_cross-species_integration/    # Cross-species data integration with S. mansoni scRNA-seq data (Wendt et al 2020)
├── Stereo-seq_spatial/            # Data processing and annotation for Stereo-seq spatial transcriptomics data.
└── README.md

📜 Citation
To be added.

*Scripts subjected to update*
