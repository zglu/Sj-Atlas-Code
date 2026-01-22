library(Seurat)
library(dplyr)
library(ggplot2)

Female16<-Read10X("F16/filtered_feature_bc_matrix/", gene.column=1) # default is col2: gene names. col1 is gene ids
Female20<-Read10X("F20/filtered_feature_bc_matrix/", gene.column=1)
Female26<-Read10X("F26/filtered_feature_bc_matrix/", gene.column=1)
Male16<-Read10X("M16/filtered_feature_bc_matrix/", gene.column=1)
Male20<-Read10X("M20/filtered_feature_bc_matrix/", gene.column=1)
Male26<-Read10X("M26/filtered_feature_bc_matrix/", gene.column=1)

Female16<-CreateSeuratObject(Female16, project = "Female16", min.cells = 5, min.features = 500) 
Female20<-CreateSeuratObject(Female20, project = "Female20", min.cells = 5, min.features = 500)
Female26<-CreateSeuratObject(Female26, project = "Female26", min.cells = 5, min.features = 500)
Male16<-CreateSeuratObject(Male16, project = "Male16", min.cells = 5, min.features = 500)
Male20<-CreateSeuratObject(Male20, project = "Male20", min.cells = 5, min.features = 500)
Male26<-CreateSeuratObject(Male26, project = "Male26", min.cells = 5, min.features = 500)


Female16[["percent.mt"]] <- PercentageFeatureSet(object = Female16, pattern = "^Sj-")
Female20[["percent.mt"]] <- PercentageFeatureSet(object = Female20, pattern = "^Sj-")
Female26[["percent.mt"]] <- PercentageFeatureSet(object = Female26, pattern = "^Sj-")
Male16[["percent.mt"]] <- PercentageFeatureSet(object = Male16, pattern = "^Sj-")
Male20[["percent.mt"]] <- PercentageFeatureSet(object = Male20, pattern = "^Sj-")
Male26[["percent.mt"]] <- PercentageFeatureSet(object = Male26, pattern = "^Sj-")

Female16<-subset(Female16, subset=percent.mt<5)
Female20<-subset(Female20, subset=percent.mt<5)
Female26<-subset(Female26, subset=percent.mt<5)
Male16<-subset(Male16, subset=percent.mt<5)
Male20<-subset(Male20, subset=percent.mt<5)
Male26<-subset(Male26, subset=percent.mt<5)

# doublet detection for each sample
suppressPackageStartupMessages({
  library(scDblFinder)
})

F16.sce<-as.SingleCellExperiment(Female16)
F16.sce <- scDblFinder(F16.sce)

Female16@meta.data<-as.data.frame(colData(F16.sce))
Female16[["ident"]]<-NULL
saveRDS(Female20, file="F16_mt5_gene500_raw.rds")

### same for the other samples
