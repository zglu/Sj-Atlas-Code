# Rscript 1_integrate_SCT-RPCA_all.R [F16 F20 F26 M16 M20 M26]

library(Seurat)
library(dplyr)
args<-commandArgs(T)

# each data has been normalised using sctransform: Rscript 1_single-sct-clustering.R
Data1<-readRDS(paste0(args[1], "_Fil_SCT_PC50.rds"))
Data2<-readRDS(paste0(args[2], "_Fil_SCT_PC50.rds"))
Data3<-readRDS(paste0(args[3], "_Fil_SCT_PC50.rds"))
Data4<-readRDS(paste0(args[4], "_Fil_SCT_PC50.rds"))
Data5<-readRDS(paste0(args[5], "_Fil_SCT_PC50.rds"))
Data6<-readRDS(paste0(args[6], "_Fil_SCT_PC50.rds"))


Data1[["ident"]]<-NULL
Data2[["ident"]]<-NULL
Data3[["ident"]]<-NULL
Data4[["ident"]]<-NULL
Data5[["ident"]]<-NULL
Data6[["ident"]]<-NULL


int_list<-list(Data1, Data2, Data3, Data4, Data5, Data6)
features <- SelectIntegrationFeatures(object.list = int_list, nfeatures = 3000)

int_list <- PrepSCTIntegration(object.list = int_list, anchor.features = features)

## integration using RPCA
int_list <- lapply(X = int_list, FUN = RunPCA, features = features)

int.anchors <- FindIntegrationAnchors(object.list = int_list, normalization.method = "SCT",
    anchor.features = features, reduction = "rpca")
int.combined.sct <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT")

message("running PCA on integrated data...")
int.combined.sct <- RunPCA(int.combined.sct, verbose = FALSE)
int.combined.sct <- RunUMAP(int.combined.sct, reduction = "pca", dims = 1:50)
int.combined.sct <- FindNeighbors(int.combined.sct, reduction = "pca", dims = 1:50)
int.combined.sct <- FindClusters(int.combined.sct, resolution = 2)

message("Saving RDS file...")
saveRDS(int.combined.sct, file=paste0("all","_integrated_SCT-RPCA.rds"))


