library(anndata)
library(Seurat)
ann<-read_h5ad("G214.adjusted.cellbin.h5ad")
SeuObj <- CreateSeuratObject(counts = t(as.matrix(ann$Raw$X)), meta.data = ann$obs)

colnames(SeuObj@meta.data)[9:10]<-c("coord_x", "coord_y")
SeuObj$coord_y<--1*SeuObj$coord_y

## add spatial embedding
sp.embedding<-SeuObj@meta.data[,c("coord_x", "coord_y")]
colnames(sp.embedding)<-c("Spatial_1", "Spatial_2") # prefix with number
sp.embedding<-as.matrix(sp.embedding) # needs to be a matrix
SeuObj$spatial<-CreateDimReducObject(embeddings=sp.embedding, key='Spatial_', assay='RNA')


DimPlot(SeuObj, reduction="spatial", group.by = "orig.ident")
SeuObj<-NormalizeData(SeuObj)
SeuObj[["percent.mt"]] <- PercentageFeatureSet(object = SeuObj, pattern = "^Sj-")

saveRDS(SeuObj, file="G214.adjusted.cellbin.rds")
