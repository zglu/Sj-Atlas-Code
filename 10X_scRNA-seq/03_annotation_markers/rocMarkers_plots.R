# Rscript rocMarkers.R all_integrated_SCT-RPCA_RNA.rds cell_type

args<-commandArgs(T)
if (length(args)<2) {
  stop("Usage: Rscript rocMarkers.R [rds] [group]")
}

suppressPackageStartupMessages({
library(Seurat)
library(dplyr)
library(ggplot2)
library(SingleCellExperiment)
library(RColorBrewer)
library(gridExtra)
})

myCol<-unique(c(brewer.pal(8, "Set1"), brewer.pal(12, "Set3"),brewer.pal(12, "Paired"), brewer.pal(8, "Accent"), brewer.pal(8, "Dark2"), brewer.pal(8, "Set2")))

wd<-"/data/work/test/"
seuData<-paste0(wd, args[1])
myGrp<-args[2] #cell_type, main_type, anno

merged<-readRDS(seuData)
merged

Idents(merged)<-myGrp

## markers
deg_roc <- FindAllMarkers(merged,only.pos = TRUE, test.use="roc", assay="RNA") 
table(deg_roc$cluster)
rownames(deg_roc)<-gsub("EWB00-", "EWB00_", rownames(deg_roc))
deg_roc$gene<-gsub("EWB00-", "EWB00_",deg_roc$gene)

deg_roc$id<-rownames(deg_roc)

write.csv(deg_roc, file=paste0(seuData, "_allMarkers_roc.", myGrp, ".csv"))

geneprod<-read.delim("HuSjv2_prod.txt", sep="\t", header=T)
SmOrth<-read.delim("SmSj-ortho.txt", sep="\t", header=F)
SmClusters<-read.delim("Sj_marker_wendt2020.txt", sep=" ", header=F)
library(data.table)
#setDT(deg_roc, keep.rownames="gene")
deg_roc_comb<-merge(deg_roc, SmOrth, all.x=TRUE, by.x="gene", by.y="V2")
deg_roc_comb<-merge(deg_roc_comb, SmClusters, all.x=TRUE, by.x="gene", by.y="V1")
deg_roc_comb<-merge(deg_roc_comb, geneprod, all.x=TRUE, by.x="gene", by.y="X.ID")

deg_roc_comb<-deg_roc_comb[order(deg_roc_comb$cluster, -deg_roc_comb$myAUC),]
colnames(deg_roc_comb)[10:11]<-c("Sm_ortholog", "Wendt_clusters")

write.csv(deg_roc_comb, file=paste0(seuData, "_combMarkers_roc.", myGrp, ".csv"))

topmarkers<-deg_roc%>% group_by(cluster) %>% top_n(n = 3, wt = myAUC) # n=5 for certain cell types
topfeatures<-as.vector(topmarkers$gene)
topfeatures<-unique(topfeatures)
topfeatures<-gsub("EWB00_", "EWB00-", topfeatures)

pdf(paste0(seuData, "_topMarkers.", myGrp, ".dot.pdf"), height=9, width=21)
DotPlot(merged, features = topfeatures) + RotatedAxis()
dev.off()

pdf(paste0(seuData, "_topMarkers.", myGrp, ".heat.pdf"), height=12, width=15)
DoHeatmap(merged, features=topfeatures)+NoLegend() 
dev.off()
