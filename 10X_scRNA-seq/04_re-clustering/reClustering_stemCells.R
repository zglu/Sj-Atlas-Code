library(Seurat)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(gridExtra)

myCol<-unique(c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"), brewer.pal(12, "Set3"), brewer.pal(8, "Set1"), brewer.pal(8, "Accent"))) 

## get stem cells
SeuObj$ago2_data<-SeuObj@assays$RNA@data["EWB00-004346",]
FeaturePlot(SeuObj, features = "ago2_data", coord.fixed = T, order=T)+ggtitle("ago2-1")

stem<-subset(SeuObj, subset=ago2_data>1) # stems with ago2-1 expr data > 1

SeuObj$stem<-SeuObj$ago2_data
SeuObj$stem[SeuObj$stem>1]<-"Stem cells"
SeuObj$stem[SeuObj$stem<=1]<-"Other cells"
DimPlot(SeuObj, group.by="stem")+coord_fixed()+labs(title="")


### clustering of different stages
stem$orig_clusters<-stem$seurat_clusters

phase<- "16" # "20" and "26"

stem<-stem[,grepl(phase, stem$orig.ident)]


### sct
stem %>% SCTransform(verbose = FALSE,variable.features.n = 3000) %>%
  RunPCA(verbose = FALSE,assay="SCT") %>%
  RunUMAP( dims = 1:20, verbose = FALSE)%>%
  FindNeighbors( dims = 1:20, verbose = FALSE)%>%
  FindClusters(res=0.3, verbose = FALSE) -> stem

saveRDS(stem, paste0("Sj_stem_D", phase, "_sct.clustering.rds"))

p1<-DimPlot(stem, label=T, repel = T, cols=myCol)+coord_fixed()+NoLegend()+labs(title=paste0("Clustering of D", phase, " ago2-1+ cells"))
p2<-DimPlot(stem, group.by="cell_type", label=F, repel = T, cols=myCol)+coord_fixed()+labs(title=paste0("Original cell type in D", phase, " ago2-1+ cells"))
#DimPlot(stem, group.by="orig.ident", label=F, repel = T, cols=myCol)+coord_fixed()+labs(title="Clustering of ago2-1+ cells")
p3<-FeaturePlot(stem, features = "EWB00-002258", order=T, coord.fixed = T, label=T)+labs(title="nanos1 [2258]")


deg_roc <- FindAllMarkers(stem,only.pos = TRUE, test.use="roc", assay = "RNA") # method: roc
sjtf<-read.table("SjTF.ids", header=F)
#deg_roc$gene<-rownames(deg_roc)
deg_roc$gene<-gsub("EWB00-", "EWB00_", deg_roc$gene)
deg_roc$is.TF<-deg_roc$gene %in% sjtf$V1
geneprod<-read.delim("HuSjv2_prod.txt", sep="\t", header=T)
SmOrth<-read.delim("SjSm-ortho_comb.txt", sep=" ", header=F)
library(data.table)
dyna_comb<-merge(deg_roc, SmOrth, all.x=TRUE, by.x="gene", by.y="V1")
dyna_comb<-merge(dyna_comb, geneprod, all.x=TRUE, by.x="gene", by.y="X.ID")
dyna_comb<-dyna_comb[order(dyna_comb$cluster, -dyna_comb$myAUC),]
library(openxlsx)
write.xlsx(dyna_comb, file=paste0("Sj_stem_D", phase, "_clustering_markers.xlsx"))


## heatmap
deg_roc$gene<-gsub("EWB00_", "EWB00-", deg_roc$gene)
deg.top <- deg_roc %>% group_by(cluster) %>% top_n(n = 10, wt = myAUC)
deg.topmarkers<-as.vector(deg.top$gene)
deg.topmarkers<-unique(deg.topmarkers)

p4<-pheat <- DoHeatmap(stem,features = deg.top$gene)

pdf(paste0("Sj_stem_D", phase, "_clustering.pdf"), width=12, height=8)
grid.arrange(p1, p2, ncol=2)
grid.arrange(p3, p4, ncol=2)
dev.off()


p5<-pheat <- DoHeatmap(stem,features = deg_roc$gene)
pdf(paste0("Sj_stem_D", phase, "_clustering_heatmap_all.pdf"), width=8, height=12)
print(p5)
dev.off()
