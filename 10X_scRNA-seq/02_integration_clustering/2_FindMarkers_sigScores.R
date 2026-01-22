#Rscript 2_FindMarker_sigScores.R [SCT rds]

args<-commandArgs(T)


suppressPackageStartupMessages({
  library(RColorBrewer)
  library(Seurat)
  library(gridExtra)
  library(ggplot2)
  library(SCINA)
  library(AUCell)
  library(ggsci)
  library(UCell)
  library(dplyr)
  library(SingleCellExperiment)
})


myCol<-unique(c(pal_d3("category10")(10),pal_rickandmorty("schwifty")(12), pal_lancet("lanonc")(9),
                pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),
                pal_jama("default")(7),pal_jco("default")(10),
                pal_locuszoom("default")(7),pal_startrek("uniform")(7),
                pal_tron("legacy")(7),pal_futurama("planetexpress")(12),
                pal_simpsons("springfield")(16),
                pal_gsea("default")(12)))


data<-readRDS(args[1])

table(data@active.ident)

data<-FindClusters(res=2)

## dim and feature plot
pdf(paste0(args[1],"_dimplot.pdf"), width=18, height=8)
DimPlot(data, reduction = "umap",label=TRUE, label.color='black', cols=myCol)+NoLegend()+coord_fixed()
DimPlot(data, reduction = "umap",label=TRUE, label.color='black',split.by = "orig.ident", cols=myCol)+NoLegend()+coord_fixed()
dev.off()


DefaultAssay(data) <- "RNA"
message("Calculating scores using addModuleScore:")

allgenes<-rownames(data)

markers<-read.csv("allcombined_markers_1e-20_AUC0.7_grouped.csv", header=T)

p1 <- lapply(colnames(markers), function(i) {
  m <- markers[, i]
  m <- intersect(m, allgenes)
  data<-AddModuleScore(data, features=list(m), name=i)
  FeaturePlot(data, features = paste0(i, "1"), label = TRUE, repel = TRUE, raster=TRUE) + coord_fixed()+theme_void()+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlOrRd"))+labs(color="SeuratScore")
})

ggsave(
  filename = paste0(args[1], "_SeuratScores.pdf"),
  plot = marrangeGrob(p1, nrow=2, ncol=4),
  width = 18, height = 8
)

## markers
deg_roc <- FindAllMarkers(data,only.pos = TRUE, test.use="roc", verbose = F) # method: roc
table(deg_roc$cluster)
rownames(deg_roc)<-gsub("EWB00-", "EWB00_", rownames(deg_roc))
deg_roc$gene<-gsub("EWB00-", "EWB00_",deg_roc$gene)

deg_roc$id<-rownames(deg_roc)

#write.csv(deg_roc, file=paste0(args[1], "_allMarkers_roc.csv"))

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

library(openxlsx)
write.xlsx(deg_roc_comb, paste0(args[1], "_allMarkers_roc.xlsx"))

