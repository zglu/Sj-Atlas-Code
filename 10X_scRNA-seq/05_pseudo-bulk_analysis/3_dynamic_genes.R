### single comparisions
F2016<-results(dds, cooksCutoff=TRUE, independentFiltering=TRUE, contrast=c("groups","Female20","Female16")); 
F2016<-subset(F2016, subset=abs(log2FoldChange>=0.585) & padj<0.05) #190
F2620<-results(dds, cooksCutoff=TRUE, independentFiltering=TRUE, contrast=c("groups","Female26","Female20")); 
F2620<-subset(F2620, subset=abs(log2FoldChange>=0.585) & padj<0.05) #881
F2616<-results(dds, cooksCutoff=TRUE, independentFiltering=TRUE, contrast=c("groups","Female26","Female16")); 
F2616<-subset(F2616, subset=abs(log2FoldChange>=0.585) & padj<0.05) #1450
F2016$comp<-"Female20_vs_Female16"; F2016$gene<-rownames(F2016)
F2620$comp<-"Female26_vs_Female20"; F2620$gene<-rownames(F2620)
F2616$comp<-"Female26_vs_Female16"; F2616$gene<-rownames(F2616)

M2016<-results(dds, cooksCutoff=TRUE, independentFiltering=TRUE, contrast=c("groups","Male20","Male16")); 
M2016<-subset(M2016, subset=abs(log2FoldChange>=0.585) & padj<0.05) #289
M2620<-results(dds, cooksCutoff=TRUE, independentFiltering=TRUE, contrast=c("groups","Male26","Male20")); 
M2620<-subset(M2620, subset=abs(log2FoldChange>=0.585) & padj<0.05) #699
M2616<-results(dds, cooksCutoff=TRUE, independentFiltering=TRUE, contrast=c("groups","Male26","Male16")); 
M2616<-subset(M2616, subset=abs(log2FoldChange>=0.585) & padj<0.05) # 1223
M2016$comp<-"Male20_vs_Male16"; M2016$gene<-rownames(M2016)
M2620$comp<-"Male26_vs_Male20"; M2620$gene<-rownames(M2620)
M2616$comp<-"Male26_vs_Male16"; M2616$gene<-rownames(M2616)

alldegs<-rbind(F2016, F2620, F2616, M2016, M2620, M2616)
alldegs<-as.data.frame(alldegs)
write.csv(alldegs, file="pseudo_main_type_Sample_allDEGs005-1.5.csv")

### add comparison to cols
test<- alldegs[,c("log2FoldChange","comp", "gene")] %>% group_by(comp)
library(reshape2)
casted<-dcast(test, comp ~ gene, value.var="log2FoldChange")
comb<-t(casted)
comb[is.na(comb)]<-0
comb<-as.data.frame(comb)
write.csv(comb, file="pseudo_main_type_Sample_allDEGs_comb.csv")
# 2670 genes # remove header

### combine with product, TF, and cluster marker information
dyna<-read.csv("pseudo_main_type_Sample_allDEGs_comb_dynamic.csv", row.names = "gene")
sjtf<-read.table("SjTF.ids", header=F)
dyna$is.TF<-rownames(dyna) %in% sjtf$V1
dyna$gene<-rownames(dyna)
geneprod<-read.delim("HuSjv2_prod.txt", sep="\t", header=T)
SmOrth<-read.delim("SmSj-ortho.txt", sep="\t", header=F)
SjClusters<-read.delim("Sj_gene-cell_type.txt", sep=" ", header=F)
library(data.table)
dyna_comb<-merge(dyna, SmOrth, all.x=TRUE, by.x="gene", by.y="V2")
dyna_comb<-merge(dyna_comb, geneprod, all.x=TRUE, by.x="gene", by.y="X.ID")
dyna_comb<-merge(dyna_comb, SjClusters, all.x=TRUE, by.x="gene", by.y="V1")
colnames(dyna_comb)[9:11]<-c("Sm_ortholog", "description", "Sj_cluster_marker")
library(openxlsx)
write.xlsx(dyna_comb, file="pseudo_main_type_Sample_allDEGs_comb_dynamic.xlsx")


### plot log2FC values in different comparisons
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(dyna), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(dyna)/paletteLength, max(dyna), length.out=floor(paletteLength/2)))
pheatmap(dyna, scale = "none", color=myColor, breaks=myBreaks, show_rownames  = F, main="dynamic genes in comparison")

### how many are TF
dynatf<-subset(dyna, subset=rownames(dyna) %in% sjtf$V1)
pheatmap(dynatf, scale = "none", color=myColor, file="dynamic_TFs_log2FC.pdf", cellheight = 12, cellwidth = 12, breaks=myBreaks, show_rownames  = T, main="dynamic TFs in comparison")

### clustered dot plot
library(Seurat)
library(scCustomize)
pseudo<-readRDS("all_integrated_dpi-sex_pseudo.rds")
pdf("dynamic_TF_dot.pdf", height=12, width=5.5)
Clustered_DotPlot(pseudo, features = gsub("EWB00_", "EWB00-", rownames(dynatf)), group.by="orig.ident", row_label_size = 8, column_label_size = 10)
dev.off()

pdf("dynamic_TF_mainType_dot.pdf", height=12, width=10)
Clustered_DotPlot(SeuObj, features = gsub("EWB00_", "EWB00-", rownames(dynatf)), group.by="main_type", row_label_size = 8, column_label_size = 10)
dev.off()
