# Rscript singleComp_combInfo.R Male26 Male20

args<-commandArgs(T) 

grp1<-args[1]
grp2<-args[2]

library(DESeq2)
library(ggplot2)
load("pseudo_main_type_Sample.dds.rda")
dge1<-results(dds, cooksCutoff=TRUE, contrast=c("groups", grp1, grp2))
dge1<-as.data.frame(dge1)
dge1$gene<-rownames(dge1)
# volcano
dge1$threshold<-ifelse(dge1$log2FoldChange>=1 & dge1$padj<0.01,"UP", ifelse(dge1$log2FoldChange<=-1 & dge1$padj<0.01, "DOWN", "NO"))

# volcano plot with colors and labels
dge1$delabel <- dge1$gene
#dge1$delabel[dge1$threshold=="NO"]<-NA
dge1$delabel[abs(dge1$log2FoldChange)<5 | dge1$padj>0.01]<-NA

pdf(paste0(grp1, "_", grp2, "_volcano.pdf"), width=9, height=6)

ggplot(as.data.frame(dge1), aes(x=log2FoldChange, y=-log10(padj), color=threshold, label=delabel))+
  geom_point()+scale_colour_manual(values=c("UP"="red","DOWN"="blue", "NO"="grey"))+
  geom_text(size=3, nudge_y = 0.3)+labs(title="Female26 vs Female20")

ggplot(as.data.frame(dge1), aes(x=log2FoldChange, y=-log10(padj), color=threshold))+
  geom_point()+scale_colour_manual(values=c("UP"="red","DOWN"="blue", "NO"="grey"))+
  labs(title=paste0(grp1, " vs ", grp2))
dev.off()

table(dge1$threshold)
sjtf<-read.table("SjTF.ids", header=F)
dge1$gene<-rownames(dge1)
dge1$gene<-gsub("EWB00-", "EWB00_", dge1$gene)
dge1$is.TF<-dge1$gene %in% sjtf$V1
geneprod<-read.delim("HuSjv2_prod.txt", sep="\t", header=T)
SmOrth<-read.delim("SjSm-ortho_comb.txt", sep=" ", header=F)
SjClusters<-read.delim("Sj_gene-cell_type.txt", sep=" ", header=F)
library(data.table)
dyna_comb<-merge(dge1, SmOrth, all.x=TRUE, by.x="gene", by.y="V1")
dyna_comb<-merge(dyna_comb, geneprod, all.x=TRUE, by.x="gene", by.y="X.ID")
dyna_comb<-merge(dyna_comb, SjClusters, all.x=TRUE, by.x="gene", by.y="V1")
colnames(dyna_comb)[11:13]<-c("Sm_ortholog", "description", "Sj_cluster_marker")
dyna_comb$delabel<-NULL
library(openxlsx)
write.xlsx(dyna_comb, file=paste0(grp1, "_vs_", grp2, "_deg.xlsx"))

