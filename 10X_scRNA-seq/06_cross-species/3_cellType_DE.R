library(Seurat)
library(DESeq2)

## pseudo bulk

merged<-readRDS("merged_SmSj_ortho.rds_seu5log80harmony.rds")
pseudo <- AggregateExpression(merged, assays = "RNA", return.seurat = T, group.by = c("species", "orig.anno")) # main_type

mymeta<-merged@meta.data
mymeta$orig.anno<-gsub("_", "-", mymeta$orig.anno)
meta1<-mymeta[,c("species", "orig.anno", "main_type")]
meta1$spcl<-paste0(meta1$species, "_", meta1$orig.anno)
write.table(meta1, file="merged_SmSj_ortho_main_type_SpCl.txt", quote=F, row.names = F, col.names = T, sep="\t")
## vim print $2 and $2 and uniq

meta2<-read.csv("merged_SmSj_ortho_main_type_SpCl.csv", header=T)
pseudo$main_type<-meta2$main_type[match(rownames(pseudo@meta.data), meta2$spcl)]
dim(pseudo@meta.data) # in total 130 original anno for clusters
[1] 130   4
saveRDS(pseudo, file="Pseudo_species-anno.rds")

# DESeq2
pseudo<-readRDS("Pseudo_species-anno.rds")

counts<-pseudo@assays$RNA@layers$counts
rownames(counts)<-rownames(pseudo)
colnames(counts)<-colnames(pseudo)

species<-as.factor(pseudo$species)
maintype<-as.factor(pseudo$main_type)
spcl<-as.factor(paste0(pseudo$species, "_", pseudo$main_type))

coldata<-data.frame(row.names=colnames(counts), species, maintype, spcl) # can add multiple factors

                species maintype        spcl
Sj_C0-neoblast       Sj neoblast Sj_neoblast
Sj_C1-muscle         Sj   muscle   Sj_muscle
Sj_C10-neoblast      Sj neoblast Sj_neoblast
Sj_C11-muscle        Sj   muscle   Sj_muscle

dds<-DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~spcl)
dds<-DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~0+spcl) # to avaoid setting intercept
# design =~0+maintype
dds<-DESeq(dds)
save(dds, file="pseudo_species-anno_spcl.dds.rda")


## calculate DEG between different groups; ie S1, neoblast, tegument etc
library(DESeq2)
load("pseudo_species-anno_spcl.dds.rda")

dge1<-results(dds, cooksCutoff=TRUE, contrast=c("spcl","Sj_s1","Sm_S1")) ## here change comparison to S1
dge1<-as.data.frame(dge1)
dge1$gene<-rownames(dge1)
# volcano
dge1$threshold<-ifelse(dge1$log2FoldChange>=1 & dge1$padj<0.01,"UP", ifelse(dge1$log2FoldChange<=-1 & dge1$padj<0.01, "DOWN", "NO"))

table(dge1$threshold)

dge1$delabel <- dge1$gene
#dge1$delabel[dge1$threshold=="NO"]<-NA
dge1$delabel[abs(dge1$log2FoldChange)<5 | dge1$padj>0.01]<-NA

# combine prod and kegg info
dge1$Sj_gene<-gsub("Sm.*-Sj", "EWB00_", dge1$gene)
geneprod<-read.delim("sjscsp/HuSjv2_prod.txt", sep="\t", header=T)
keggpath<-read.delim("sjscsp/SjID_Knames.txt", sep="\t", header=F)
SjClusters<-read.delim("sjscsp/Sj_gene-cell_type.txt", sep=" ", header=F)
deg_comb<-merge(dge1, geneprod, all.x=TRUE, by.x="Sj_gene", by.y="X.ID")
deg_comb<-merge(deg_comb, keggpath, all.x=TRUE, by.x="Sj_gene", by.y="V1")
deg_comb<-merge(deg_comb, SjClusters, all.x=TRUE, by.x="Sj_gene", by.y="V1")
deg_comb<-deg_comb[order(-deg_comb$log2FoldChange, deg_comb$padj),]

colnames(deg_comb)[12:13]<-c("KEGG_name","Cluster_marker")

write.csv(deg_comb, file="pseudo_s1_Sj-vs-Sm.csv")

## remove NA for plotting
dge1<-subset(dge1, subset=!is.na(dge1$threshold))

## modify labels to keep Sj for UP and Sm for down
dge1$delabel[dge1$threshold=="UP"]<-gsub(".*[-](.+)", "\\1", dge1$delabel[dge1$threshold=="UP"])
dge1$delabel[dge1$threshold=="DOWN"]<-gsub("(.+)*[-].*", "\\1", dge1$delabel[dge1$threshold=="DOWN"])


p<-ggplot(as.data.frame(dge1), aes(x=log2FoldChange, y=-log10(padj), color=threshold))+geom_point()+
  scale_colour_manual(values=c("UP"="#1f77b4","DOWN"="#ff7f0e", "NO"="grey"))+labs(title="Sj s1 vs Sm S1")#+geom_text(size=3, nudge_y = 0.3)

pdf("Sj-Sm_pseudo_volcano_s1.pdf", height=8, width=10)
p+geom_text_repel(aes(label=dge1$delabel))+theme(legend.title=element_blank())
dev.off()

