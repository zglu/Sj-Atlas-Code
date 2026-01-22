## aggregate expression on main cell type to create pseudo-bulk samples for each stage/sex
merged<-readRDS("all_integrated_SCT-RPCA_RNA.rds")
pseudo <- AggregateExpression(merged, assays = "RNA", return.seurat = T, group.by = c("orig.ident", "main_type")) # main_type

												 orig.ident    cell_type
Female16_ambiguous       Female16_ambiguous    ambiguous
Female16_flame               Female16_flame        flame
Female16_GSC                   Female16_GSC          GSC

pseudo$sample<-gsub("(.+)*[_].*", "\\1", pseudo$orig.ident)
pseudo$dpi<-paste0("D", as.numeric(gsub("\\D", "", pseudo$sample)))
Idents(pseudo)<-"sample"
saveRDS(pseudo, file="all_integrated_main_type_pseudo.rds")

## DE analysis between samples
# DESeq2
library(DESeq2)
counts<-pseudo@assays$RNA@layers$counts
rownames(counts)<-rownames(pseudo)
colnames(counts)<-colnames(pseudo)
groups<-as.factor(pseudo$sample)  # or use pseudo$dpi

coldata<-data.frame(row.names=colnames(counts), groups)
dds<-DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~0+groups)
dds<-DESeq(dds)
ResultNames(dds)
[1] "groupsFemale16" "groupsFemale20" "groupsFemale26" "groupsMale16"   "groupsMale20"   "groupsMale26"

save(dds, file="pseudo_main_type_Sample.dds.rda")



## DE results; PCA
dge1<-results(dds, cooksCutoff=TRUE, contrast=c("groups","Female20","Female16"))
dge1<-as.data.frame(dge1)
dge1$gene<-rownames(dge1)
# volcano
dge1$threshold<-ifelse(dge1$log2FoldChange>=1 & dge1$padj<0.01,"UP", ifelse(dge1$log2FoldChange<=-1 & dge1$padj<0.01, "DOWN", "NO"))

# volcano plot with colors and labels
dge1$delabel <- dge1$gene
#dge1$delabel[dge1$threshold=="NO"]<-NA
dge1$delabel[abs(dge1$log2FoldChange)<5 | dge1$padj>0.01]<-NA
ggplot(as.data.frame(dge1), aes(x=log2FoldChange, y=-log10(padj), color=threshold, label=delabel))+geom_point()+scale_colour_manual(values=c("UP"="red","DOWN"="blue", "NO"="grey"))+geom_text(size=3, nudge_y = 0.3)+labs(title="Female20 vs Female16")

### pca
myCol<-unique(c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"), brewer.pal(12, "Set3"), brewer.pal(8, "Set1"), brewer.pal(8, "Accent")))

vst<-vst(dds)
plotPCA(vst, intgroup=c("groups"))+scale_color_manual(values = myCol)+geom_text(label=pseudo$cell_type, size=2.2, col="black", nudge_y=1)


