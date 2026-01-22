## Create Sm ortholog-based object
# read Wendt et al data
Sm<-readRDS("Wendt2020_SmSC_RNA.rds")

# slice count matrix with 1:1 orthologs
SmCounts<-Sm@assays$RNA@counts
allSm<-rownames(Sm) # all genes detected in Sm
ortholist<-read.table("SmSj_1to1.txt", header=F, sep="\t")

# genes in both Sm data and 1:1 ortholog list
uniques<-intersect(allSm, ortholist$V1) 
SmCounts<-SmCounts[uniques, ]
## change gene names
rownames(SmCounts)<-ortholist$V3[match(rownames(SmCounts), ortholist$V1)]

# modify meta data
smmeta<-Sm@meta.data
smmeta$orig.ident<-smmeta$Group
smmeta<-smmeta[,c(1,12,11)] # changed "anno" to "cell_type"
colnames(smmeta)<-c("orig.ident", "lineage", "cell_type")
## cells are the same; only genes reduced
                      orig.ident   lineage cell_type
AACTGGTCACCAGATT_1     Female Neoblasts  neoblast
ATAAGAGCATAAGACA_1     Female Neoblasts  neoblast
ATTCTACAGTACTTGC_1     Female Neoblasts  neoblast

# create SeuObj
Smorth<-CreateSeuratObject(counts = SmCounts, meta.data=smmeta, min.cells = 2)
saveRDS(Smorth, file="Sm_Orthologs_only.rds")

# 6860 features in 43462 cells



## Create Sj ortholog based object
# read Wendt et al data
Sj<-readRDS("all_integrated_SCT-RPCA_RNA.rds")

# slice count matrix with 1:1 orthologs
SjCounts<-Sj@assays$RNA@counts
allSj<-rownames(Sj) # all genes detected in Sm
ortholist<-read.table("SmSj_1to1.txt", header=F, sep="\t")

# genes in both Sm data and 1:1 ortholog list
uniques<-intersect(allSj, ortholist$V2) 
SjCounts<-SjCounts[uniques, ]
## change gene names
rownames(SjCounts)<-ortholist$V3[match(rownames(SjCounts), ortholist$V2)]

# modify meta data
sjmeta<-Sj@meta.data
sjmeta<-sjmeta[,c(1,18:19)]

## cells are the same; only genes reduced
orig.ident        main_type        cell_type
AAACCCAAGACAACTA-1_1   Female16 neoblast progeny neoblast progeny
AAACCCACAATACCTG-1_1   Female16         tegument         tegument
AAACCCACAGTTCTAG-1_1   Female16           neuron           neuront

# create SeuObj
Sjorth<-CreateSeuratObject(counts = SjCounts, meta.data=sjmeta, min.cells = 2)
saveRDS(Sjorth, file="Sj_Orthologs_only.rds")

#6876 features across 60029 samples within 1 assay



## merged and HarmonyIntegration
merged<-merge(Smorth, y=Sjorth, add.cell.ids=c("Sm", "Sj"))
merged$species<-substr(rownames(merged@meta.data), 1, 2)


suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(RColorBrewer)
  library(gridExtra)
  library(dplyr)
  library(cowplot)
  library(ggpubr)
})

wd<-"/data/work/test/"
seuData<-paste0(wd, "merged_SmSj_ortho.rds")
merged<-readRDS(seuData)

ndims=80

merged <- NormalizeData(merged) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(npcs = ndims, verbose=F)

merged <- IntegrateLayers(object = merged, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE)
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:ndims)
merged <- FindClusters(merged, res=1)
merged <- RunUMAP(merged, dims = 1:ndims, reduction = "harmony")
merged <- JoinLayers(merged) # for DE analysis

saveRDS(merged, file=paste0(seuData, "_seu5log", ndims, "harmony.rds"))
mymeta<-merged@meta.data
mymeta<-cbind(mymeta, merged$umap@cell.embeddings)
write.csv(mymeta, file=paste0(seuData, "_seu5log", ndims, "harmony.meta.csv"))

library(gridExtra)
pdf(paste0(wd, seuData, "_seu5log", ndims, "harmony.Dimplot.pdf"), width=18, height=9)
p1<-DimPlot(merged, group.by="orig.ident", cols=myCol, raster=F)+coord_fixed()
p2<-DimPlot(merged, group.by="species", cols=myCol, raster=F)+coord_fixed()
p3<-DimPlot(merged, group.by="lineage", cols=myCol, raster=F)+coord_fixed()
LabelClusters(p3, id="lineage", color="black", size = 3, repel = T, box.padding = 0.5)
p4<-DimPlot(merged, group.by="main_type", cols=myCol, raster=F)+coord_fixed()
LabelClusters(p4, id="main_type", color="black", size = 3, repel = T, box.padding = 0.5)
p5<-DimPlot(merged, label=T, repel=T, raster=F)+coord_fixed()+NoLegend()
grid.arrange(p1, p2, ncol=2)
grid.arrange(p3, p4, ncol=2)
grid.arrange(p5)
dev.off()

## seurat scores
markers<-read.csv("allcombined_markers_orthologs.csv", header=T, na.string="NA")
allgenes<-rownames(merged)
p1 <- lapply(colnames(markers), function(i) {
	m <- markers[, i]
	m <- intersect(m, allgenes)
	#print(i)
	merged<-AddModuleScore(merged, features=list(m), name=i)
	FeaturePlot(merged, features = paste0(i, "1"), label = T, repel = TRUE, order=T, raster=T) + coord_fixed()+theme_void()+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "RdPu"))+labs(color="SeuratScore") #pt.size=0.05
})

ggsave(
   filename = paste0(seuData, "_seu5log", ndims, "harmony.SeuratScores.pdf"),
   plot = marrangeGrob(p1, nrow=3, ncol=4),
   width = 18, height = 15
)



## add original annotation
# add orginal seurat cluster number
# aggregate by original clusters
smmeta<-read.csv("Wendt2020_SmSC.meta.csv", row.names = "X")
rownames(smmeta)<-paste0("Sm_", rownames(smmeta))
smmeta$name<-rownames(smmeta)
sjmeta<-read.csv("all_integrated_SCT-RPCA.meta.csv", row.names="X")
rownames(sjmeta)<-paste0("Sj_", rownames(sjmeta))
sjmeta$name<-rownames(sjmeta)

combmeta<-rbind(smmeta[,c("name","anno")], sjmeta[,c("name","anno")])

merged$orig.anno<-combmeta$anno[match(rownames(merged@meta.data), rownames(combmeta))]


## Unify main cell types (lower/upper case)
library(ggplot2)
library(plotly)
library(RColorBrewer)

myCol<-unique(c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"), brewer.pal(12, "Set3"), brewer.pal(9, "Set1"), brewer.pal(8, "Accent"), brewer.pal(8, "Set2"), brewer.pal(8, "Pastel2"),brewer.pal(9, "Pastel1")))

########## Wendt et al metadata
smmeta<-read.csv("Wendt2020_SmSC.meta.csv", row.names="X")
ggplot(smmeta,aes(x=lineage, y=1, fill=anno))+geom_col(position="fill")+
  scale_y_continuous(labels = scales::percent)+scale_fill_manual(values = myCol)+
  labs(x="Lineage", y="Proportion", fill="Cell type")+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# split Germline into GSCs, GSC progeny, male gametes and female gametes
smmeta$main_type<-smmeta$lineage
smmeta$main_type[which(smmeta$main_type=="Germline")]<-smmeta$anno[which(smmeta$main_type=="Germline")]

# split Vetellaria into S1, S1 progeny, early vitellocytes, late vitellocytes, mature vitellocytes
smmeta$main_type[which(smmeta$main_type=="Vitellaria")]<-smmeta$anno[which(smmeta$main_type=="Vitellaria")]
table(smmeta$main_type)

############# Harmony integrated metadata
mymeta<-read.csv("merged_SmSj_ortho.rds_seu5log80harmony.meta.csv", row.names="X")
mymeta$main_type<-gsub("GSC", "gsc", mymeta$main_type)
mymeta$main_type<-gsub("GSC progeny", "gsc progeny", mymeta$main_type)
mymeta$main_type<-gsub("S1", "s1", mymeta$main_type)
mymeta$main_type<-gsub("S1 progeny", "s1 progeny", mymeta$main_type)

# add sm main type to merged main type
mymeta$main_type[is.na(mymeta$main_type)]<-smmeta$main_type[match(rownames(mymeta), paste0("Sm_", rownames(smmeta)))]

mymeta$main_type<-gsub("early vitellocytes", "Vitellocytes", mymeta$main_type)
mymeta$main_type<-gsub("late vitellocytes", "Vitellocytes", mymeta$main_type)
mymeta$main_type<-gsub("mature vitellocytes", "Vitellocytes", mymeta$main_type)
mymeta$main_type<-gsub("male gametes", "Male gametes", mymeta$main_type)
mymeta$main_type<-gsub("female gametes", "Female gametes", mymeta$main_type)

mymeta<-cbind(mymeta, merged$umap@cell.embeddings)
table(mymeta$main_type)


p4<-DimPlot(merged, group.by="main_type", cols=myCol, raster=F)+theme_void()+coord_fixed()+NoLegend()
LabelClusters(p4, id="main_type", color="black", size = 3, repel = T, box.padding = 0.5)

p<-ggplot(mymeta, aes(umap_1, umap_2))+geom_point(aes(color=main_type), size=0.2)+scale_color_manual(values = myCol)+theme_void()+coord_fixed()
library(plotly)
ggplotly(p)


