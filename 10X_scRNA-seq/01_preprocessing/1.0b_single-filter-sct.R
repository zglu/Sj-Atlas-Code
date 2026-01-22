# Rscript 1.0b_single-filter-sct.R [F16/F26/M16/M26]

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(gridExtra)
  library(ggplot2)
  library(ggsci)
  library(ggpubr)
  library(RColorBrewer)
})

args<-commandArgs(T)


myCol<-unique(c(pal_d3("category10")(10),pal_rickandmorty("schwifty")(12), pal_lancet("lanonc")(9),
                pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),
                pal_jama("default")(7),pal_jco("default")(10),
                pal_locuszoom("default")(7),pal_startrek("uniform")(7),
                pal_tron("legacy")(7),pal_futurama("planetexpress")(12),
                pal_simpsons("springfield")(16),
                pal_gsea("default")(12)))



# raw data already pre-filtered with min.cells=5, min.features=500, and percent.mt<5
# read from rds with doublet detection
rawData<-readRDS(paste0("rawRDS/", args[1], "_mt5_gene500_raw.rds"))

rawData[["ident"]]<-NULL

rawData

head(rawData@meta.data)

# Distribution plot
# novelty score
rawData$log10GenesPerUMI <- log10(rawData$nFeature_RNA) / log10(rawData$nCount_RNA)

p3<-ggplot(rawData@meta.data,aes(x=nCount_RNA, y=nFeature_RNA, color=log10GenesPerUMI)) +
  geom_point() +scale_colour_gradient(low = "gray90", high = "skyblue") +theme_classic()+labs(title = args[1])+stat_cor(method="pearson")+
  geom_vline(xintercept=25000, linetype="dashed", color = "lightblue", linewidth=0.5)+
  geom_hline(yintercept=3000, linetype="dashed", color = "lightblue", linewidth=0.5)

p4<-ggplot(rawData@meta.data,aes(x=nFeature_RNA, y=percent.mt, color=nCount_RNA)) +
  geom_point() +scale_colour_gradient(low = "gray90", high = "red") +theme_classic()+labs(title = args[1])+stat_cor(method="pearson")+
  geom_vline(xintercept=1000, linetype="dashed", color = "red", linewidth=0.5)+
  geom_hline(yintercept=5, linetype="dashed", color = "red", linewidth=0.5)

p5<-ggplot(rawData@meta.data, aes(x=log10GenesPerUMI)) +
  geom_density(alpha = 0.4, fill="khaki1") + labs(title = "log10 Genes per UMI")+
  theme_classic()

p6<-ggplot(rawData@meta.data, aes(x=percent.mt)) +
  geom_histogram(alpha=0.4, color="black",fill="hotpink")+theme_classic()+xlim(0, 10)+labs(title="mito ratio")
#  geom_density(alpha = 0.4, fill="hotpink") + labs(title = "mito ratio")+theme_classic()

p7<-ggplot(rawData@meta.data, aes(x=nFeature_RNA)) +
  geom_density(alpha = 0.4, fill="lightblue")+labs(title = "distribution of nGene")+
  theme_classic()

p8<-ggplot(rawData@meta.data, aes(x=nCount_RNA)) +
  geom_density(alpha = 0.4, fill="bisque")+labs(title = "distribution of nUMI")+
  theme_classic()


### the mean-variance relationship for each gene
rawcounts<-rawData@assays$RNA@counts

gene_attr <- data.frame(mean = rowMeans(rawcounts), detection_rate = rowMeans(rawcounts>0), var = apply(rawcounts, 1, var))
gene_attr$log_mean <- log10(gene_attr$mean) 
gene_attr$log_var <- log10(gene_attr$var)
rownames(gene_attr) <- rownames(rawcounts)

cell_attr <- data.frame(n_umi = colSums(rawcounts), n_gene = colSums(rawcounts>0))
rownames(cell_attr) <- colnames(rawcounts)

p_mean_var.raw<-ggplot(gene_attr, aes(log_mean, log_var)) + geom_point(alpha = 0.3, shape = 16) + 
            geom_density_2d(size = 0.3) + geom_abline(intercept = 0, slope = 1, color = "red")+labs(title="raw counts")+theme_classic()


p.dbl<-ggplot(rawData@meta.data, aes(x=nCount_RNA, y=nFeature_RNA, color=scDblFinder.class))+
  geom_point(size=0.5)+theme_classic()+labs(title=args[1])

pdf(paste0(args[1], "_rawQC.pdf"), width=18, height=8)
grid.arrange(p3, p4, p7, p8, p5, p6, p_mean_var.raw, p.dbl, nrow=2, ncol=4, top="raw data (min.features=500; mt<5)")
dev.off()

message("QC plot of raw data exported.")

## filtering: take only singlets
rawData<-rawData[, rawData@meta.data$scDblFinder.class=="singlet"]

# set thresholds for individual metrics to be as permissive as possible, and always consider the joint effects
minFeature<-500
maxFeature<-5000
minCount<-1000
maxCount<-30000 #F26 50000
mtRatio<-3
# 

message("Filtering conditions:", minFeature, " < nFeature < ", maxFeature, " & ", minCount, " < nCount < ", maxCount, " & mito ratio < ", mtRatio)

Fil_sct <- subset(x = rawData, subset = nFeature_RNA > minFeature & nFeature_RNA < maxFeature & nCount_RNA > minCount & nCount_RNA < maxCount & percent.mt < mtRatio) 
Fil_sct

# Distribution plot
Fil_sct$log10GenesPerUMI <- log10(Fil_sct$nFeature_RNA) / log10(Fil_sct$nCount_RNA)

p3<-ggplot(Fil_sct@meta.data,aes(x=nCount_RNA, y=nFeature_RNA, color=log10GenesPerUMI)) +
  geom_point() +scale_colour_gradient(low = "gray90", high = "skyblue") +theme_classic()+labs(title = paste0(args[1], "_filtered"))+stat_cor(method="pearson")
p4<-ggplot(Fil_sct@meta.data,aes(x=nFeature_RNA, y=percent.mt, color=nCount_RNA)) +
  geom_point() +scale_colour_gradient(low = "gray90", high = "red") +theme_classic()+labs(title = paste0(args[1], "_filtered"))+stat_cor(method="pearson")

p5<-ggplot(Fil_sct@meta.data, aes(x=log10GenesPerUMI)) +
  geom_density(alpha = 0.4, fill="khaki1") + labs(title = "log10 Genes per UMI")+
  theme_classic()


p7<-ggplot(Fil_sct@meta.data, aes(x=nFeature_RNA)) +
  geom_density(alpha = 0.4, fill="lightblue")+labs(title = "distribution of nGene")+
  theme_classic()

p8<-ggplot(Fil_sct@meta.data, aes(x=nCount_RNA)) +
  geom_density(alpha = 0.4, fill="bisque")+labs(title = "distribution of nUMI")+
  theme_classic()


## normalisation
message("Normalisation:")

Fil_sct %>%  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() -> Fil_sct 

Fil_sct <- SCTransform(Fil_sct, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = F) %>%
    RunPCA() %>%
    RunUMAP(dims = 1:50) %>%
    FindNeighbors(dims = 1:50) %>%
    FindClusters(res=2)

table(Fil_sct@active.ident)

# feature plot as separated by clusters
p1<-FeatureScatter(Fil_sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols=myCol);
p2<-FeatureScatter(Fil_sct, feature1 = "nFeature_RNA", feature2 = "percent.mt", cols=myCol);

#mean-var for SCT
p_mean_var.sct<-ggplot(Fil_sct@assays$SCT@SCTModel.list$model1@feature.attributes, aes(log10(gmean), residual_variance))+
      geom_point(alpha=0.5)+theme_classic()+labs(title="after SCT")+geom_hline(yintercept = 1, color="red")

pelbow<-ElbowPlot(Fil_sct, ndims = 50, reduction = "pca")
pdim<-DimPlot(Fil_sct, reduction = "umap", label = T, repel = T, cols=myCol)+coord_fixed()+NoLegend()+labs(title = args[1])

pdf(paste0(args[1], "_filQC.pdf"), width=18, height=8)
grid.arrange(p3, p4, p7, p8, pelbow, p_mean_var.sct, nrow=2, ncol=4, top=paste0("filtered data: ", minFeature, " < nFeature < ", maxFeature, " & ", minCount, " < nCount < ", maxCount, " & mito ratio < ", mtRatio))
grid.arrange(p1, p2, pdim, ncol=3)
dev.off()


saveRDS(Fil_sct, file=paste0(args[1],"_Fil_SCT_PC50.rds"))


#### Seurat addModuleScore ####
DefaultAssay(Fil_sct) <- "RNA"
message("Calculating scores using addModuleScore:")

allgenes<-rownames(Fil_sct)

markers<-read.csv("allcombined_markers_1e-20_AUC0.7_grouped.csv", header=T)

p1 <- lapply(colnames(markers), function(i) {
  m <- markers[, i]
  m <- intersect(m, allgenes)
  Fil_sct<-AddModuleScore(Fil_sct, features=list(m), name=i)
  FeaturePlot(Fil_sct, features = paste0(i, "1"), label = TRUE, repel = TRUE, raster=TRUE) + coord_fixed()+theme_void()+scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlOrRd"))+labs(color="SeuratScore")
})

ggsave(
  filename = paste0(args[1], "_SeuratScores.pdf"),
  plot = marrangeGrob(p1, nrow=2, ncol=4),
  width = 18, height = 8
)

library(qpdf)
qpdf::pdf_combine(input = c(paste0(args[1], "_rawQC.pdf"), paste0(args[1], "_filQC.pdf"), paste0(args[1], "_SeuratScores.pdf")), output = paste0(args[1], "_QC-sigScores.pdf"))

system("rm -f *QC.pdf *SeuratScores.pdf")
