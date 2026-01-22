library(Seurat)
library(spacexr)
library(SingleR)
library(SCINA)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

myCol<-unique(c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"), brewer.pal(12, "Set3"),
                brewer.pal(8, "Set1"), brewer.pal(8, "Accent")))

scData<-"all_integrated_SCT-RPCA_RNA.rds"
spData<-"G214.adjusted.cellbin.rds"

sc.data<-readRDS(scData)
sp.data<-readRDS(spData)

message("RCTD:")
sc.counts<-sc.data@assays$RNA@counts
meta_data<-sc.data@meta.data

# set to default idents
cell_types<-as.factor(sc.data$cell_type); 
names(cell_types)<-rownames(meta_data)
nUMI<-meta_data$nCount_RNA; names(nUMI)<-rownames(meta_data)

reference<-Reference(sc.counts, cell_types, nUMI)

sp.counts<-sp.data@assays$RNA@layers$counts
rownames(sp.counts)<-rownames(sp.data)
colnames(sp.counts)<-colnames(sp.data)
coords<-sp.data@meta.data[, c("coord_x", "coord_y")]
sp.nUMI<-colSums(sp.counts)
puck<-SpatialRNA(coords, sp.counts, sp.nUMI)

barcodes <- colnames(puck@counts)
plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), title ='plot of nUMI') 

myRCTD<-create.RCTD(puck, reference, max_cores = 1) # for parallel processing, the number of cores used. 1 means no parallel processing.
myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet") 
save(myRCTD, file=paste0(spData, "_myRCTD.rda"))


results<-myRCTD@results
# normalize the cell type proportions to sum to 1.
norm_weights = normalize_weights(results$weights) 
cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- myRCTD@spatialRNA

## combine max weight with meta data
norm_weights = normalize_weights(results$weights); norm_weights<-as.data.frame(norm_weights)
#noReject<-rownames(subset(results$results_df, subset=spot_class!="reject"))
#norm_weights<-norm_weights[which(rownames(norm_weights) %in% noReject),]

maxcol<-colnames(norm_weights)[max.col(norm_weights, ties.method = "first")]
maxweight<-cbind(rownames(norm_weights), maxcol); 
maxweight<-as.data.frame(maxweight)

table(maxweight$maxcol)

sp.data@meta.data$RCTD<-maxweight$maxcol[match(rownames(sp.data@meta.data), maxweight$V1)]
sp.data$RCTD[which(is.na(sp.data$RCTD))]<-"ambiguous"


sp.data$cell_type<-gsub("ambiguous1", "esophageal gland", sp.data$RCTD)
sp.data$cell_type<-gsub("ambiguous2", "ambiguous", sp.data$cell_type)
Idents(sp.data)<-"cell_type"

sp.data$orig.ident<-"SjD26"


saveRDS(sp.data, file=paste0(spData, "_RCTD-cell_type.rds"))


