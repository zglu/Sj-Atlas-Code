## cross-species seurat
merged<-readRDS("merged_SmSj_ortho.rds_seu5log80harmony.rds")
merged

pdf("merged_pca.pdf", width=12, height=9) # doesn't tell much
DimPlot(merged, reduction = "pca", group.by = "main_type", label = T)+NoLegend()
dev.off()

## Pearson correlation heatmap
merged$celltype_species <- paste0(merged$main_type, "_", merged$species)
Idents(merged) <- "celltype_species"

# Calculate average expression for a set of variable genes (or all genes)
avg_expression <- log1p(AverageExpression(merged, assays = "RNA", slot = "data")$RNA)

# Calculate the Pearson correlation matrix
cor_matrix <- cor(as.matrix(avg_expression))

# Convert correlation to a distance (1 - correlation) for plotting
dist_matrix <- 1 - cor_matrix

library(pheatmap)
pdf("merged_pearson_mainType.pdf", width=13, height=9)
p<-pheatmap(
  cor_matrix,
  clustering_distance_rows = as.dist(dist_matrix),
  clustering_distance_cols = as.dist(dist_matrix),
  display_numbers = TRUE, # Set to TRUE if you want to see the correlation values
  main = "Cell type correlation across species"
)
print(p)
dev.off()

########################
#### selected cell types
## select correlation heatmap
# take only corresponding main types
subDist<-cor_matrix[grepl("Sj", rownames(cor_matrix)), grepl("Sm", colnames(cor_matrix))]
subDist<-subDist[!grepl("ambiguous", rownames(subDist)),!grepl("Unknown", colnames(subDist))]
subDist<-subDist[,!grepl("Oesopha", colnames(subDist))]
# the neoblast and progeny order differrnt in Sm and Sj
subCor<- subDist[ order(row.names(subDist)), ]
subCor<- subCor[ , order(colnames(subCor))]

p2<-pheatmap(subCor, cluster_rows = F, cluster_cols = F, display_numbers = T, legend=F)

pdf("merged_pearson_mainType_selected.pdf", width=6, height=6)
print(p2)
dev.off()

