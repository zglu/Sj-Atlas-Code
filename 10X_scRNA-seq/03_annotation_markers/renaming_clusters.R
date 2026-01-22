wd<-"./"
merged<-readRDS("all_integrated_SCT-RPCA_noSCT.rds")
Idents(merged)<-"seurat_clusters"
cluster.ids<-c(
  "C0_neoblast", 
  "C1_muscle",
  "C2_tegument progenitor",
  "C3_tegument",
  "C4_neoblast progeny",
  "C5_muscle",
  "C6_parenchyma",
  "C7_neoblast progeny",
  "C8_muscle",
  "C9_neuron",
  "C10_neoblast", 
  "C11_muscle",
  "C12_neuron",
  "C13_neoblast",
  "C14_neuron",
  "C15_neoblast progeny",
  "C16_neoblast",
  "C17_neoblast progeny",
  "C18_parenchyma",
  "C19_vitellocyte.early",
  "C20_tegument",
  "C21_muscle",
  "C22_neoblast progeny",
  "C23_S1",
  "C24_ambiguous1", 
  "C25_tegument",
  "C26_S1 progeny",
  "C27_gut",
  "C28_neoblast progeny",
  "C29_vitellocyte.mature",
  "C30_neoblast progeny",
  "C31_tegument",
  "C32_muscle",
  "C33_tegument progenitor",
  "C34_neoblast",
  "C35_muscle",
  "C36_neuron",
  "C37_parenchyma",
  "C38_tegument progenitor",
  "C39_neuron",
  "C40_muscle",
  "C41_GSC progeny",
  "C42_GSC",
  "C43_flame",
  "C44_tegument",
  "C45_neuron",
  "C46_flame",
  "C47_female late",
  "C48_neoblast progeny",
  "C49_neuron",
  "C50_vitellocyte.late",
  "C51_male late",
  "C52_tegument",
  "C53_mehlis gland",
  "C54_tegument progenitor",
  "C55_ambiguous2", 
  "C56_tegument",
  "C57_neoblast progeny",
  "C58_neuron",
  "C59_neuron",
  "C60_flame",
  "C61_neuron"

)


names(x=cluster.ids)<-levels(x=merged)
merged <- RenameIdents(object = merged, cluster.ids)
merged$anno<-Idents(merged)

# main cell type
merged$main_type<-gsub(".*[_]([^.]+).*", "\\1", merged$anno)
# cell type: after _
merged$cell_type<-gsub(".*[_](.+).*", "\\1",merged$anno)


# add umap embedding to metadata
mymeta<-merged@meta.data
mymeta<-cbind(mymeta, merged$umap@cell.embeddings)
write.csv(mymeta, file=paste0(wd, "all_integrated_SCT-RPCA.meta.csv"))

# remove integrate assay
DefaultAssay(merged)<-"RNA"
merged[['integrated']]<-NULL

saveRDS(merged, file=paste0(wd, "all_integrated_SCT-RPCA_RNA.rds"))
