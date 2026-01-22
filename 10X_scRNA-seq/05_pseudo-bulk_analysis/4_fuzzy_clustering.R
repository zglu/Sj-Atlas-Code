dyna<-read.csv("pseudo_main_type_Sample_allDEGs_comb_dynamic.csv", row.names = "gene")

# pseudo count
pseudo<-readRDS("all_integrated_dpi-sex_pseudo.rds")

#Normalization steps : From the input count tables, the Mfuzz_RNAseq.R script performs a library size normalization with DESeq method and then adjust these normalized data for gene length (normalized data / gene length). These normalization steps are carried out to make all the samples comparable, which is required by Mfuzz package.
# in single cell can we compare genes?
psdata<-pseudo@assays$RNA@layers$data
colnames(psdata)<-colnames(pseudo)
rownames(psdata)<-rownames(pseudo)
rownames(psdata)<-gsub("EWB00-", "EWB00_", rownames(psdata))
psdata<-as.data.frame(psdata)
# select only dynamic genes
dydata<-subset(psdata, subset=rownames(psdata)%in%rownames(dyna))

# select only male/female samples
subdata<-dydata[,c(1,3,5,2,4,6)]
# change to cell x gene data frame
expre<-as.matrix(subdata)
# build mfuzz object
exprSet=ExpressionSet(assayData=expre)
# filter NAs
exprSet.r <- filter.NA(exprSet,thres = 0.25)
# fill NAs
exprSet.f <- fill.NA(exprSet.r,mode = 'mean')

# mfuzz
exprSet.s=standardise(exprSet.f)
m1=mestimate(exprSet.s)
cl=mfuzz(exprSet.s,c=8,m=m1) # define the number of clusters

pdf("mfuzz-female-male.pdf")
mfuzz.plot2(exprSet.s, cl=cl,mfrow=c(4,4),centre=TRUE,x11=F,centre.lwd=0.2)
dev.off()

#export genes in each cluster
dir.create(path="mfuzz",recursive = TRUE)
for(i in 1:8){
  potname<-names(cl$cluster[unname(cl$cluster)==i])
  write.csv(cl[[4]][potname,i],paste0("mfuzz","/mfuzz_",i,".csv"))
}
