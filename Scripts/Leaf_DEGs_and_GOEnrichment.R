setwd("/home/taslima/data/JuengerLab/Research_Article_Preps/NatVariation_salinity_adaptation_Phallii/")
#Clean WorkPlace
rm(list=ls())

#Load Libraries
library(DESeq2)
library(ggplot2)
library("RColorBrewer")
#library("pheatmap")
#library(ggfortify)
library(devtools)
library("BiocParallel")
library("vsn")
library(gridExtra)
library(tidyverse)
library(xlsx)
source("Scripts/RuntopGo.R")

#Assign 4-cores
register(MulticoreParam(workers=4))

#Load Datafile
Design<-read.csv("Data/PHSALT_Leaf_TAGSeq_Design.csv",row.names = 1)
dat<-read.csv("Data/PHSALT_Leaf_TAGSeq_Count.csv",check.names = F,row.names =1)

Design<-Design[which(rownames(Design) %in% colnames(dat)),]

countTable<-dat[,rownames(Design)]

#Filter genes with a mean count of >=4
countTable<-countTable[-which(as.vector(rowMeans(countTable,na.rm = F)) < 4),]

dds <- DESeqDataSetFromMatrix(countData = countTable,
                              colData = Design,
                              design = ~ Treatment + Genotype +Treatment*Genotype )
dds<-DESeq(dds)

dds = estimateSizeFactors(dds)
sizeFactors(dds)


seC<-SummarizedExperiment(assays = data.matrix(countTable),
                          colData = DataFrame(Design))


flist = list(c("~ Genotype + Treatment","~ Treatment","Genotype_HAL_vs_FIL"),
             c("~ Genotype + Treatment","~ Genotype","Treatment_Salinity_vs_Control"),
             c("~  Genotype + Treatment + Treatment*Genotype","~Genotype + Treatment ","Treatment_Salinity_vs_Control"))

#names= c("Var_hallii_vs_filipes","Treatment_Salt_vs_Control","Treatment_Salt_vs_Control")

geneIDsC = rownames(countTable)

outC= lapply(flist, function(x){
  print(x)
  dds<- DESeqDataSet(se = seC, design = as.formula(x[1]))
  des<-DESeq(dds,test="LRT", reduced= as.formula(x[2]), parallel = T)
  print(results(des,name = x[3]))
  resAll<-data.frame(gene=geneIDsC, results(des,name = x[3]))
  return(resAll)
})

ids = c("_G","_T","_GxT")
for(i in 1:3){
  colnames(outC[[i]])[-c(1:2)]<-paste0(colnames(outC[[i]])[-c(1:2)],ids[i])
}

#length(which(outC[[1]][,7] < (0.05/6)))

DEGresultLeaf<-merge(
  merge(outC[[1]][,c('gene','log2FoldChange_G','pvalue_G','padj_G')],
      outC[[2]][,c('gene','log2FoldChange_T','pvalue_T','padj_T')] ,by="gene",all=T),
  outC[[3]][,c('gene','log2FoldChange_GxT','pvalue_GxT','padj_GxT')] ,by="gene",all=T)

write.xlsx(DEGresultLeaf,file = "Results/Supporting_Leaf_DEG.xlsx",row.names = F)

######## Write DEG
library(xlsx)
info<-read.delim("DBs/PhalliiHAL_496_v2.1.annotation_info.txt",na.strings = "")

dim(info)
info<-info[,c(2,13,16)]
info<-info[!duplicated(info$locusName),]
dim(info)

# G
GGene<-as.data.frame(outC[[1]][which(outC[[1]]$padj_G < (0.05/1)),])

#sigGElist<-outC[[3]][which(outC[[3]][,7] < (0.05/6)),]
GGeneAnnot<-merge(GGene, info, by.x="gene", by.y="locusName",all.y=F)

write.xlsx(GGeneAnnot,file="Results/Salt_Leaf_DEG.xlsx",sheetName = "GGene",append = F)

#E Gene
EGene<-as.data.frame(outC[[2]][which(outC[[2]]$padj_T < (0.05/1)),])

#sigGElist<-outC[[3]][which(outC[[3]][,7] < (0.05/6)),]
EGeneAnnot<-merge(EGene, info, by.x="gene", by.y="locusName",all.y=F)

write.xlsx(EGeneAnnot,file="Results/Salt_Leaf_DEG.xlsx",sheetName = "EGene",append = T)

# GE Genes
GEGene<-as.data.frame(outC[[3]][which(outC[[3]]$padj_GxT < (0.05/1)),])

#sigGElist<-outC[[3]][which(outC[[3]][,7] < (0.05/6)),]
GEGeneAnnot<-merge(GEGene, info, by.x="gene", by.y="locusName",all.y=F)

write.xlsx(GEGeneAnnot,file="Results/Salt_Leaf_DEG.xlsx",sheetName = "GEGene",append = T)


####################################
DEG_Shoot_G<-openxlsx::read.xlsx("Results/Salt_Leaf_DEG.xlsx",sheet = 1)
DEG_Shoot_G<-DEG_Shoot_G[,c('gene','baseMean')]
DEG_Shoot_G$Leaf_G<-rep(1,dim(DEG_Shoot_G)[1])

DEG_Shoot_E<-openxlsx::read.xlsx("Results/Salt_Leaf_DEG.xlsx",sheet = 2)
DEG_Shoot_E<-DEG_Shoot_E[,c('gene','baseMean')]
DEG_Shoot_E$Leaf_E<-rep(1,dim(DEG_Shoot_E)[1])

DEG_Shoot_GE<-openxlsx::read.xlsx("Results/Salt_Leaf_DEG.xlsx",sheet = 3)
DEG_Shoot_GE<-DEG_Shoot_GE[,c('gene','baseMean')]
DEG_Shoot_GE$Leaf_GE<-rep(1,dim(DEG_Shoot_GE)[1])

GXEShoot<-DEG_Shoot_GE$gene
GplusEShoot<-setdiff(intersect(DEG_Shoot_G$gene,DEG_Shoot_E$gene),
                     GXEShoot)
GShoot<-setdiff(DEG_Shoot_G$gene,unique(GplusEShoot,GXEShoot))
EShoot<-setdiff(DEG_Shoot_E$gene,unique(GplusEShoot,GXEShoot))
length(unique(c(GShoot,EShoot,GplusEShoot,GXEShoot)))

datS<-as.data.frame(cbind( Gene=as.vector(as.character(
  unique(c(GShoot,EShoot,GplusEShoot,GXEShoot))
)) ))

datS$Category<-ifelse(datS$Gene %in% GXEShoot,"GXE" ,ifelse(
  datS$Gene %in% GplusEShoot, "G+E", ifelse(
    datS$Gene %in% GShoot, "G","E"
  )
))

GFFAnnot<-read.delim(file="DBs/PhalliiHAL_496_v2.1.Annotation.tab",
                     sep = "\t")

#### FOR GO
genedb<-read.delim("DBs/PhalliiHAL_496_v2.1.Onlygene.mod.gff3",
                   header = F,sep = "\t")
geneID2GO <- "DBs/PH_HAL_v2_GO_Sal_Shoot_detected.tab"

info<-read.delim("DBs/PhalliiHAL_496_v2.1.annotation_info.txt",na.strings = "")
head(info)
#Select only columns for geneID and GO
info<-info[,c(2,10,16)]
#Get only genes that have GO annotation
length(which(!is.na(info$GO)))
info<-info[which(!is.na(info$GO)),]

for (i in c("G","E","G+E","GXE")) {
  if (i=="G") {
    GeneID<-datS$Gene[which(datS$Category==i)]
    SigGO<-RunTopGO(geneID2GO,info,GeneID)
    if (nrow(SigGO)==0) {
      SigGO[1,]<-rep(NA,8)
    }
    write.xlsx(SigGO,file="Results/Salt_Leaf_4Cat_GO_enrichment.xlsx",sheetName = i,append = F)
    
  } else if (i=="E") {
    DEG_Shoot_E<-openxlsx::read.xlsx("Results/Salt_Leaf_DEG.xlsx",sheet = 2)
    DEG_Shoot_E<-DEG_Shoot_E[,c('gene','baseMean','log2FoldChange_T')]
    DEG_Shoot_E<-DEG_Shoot_E[which(DEG_Shoot_E$gene %in% datS$Gene[datS$Category=="E"]),]
    
    #### UP
    DEG_Shoot_E_UP<-DEG_Shoot_E[which(DEG_Shoot_E$log2FoldChange_T > 0),]
    GeneID<-DEG_Shoot_E_UP$gene
    SigGO<-RunTopGO(geneID2GO,info,GeneID)
    if (nrow(SigGO)==0) {
      SigGO[1,]<-rep(NA,8)
    }
    write.xlsx(SigGO,file="Results/Salt_Leaf_4Cat_GO_enrichment.xlsx",
               sheetName = paste(i,"UP",sep = "_"),append = T)
    
    ###DOWN
    DEG_Shoot_E_DOWN<-DEG_Shoot_E[which(DEG_Shoot_E$log2FoldChange_T < 0),]
    GeneID<-DEG_Shoot_E_DOWN$gene
    SigGO<-RunTopGO(geneID2GO,info,GeneID)
    if (nrow(SigGO)==0) {
      SigGO[1,]<-rep(NA,8)
    }
    
    write.xlsx(SigGO,file="Results/Salt_Leaf_4Cat_GO_enrichment.xlsx",
               sheetName = paste(i,"DOWN",sep = "_"),append = T)
    
  } else if (i =="G+E") {
    DEG_Shoot_E<-openxlsx::read.xlsx("Results/Salt_Leaf_DEG.xlsx",sheet = 2)
    DEG_Shoot_E<-DEG_Shoot_E[,c('gene','baseMean','log2FoldChange_T')]
    DEG_Shoot_E<-DEG_Shoot_E[which(DEG_Shoot_E$gene %in% datS$Gene[datS$Category=="G+E"]),]
    
    #### UP
    DEG_Shoot_E_UP<-DEG_Shoot_E[which(DEG_Shoot_E$log2FoldChange_T > 0),]
    GeneID<-DEG_Shoot_E_UP$gene
    SigGO<-RunTopGO(geneID2GO,info,GeneID)
    if (nrow(SigGO)==0) {
      SigGO[1,]<-rep(NA,8)
    }
    write.xlsx(SigGO,file="Results/Salt_Leaf_4Cat_GO_enrichment.xlsx",
               sheetName = paste(i,"UP",sep = "_"),append = T)
    
    ###DOWN
    DEG_Shoot_E_DOWN<-DEG_Shoot_E[which(DEG_Shoot_E$log2FoldChange_T < 0),]
    GeneID<-DEG_Shoot_E_DOWN$gene
    SigGO<-RunTopGO(geneID2GO,info,GeneID)
    if (nrow(SigGO)==0) {
      SigGO[1,]<-rep(NA,8)
    }
    
    write.xlsx(SigGO,file="Results/Salt_Leaf_4Cat_GO_enrichment.xlsx",
               sheetName = paste(i,"DOWN",sep = "_"),append = T)
    
  } else if (i=="GXE") {
    GeneID<-datS$Gene[which(datS$Category==i)]
    SigGO<-RunTopGO(geneID2GO,info,GeneID)
    if (nrow(SigGO)==0) {
      SigGO[1,]<-rep(NA,8)
    }
    write.xlsx(SigGO,file="Results/Salt_Leaf_4Cat_GO_enrichment.xlsx",sheetName = i,append = T)
  }

}


