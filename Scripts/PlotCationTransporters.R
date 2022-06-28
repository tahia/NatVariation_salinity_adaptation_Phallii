
setwd("/home/taslima/data/JuengerLab/Research_Article_Preps/NatVariation_salinity_adaptation_Phallii//")
rm(list=ls())
library(DESeq2)
library(ggplot2)
library(edgeR)
library(openxlsx)
library(xlsx)
library('knitr')
library('limma')
library('reshape2')
library("RColorBrewer")
library(devtools)
library("BiocParallel")
library("vsn")
library('WGCNA')
library(tidyverse)
source("Scripts/RuntopGo.R")


register(MulticoreParam(workers=4))

cpucores <- 4
require(parallel)
options("mc.cores"=cpucores)


dat<-read.csv("Data/PHSALT_ConcatTissues_TAGSeq_Count_V4.csv",row.names = 1)

####### FILTER for GXE genes in either Leaf or Root Tissue
library(xlsx)
LGE<-xlsx::read.xlsx("Results/Salt_Leaf_DEG.xlsx",sheetIndex = 3)
RGE<-xlsx::read.xlsx("Results/Salt_Root_DEG.xlsx",sheetIndex = 3)

#ConcatGE<-unique(c(LGE$gene,RGE$gene))
#dat<-dat[which(rownames(dat) %in% ConcatGE),]

Design<-read.csv("Data/Concat_Shoot_Root_Design_V3.csv",row.names = 3)
countTable<-dat[,rownames(Design)]
countTable<-countTable[-which(as.vector(rowMeans(countTable,na.rm = F)) <4),]
#countTable<-countTable[-which(as.vector(rowMeans(countTable,na.rm = F)) <1),] #mean count >5
## None has mean less than 1

dds <- DESeqDataSetFromMatrix(countData = countTable,
                              colData = Design,
                              design = as.formula("~ Treatment + Ecotype + Tissue + 
                                                   Treatment:Ecotype+ Treatment:Tissue+ Ecotype:Treatment +
                                                   Treatment:Ecotype:Tissue"))
dds

dds = estimateSizeFactors(dds)
sizeFactors(dds)

# Variance Stabilizing Transformation
vst<-vst(dds,blind = F)
dat<-as.data.frame(assay(vst))
tdat<-as.data.frame(t(dat))
tdat$LIBID<-rownames(tdat)
Design$LIBID<-rownames(Design)
VSTDes<-merge(Design,tdat,by="LIBID")

CTs<-read.csv("Results/Leaf_Root_GE_CationTransporter.csv")
pdf("Plots/Supporting_Information_Figure_S4.pdf",onefile = TRUE)

for (i in 1:nrow(CTs)) {
  gene<-CTs$GeneID[i]
  Descript<-paste(gene,CTs$Description[i],sep = "\n")
  p<-as_tibble(VSTDes) %>% dplyr::select(Ecotype, Tissue,Treatment, gene) %>% group_by(Ecotype,Treatment,Tissue) %>%
    summarise_all(funs(mean,sd,se=sd(.)/sqrt(n()) )) %>%
    mutate(Group=paste(Tissue,Ecotype,sep = ":"))%>% ggplot()+
    geom_point(aes(x=Treatment, y=mean,color=Ecotype,shape=Tissue),alpha=1,size=3,show.legend = T)+
    geom_line(aes(x=Treatment, y=mean,group=Group,color=Ecotype,linetype=Tissue),alpha=0.75,size=1,show.legend = T)+
    geom_linerange(aes(x=Treatment,ymin=mean-se,ymax=mean+se,color=Ecotype,linetype=Tissue),alpha=0.75,size=0.5,show.legend = F)+
    scale_color_manual(values = c( "#3949AB","#CB4335" ))+
    scale_x_discrete(expand = c(0.1,0.1))+
    labs(x="Treatment",y= "Normalized Count",title = Descript,color="Ecotype")+
    theme_bw(base_size = 14)+
    theme(plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
          legend.title = element_text(size=12, face="bold"),
          legend.text = element_text(size=12),
          axis.title = element_text(size=12,face="bold"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text=element_text(size=12,vjust=1),
          strip.background = element_rect(fill = "transparent", color = NA),
          strip.text = element_text(size=12,face = "bold")
          #axis.ticks.x = element_blank()
    )
  print(p)

}
dev.off()

CTExp<-as_tibble(VSTDes) %>% dplyr::select(LIBID,CTs$GeneID[which(CTs$Selection=="Yes")])
CTExpDes<-merge(Design[,-c(1,3,4,6,7,9:13,17:19)], CTExp,by="LIBID")

CTExpDesLeaf<-CTExpDes %>% filter(Tissue=="Leaf")
MLeaf<-cor(CTExpDesLeaf[,-c(1:4)],use = "complete.obs")
corrplot(MLeaf)

CTExpDesRoot<-CTExpDes %>% filter(Tissue=="Root")
MRoot<-cor(CTExpDesRoot[,-c(1:4)],use = "complete.obs")
corrplot(MRoot)

TargetGene<-as.data.frame(cbind( 
  gene=c("PhHAL.3G288900","PhHAL.2G474100",
         "PhHAL.3G490000","PhHAL.4G025100",
         "PhHAL.2G242200","PhHAL.1G101800",
         "PhHAL.2G471600","PhHAL.9G082500",
         "PhHAL.2G001500"),
  Descript=c("Sodium hydrogen exchanger 2 (NHX2)",
             "Sodium hydrogen exchanger 2 (NHX2)",
             "Sodium proton exchanger (SOS1)",
             "High-affinity K+ transporter 1 (HKT1)",
             "Potassium transporter (KT)",
             "Potassium channel 1 (KAT1)",
             "Potassium transporter (KUP)",
             "Outward rectifying potassium\nchannel (KCO)",
             "High affinity K+ transporter (HAK)")
                                 ))
geneList<-TargetGene$gene
DesList<-TargetGene$Descript

for (i in 1:length(geneList)) {
  gene<-geneList[i]
  Title<-paste0(gene,"\n\n",DesList[i])
  p<-as_tibble(VSTDes) %>% dplyr::select(Ecotype, Tissue,Treatment, gene) %>% 
    group_by(Ecotype,Treatment,Tissue) %>% 
    summarise_all(funs(mean,sd,se=sd(.)/sqrt(n()) )) %>% 
    mutate(Group=paste(Tissue,Ecotype,sep = ":"))%>% 
    mutate(Genotype= ifelse(Ecotype=="FIL","Coastal","Inland")) %>% 
    ggplot()+
    geom_point(aes(x=Treatment, y=mean,color=Genotype,shape=Tissue),alpha=1,size=3,show.legend = F)+
    geom_line(aes(x=Treatment, y=mean,group=Group,color=Genotype,linetype=Tissue),alpha=0.75,size=1,show.legend = F)+
    geom_linerange(aes(x=Treatment,ymin=mean-se,ymax=mean+se,color=Genotype,linetype=Tissue),alpha=0.75,size=0.5,show.legend = F)+
    scale_color_manual(values = c( "#3949AB","#CB4335" ))+
    scale_x_discrete(expand = c(0.1,0.1))+
    labs(x="Treatment",y= "Normalized Count",title = Title,color="Genotype")+
    theme_bw(base_size = 14)+
    theme(plot.title = element_text(hjust = 0.5,size=9),
          legend.title = element_text(size=10, face="bold"),
          legend.text = element_text(size=10),
          axis.title = element_text(size=10,face="bold"),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text=element_text(size=12,vjust=1),
          strip.background = element_rect(fill = "transparent", color = NA),
          strip.text = element_text(size=12,face = "bold")
          #axis.ticks.x = element_blank()
    )
  #print(p)
  assign(paste("plot",gene,sep = "_"), p)
  
}

###### Legend
gene<-geneList[i]
Title<-paste0( gene,"\n\n",DesList[i])

plotLeg1<-as_tibble(VSTDes) %>% dplyr::select(Ecotype, Tissue,Treatment, gene) %>% 
  group_by(Ecotype,Treatment,Tissue) %>% 
  summarise_all(funs(mean,sd,se=sd(.)/sqrt(n()) )) %>% 
  mutate(Group=paste(Tissue,Ecotype,sep = ":"))%>% 
  mutate(Genotype= ifelse(Ecotype=="FIL","Coastal","Inland")) %>% 
  ggplot()+
  geom_point(aes(x=Treatment, y=mean,color=Genotype,shape=Tissue),alpha=1,size=3,show.legend = T)+
  geom_line(aes(x=Treatment, y=mean,group=Group,color=Genotype,linetype=Tissue),alpha=0.75,size=1,show.legend = T)+
  geom_linerange(aes(x=Treatment,ymin=mean-se,ymax=mean+se,color=Genotype,linetype=Tissue),alpha=0.75,size=0.5,show.legend = F)+
  scale_color_manual(values = c( "#3949AB","#CB4335" ))+
  scale_x_discrete(expand = c(0.1,0.1))+
  labs(x="Treatment",y= "Normalized Count",title = Title,color="Geotype")+
  guides(color="none")+
  theme_bw(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        axis.title = element_text(size=12,face="bold"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12,vjust=1),
        strip.background = element_rect(fill = "transparent", color = NA),
        strip.text = element_text(size=12,face = "bold")
        #axis.ticks.x = element_blank()
  )

plotLeg2<-as_tibble(VSTDes) %>% dplyr::select(Ecotype, Tissue,Treatment, gene) %>% 
  group_by(Ecotype,Treatment,Tissue) %>% 
  summarise_all(funs(mean,sd,se=sd(.)/sqrt(n()) )) %>% 
  mutate(Group=paste(Tissue,Ecotype,sep = ":"))%>% 
  mutate(Genotype= ifelse(Ecotype=="FIL","Coastal","Inland")) %>% 
  ggplot()+
  geom_point(aes(x=Treatment, y=mean,color=Genotype,shape=Tissue),alpha=1,size=3,show.legend = T)+
  geom_line(aes(x=Treatment, y=mean,group=Group,color=Genotype,linetype=Tissue),alpha=0.75,size=1,show.legend = T)+
  geom_linerange(aes(x=Treatment,ymin=mean-se,ymax=mean+se,color=Genotype,linetype=Tissue),alpha=0.75,size=0.5,show.legend = F)+
  scale_color_manual(values = c( "#3949AB","#CB4335" ))+
  scale_x_discrete(expand = c(0.1,0.1))+
  labs(x="Treatment",y= "Normalized Count",title = Title,color="Geotype")+
  guides(shape="none",linetype="none")+
  theme_bw(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        axis.title = element_text(size=12,face="bold"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12,vjust=1),
        strip.background = element_rect(fill = "transparent", color = NA),
        strip.text = element_text(size=12,face = "bold")
        #axis.ticks.x = element_blank()
  )


legendT <- get_legend(plotLeg1)
legendG <- get_legend(plotLeg2)


ConcatIon<-ggdraw() +
  draw_plot(plot_PhHAL.3G288900, x = 0, y = 0.66, width = 0.33,height = 0.33)+
  draw_label(x=0.03,y=0.97,label = "A)", color = "black", size = 12, fontface = "bold")+
  draw_plot(plot_PhHAL.2G474100, x = 0.33, y = 0.66, width = 0.33,height = 0.33)+
  draw_label(x=0.36,y=0.97,label = "B)", color = "black", size = 12, fontface = "bold")+
  draw_plot(plot_PhHAL.3G490000, x = 0.66, y = 0.66, width = 0.33,height = 0.33)+
  draw_label(x=0.69,y=0.97,label = "C)", color = "black", size = 12, fontface = "bold")+
  draw_plot(plot_PhHAL.4G025100, x = 0, y = 0.33, width = 0.33,height = 0.33)+
  draw_label(x=0.03,y=0.65,label = "D)", color = "black", size = 12, fontface = "bold")+
  draw_plot(plot_PhHAL.2G242200, x = 0.33, y = 0.33, width = 0.33,height = 0.33)+
  draw_label(x=0.36,y=0.65,label = "E)", color = "black", size = 12, fontface = "bold")+
  draw_plot(plot_PhHAL.1G101800,x = 0.66, y = 0.33,width = 0.33,height = 0.33) + 
  draw_label(x=0.69,y=0.65,label = "F)", color = "black", size = 12, fontface = "bold")+
  draw_plot(plot_PhHAL.2G471600,x = 0, y = 0,width = 0.33,height = 0.33)+
  draw_label(x=0.03,y=0.32,label = "G)", color = "black", size = 12, fontface = "bold")+
  draw_plot(plot_PhHAL.9G082500,x = 0.33, y = 0,width = 0.33,height = 0.33)+
  draw_label(x=0.36,y=0.32,label = "H)", color = "black", size = 12, fontface = "bold")+
  draw_plot(plot_PhHAL.2G001500,x = 0.66, y = 0,width = 0.33,height = 0.33)+
  draw_label(x=0.69,y=0.32,label = "I)", color = "black", size = 12, fontface = "bold")+
  draw_plot(legendG, x = 0.21, y = -0.05, width = 0.5, height = 0.5)+
  draw_plot(legendT, x = 0.03, y = -0.14, width = 0.5, height = 0.5)


tiff("Plots/Figure_IonTransporters.tiff",width=9,height=9,units="in",res=100)
ConcatIon
dev.off()
