setwd("/home/taslima/data/JuengerLab/Research_Article_Preps/NatVariation_salinity_adaptation_Phallii/")
library(DESeq2)
library(ggplot2)
library("RColorBrewer")
library(devtools)
library("BiocParallel")
library("vsn")
library(gridExtra)
library(tidyverse)
library(cowplot)
library(FactoMineR)
library(eulerr)
library(adegenet)

#source("Scripts/RunPrettyPCA.R")

#Assign 4-cores
register(MulticoreParam(workers=4))


######################### DAPC

ThemeDAPC<-theme(plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
                 legend.title = element_text(size=10, face="bold"),
                 legend.text = element_text(size=10),
                 legend.position=c(0.5, 0.65),
                 axis.title = element_text(size=10,face="bold",colour = "black"),
                 panel.border = element_blank(), 
                 #panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), 
                 axis.line = element_line(colour = "black"),
                 axis.text=element_text(size=10,vjust=1,face = "bold",colour = "black"),
                 plot.margin = margin(0.15, 0.15, 0.2, 0.15, "cm")
                 # Change margin here to add/remove extra margin
                 #axis.ticks.x = element_blank()
)


######## Leaf

#Load Datafile
DesignL<-read.csv("Data/PHSALT_Leaf_TAGSeq_Design.csv",row.names = 1) 
  
datL<-read.csv("Data/PHSALT_Leaf_TAGSeq_Count.csv",check.names = F,row.names =1)

DesignL<-DesignL[which(rownames(DesignL) %in% colnames(datL)),]

countTableL<-datL[,rownames(DesignL)]

#Filter genes with a mean count of >=1
countTableL<-countTableL[-which(as.vector(rowMeans(countTableL,na.rm = F)) <4),]

ddsL <- DESeqDataSetFromMatrix(countData = countTableL,
                               colData = DesignL,
                               design = ~ Treatment + Genotype +Treatment*Genotype )
ddsL

ddsL = estimateSizeFactors(ddsL)
sizeFactors(ddsL)

vstL<-vst(ddsL,blind = T)
adatL<-as.data.frame(assay(vstL))
tdatL<-as.data.frame(t(adatL))
tdatL$LIBID<-rownames(tdatL)
DesignL$LIBID<-rownames(DesignL)
VSTDesL<-merge(DesignL,tdatL,by="LIBID")


######## Root

#Load Datafile
DesignR<-read.csv("Data/PHSALT_Root_TAGSeq_Design.csv",row.names = 1)
datR<-read.csv("Data/PHSALT_Root_TAGSeq_Count.csv",check.names = F,row.names =1)
DesignR<-DesignR[which(rownames(DesignR) %in% colnames(datR)),]

countTableR<-datR[,rownames(DesignR)]

#Filter genes with a mean count of >=4
countTableR<-countTableR[-which(as.vector(rowMeans(countTableR,na.rm = F)) <4),]

ddsR <- DESeqDataSetFromMatrix(countData = countTableR,
                               colData = DesignR,
                               design = ~ Treatment + Genotype +Treatment*Genotype )
ddsR

ddsR = estimateSizeFactors(ddsR)
sizeFactors(ddsR)


vstR<-vst(ddsR,blind = T)
adatR<-as.data.frame(assay(vstR))
tdatR<-as.data.frame(t(adatR))
tdatR$LIBID<-rownames(tdatR)
DesignR$LIBID<-rownames(DesignR)
VSTDesR<-merge(DesignR,tdatR,by="LIBID")

######################
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
GShoot<-setdiff(DEG_Shoot_G$gene,unique(c(GplusEShoot,GXEShoot)))
EShoot<-setdiff(DEG_Shoot_E$gene,unique(c(GplusEShoot,GXEShoot)))
length(unique(c(GShoot,EShoot,GplusEShoot,GXEShoot)))



datS<-as.data.frame(cbind( Gene=as.vector(as.character(
  unique(c(GShoot,EShoot,GplusEShoot,GXEShoot))
)) ))

datS$Category<-ifelse(datS$Gene %in% GXEShoot,"GXT" ,ifelse(
  datS$Gene %in% GplusEShoot, "G+T", ifelse(
    datS$Gene %in% GShoot, "G","T"
  )
))



VennListLeaf<-list(G=DEG_Shoot_G$gene,E=DEG_Shoot_E$gene,GXT=DEG_Shoot_GE$gene)
VennLeaf<-plot(euler(VennListLeaf),quantities=T)

### Root

DEG_Root_G<-openxlsx::read.xlsx("Results/Salt_Root_DEG.xlsx",sheet = 1)
DEG_Root_G<-DEG_Root_G[,c('gene','baseMean')]
DEG_Root_G$Root_G<-rep(1,dim(DEG_Root_G)[1])

DEG_Root_E<-openxlsx::read.xlsx("Results/Salt_Root_DEG.xlsx",sheet = 2)
DEG_Root_E<-DEG_Root_E[,c('gene','baseMean')]
DEG_Root_E$Root_E<-rep(1,dim(DEG_Root_E)[1])

DEG_Root_GE<-openxlsx::read.xlsx("Results/Salt_Root_DEG.xlsx",sheet = 3)
DEG_Root_GE<-DEG_Root_GE[,c('gene','baseMean')]
DEG_Root_GE$Root_GE<-rep(1,dim(DEG_Root_GE)[1])

GXERoot<-DEG_Root_GE$gene
GplusERoot<-setdiff(intersect(DEG_Root_G$gene,DEG_Root_E$gene),
                    GXERoot)
GRoot<-setdiff(DEG_Root_G$gene,GplusERoot)
ERoot<-setdiff(DEG_Root_E$gene,GplusERoot)
length(unique(c(GRoot,ERoot,GplusERoot,GXERoot)))

datR<-as.data.frame(cbind( Gene=as.vector(as.character(
  unique(c(GRoot,ERoot,GplusERoot,GXERoot))
)) ))

datR$Category<-ifelse(datR$Gene %in% GXERoot,"GXT" ,ifelse(
  datR$Gene %in% GplusERoot, "G+T", ifelse(
    datR$Gene %in% GRoot, "G","T"
  )
))


VennListRoot<-list(G=DEG_Root_G$gene,E=DEG_Root_E$gene,GXT=DEG_Root_GE$gene)
VennRoot<-plot(euler(VennListRoot),quantities=T)

######################### DAPC

ThemeDAPC<-theme(plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
                legend.title = element_text(size=10, face="bold"),
                legend.text = element_text(size=10),
                legend.position=c(0.5, 0.55),
                axis.title = element_text(size=10,face="bold",colour = "black"),
                panel.border = element_blank(), 
                #panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), 
                axis.line = element_line(colour = "black"),
                axis.text=element_text(size=10,vjust=1,face = "bold",colour = "black"),
                plot.margin = margin(0.15, 0.15, 0.2, 0.15, "cm")
                # Change margin here to add/remove extra margin
                #axis.ticks.x = element_blank()
)
##### Leaf
euclidean <- function(a, b) sqrt(sum((a - b)^2))
head(DesignL)
DesignL$Group<-paste(DesignL$Genotype,DesignL$Treatment,sep = ":")
DesignL$SampleID<-DesignL$LIBID

dapcL<-dapc(tdatL[,-c(ncol(tdatL))],DesignL$Group,
            var.contrib=F,scale=F,n.pca=5,n.da=5)
#scatter(dapcL,2,2, cell=0, pch=18:23, cstar=0, mstree=TRUE, lwd=2, lty=2)
#scatter(dapcL,2,1, cell=0, pch=18:23, cstar=0, mstree=TRUE, lwd=2, lty=2)
#head(dapcL$ind.coord)

datDAPCLeaf<-as.data.frame(dapcL$ind.coord)
datDAPCLeaf$SampleID<-rownames(datDAPCLeaf)
datDAPCLDesLeaf<-merge(DesignL,datDAPCLeaf,by="SampleID") %>% 
  as_tibble() %>% 
  mutate_at("Group", str_replace, "FIL", "Coastal") %>% 
  mutate_at("Group", str_replace, "HAL", "Inland") 

LeafPosLab<-as_tibble(as.data.frame(dapcL$grp.coord)) %>% 
  mutate(Group=rownames(dapcL$grp.coord)) %>% 
  mutate(Genotype = str_replace(Group, ":.*", "")) %>% 
  mutate_at("Group", str_replace, "FIL", "Coastal") %>% 
  mutate_at("Group", str_replace, "HAL", "Inland") 

euclidean(LeafPosLab[3,c(1:2)],LeafPosLab[4,c(1:2)])
euclidean(LeafPosLab[1,c(1:2)],LeafPosLab[2,c(1:2)])

LeafPlasticFIL<-filter(datDAPCLDesLeaf, Treatment=="Salinity" & Genotype=="FIL") %>% 
  dplyr::select(Genotype,LD1,LD2) %>% 
  rowwise() %>% 
  mutate(Plasticity= euclidean(LeafPosLab[1,c(1:2)],c(LD1,LD2)))

LeafPlasticHAL<-filter(datDAPCLDesLeaf, Treatment=="Salinity" & Genotype=="HAL") %>% 
  dplyr::select(Genotype,LD1,LD2) %>% 
  rowwise() %>% 
  mutate(Plasticity= euclidean(LeafPosLab[3,c(1:2)],c(LD1,LD2)))

LeafPlastic<-as.data.frame(rbind(LeafPlasticFIL,LeafPlasticHAL))
modleaf<-aov(Plasticity ~ Genotype, data = LeafPlastic)
summary(modleaf)

LeafDAPCPlot<-ggplot(datDAPCLDesLeaf)+
  geom_point(aes(x=LD1,y=LD2,
                  color=Group,fill=Group),alpha=0.75,size=5,show.legend = T)+ 
  xlim(range(datDAPCLDesLeaf$LD1)*1.25)+
  geom_segment(data=LeafPosLab,aes(x=LD1[1],xend=LD1[2],y=LD2[1],yend=LD2[2]),
               arrow = arrow(length = unit(0.3, "cm")),size=1,alpha=0.8,
               colour="#0073C2FF")+
  geom_segment(data=LeafPosLab,aes(x=LD1[3],xend=LD1[4],y=LD2[3],yend=LD2[4]),
               arrow = arrow(length = unit(0.3, "cm")),size=1,alpha=0.8,
               colour="#CD534CFF")+
  scale_color_manual(name = "Leaf Tissue Groups",
                     values = c("#0073C2FF", "#56B4E9",  "#CD534CFF","lightpink2"))+
  scale_fill_manual(name = "Leaf Tissue Groups",
                    values = c("#0073C2FF", "#56B4E9",  "#CD534CFF","lightpink2"))+
  labs(x="Linear Discriminent Function 1", y="Linear Discriminent Function 2")+
  # geom_text(data=LeafPosLab, aes(x=meanLD,y=meanDensity,label=Group,color=Group),
  #           fontface="bold",size=3,show.legend = F)+
  theme_classic()+ThemeDAPC

######### Root
head(DesignR)
DesignR$Group<-paste(DesignR$Genotype,DesignR$Treatment,sep = ":")
DesignR$SampleID<-DesignR$LIBID

dapcR<-dapc(tdatR[,-c(ncol(tdatR))],DesignR$Group,
            var.contrib=F,scale=F,n.pca=10,n.da=5)
#scatter(dapcR,2,2, cell=0, pch=18:23, cstar=0, mstree=TRUE, lwd=2, lty=2)
#scatter(dapcR,2,1, cell=0, pch=18:23, cstar=0, mstree=TRUE, lwd=2, lty=2)
#head(dapcR$ind.coord)

datDAPCRoot<-as.data.frame(dapcR$ind.coord)
datDAPCRoot$SampleID<-rownames(datDAPCRoot)
datDAPCLDesRoot<-merge(DesignR,datDAPCRoot,by="SampleID") %>% 
  as_tibble() %>% 
  mutate_at("Group", str_replace, "FIL", "Coastal") %>% 
  mutate_at("Group", str_replace, "HAL", "Inland") 


RootPosLab<-as_tibble(as.data.frame(dapcR$grp.coord)) %>% 
  mutate(Group=rownames(dapcL$grp.coord)) %>% 
  mutate(Genotype = str_replace(Group, ":.*", "")) %>% 
  mutate_at("Group", str_replace, "FIL", "Coastal") %>% 
  mutate_at("Group", str_replace, "HAL", "Inland") 



euclidean(RootPosLab[3,c(1:2)],RootPosLab[4,c(1:2)]) #Inland
euclidean(RootPosLab[1,c(1:2)],RootPosLab[2,c(1:2)])

RootPlasticFIL<-filter(datDAPCLDesRoot, Treatment=="Salinity" & Genotype=="FIL") %>% 
  dplyr::select(Genotype,LD1,LD2) %>% 
  rowwise() %>% 
  mutate(Plasticity= euclidean(RootPosLab[1,c(1:2)],c(LD1,LD2)))

RootPlasticHAL<-filter(datDAPCLDesRoot, Treatment=="Salinity" & Genotype=="HAL") %>% 
  dplyr::select(Genotype,LD1,LD2) %>% 
  rowwise() %>% 
  mutate(Plasticity= euclidean(RootPosLab[3,c(1:2)],c(LD1,LD2)))

RootPlastic<-as.data.frame(rbind(RootPlasticFIL,RootPlasticHAL))
modRoot<-aov(Plasticity ~ Genotype, data = RootPlastic)
summary(modRoot)

RootDAPCPlot<-ggplot(datDAPCLDesRoot)+
  geom_point(aes(x=LD1,y=LD2,
                 color=Group,fill=Group),alpha=0.75,size=5,shape=17,show.legend = T)+ 
  xlim(range(datDAPCLDesRoot$LD1)*1.25)+
  geom_segment(data=RootPosLab,aes(x=LD1[1],xend=LD1[2],y=LD2[1],yend=LD2[2]),
               arrow = arrow(length = unit(0.3, "cm")),size=1,alpha=0.8,
               colour="#0073C2FF")+
  geom_segment(data=RootPosLab,aes(x=LD1[3],xend=LD1[4],y=LD2[3],yend=LD2[4]),
               arrow = arrow(length = unit(0.3, "cm")),size=1,alpha=0.8,
               colour="#CD534CFF")+
  scale_color_manual(name = "Root Tissue Groups",
                     values = c("#0073C2FF", "#56B4E9",  "#CD534CFF","lightpink2"))+
  scale_fill_manual(name = "Root Tissue Groups",
                    values = c("#0073C2FF", "#56B4E9",  "#CD534CFF","lightpink2"))+
  labs(x="Linear Discriminent Function 1", y="Linear Discriminent Function 2")+
  # geom_text(data=RootPosLab, aes(x=meanLD,y=meanDensity,label=Group,color=Group),
  #           fontface="bold",size=3,show.legend = F)+
  theme_classic()+ThemeDAPC


ExpFig1<-ggdraw()+
  draw_plot(LeafDAPCPlot,x=0.0,y=0.5, width = 0.5, height = 0.5)+
  draw_label(x=0.03,y=0.97,label = "A)", color = "black", size = 12, fontface = "bold")+
  draw_plot(RootDAPCPlot,x=0.5,y=0.5, width = 0.5, height = 0.5)+
  draw_label(x=0.53,y=0.97,label = "B)", color = "black", size = 12, fontface = "bold")+
  draw_plot(VennLeaf,x=0.07,y=0.0, width = 0.4, height = 0.5)+
  draw_label(x=0.03,y=0.46,label = "C)", color = "black", size = 12, fontface = "bold")+
  draw_plot(VennRoot,x=0.57,y=0.0, width = 0.4, height = 0.5)+
  draw_label(x=0.53,y=0.46,label = "D)", color = "black", size = 12, fontface = "bold")

tiff("Plots/Figure_GeneExpression.tiff",width=7.5,height=6,units="in",res=300)
ExpFig1
dev.off()

