library(tidyverse)
library(cowplot)
library(xlsx)

Map<-read.csv("Data/PHRIL2019_gmap.csv")
Map$chr<-as.numeric(as.character(sub("Chr0","",Map$chr)))

######## Constitutive
QTLCons<-read.xlsx("Results/PHRIL_QTLEffect_ALL.xlsx",sheetIndex = 1)
colnames(QTLCons)

QTLCons$Chromosome<-as.numeric(QTLCons$Chromosome)
QTLCons$lod <- cut(round(QTLCons$LOD),breaks=c(2,6,10,15),labels=c("2 - 6","6 - 10",
                                    "10 - 15"),right = F)
QTLCons$lod <- factor(QTLCons$lod,levels=c("2 - 6","6 - 10",
                                           "10 - 15"))
QTLCons$PAllele<-QTLCons$PositiveAllele
QTLCons$NAllele<-ifelse( QTLCons$PositiveAllele == "Coastal", "Inland" , "Coastal")
QTLCons$Phenotypes<-paste(QTLCons$Trait,"Constitutive",sep = " ")

(p1<-as.tibble(Map) %>% group_by(chr) %>% summarise(length=max(pos,na.rm = T)) %>% 
    ggplot(aes(x=chr))+geom_segment(aes(x=chr, xend=chr,y=0,yend=length),size=0.5)+
    geom_segment(data  = Map,aes(x=(chr-0.075),xend=(chr+0.075),y=pos,yend=pos),size=0.5)+
    geom_linerange(data=QTLCons,aes(x=ChrProxy,y=QTL_peak,ymin=UPOS,ymax=DPOS,color=Phenotypes,size=lod),
                   alpha=0.95,show.legend = F)+
    geom_linerange(data=QTLCons,aes(x=ChrProxy,ymin=QTL_peak,ymax=QTL_peak+1),
                   alpha=1,color="black",size=4)+
    geom_point(data=QTLCons,aes(x=ChrProxy,y=DPOS+2,fill=PAllele),
               size=3,alpha=0.95,shape=24)+
    # geom_point(data=QTLCons,aes(x=ChrProxy,y=DPOS+2.0,fill=NAllele),
    #            size=3,alpha=0.95,shape=25,show.legend = F)+
    #geom_point(data=QTL,aes(x=Chromosome+0.3,y=QTL_peak,shape=Allele,fill=Allele),
    #              size=4,alpha=0.95)+
    #geom_point(data = Effect,aes(x=Chromosome+0.2+FIL,y=Position),size=0.5)+
    #geom_linerange(data=CNLD,aes(x=chr+0.1,y=pos,ymin=lowposition, ymax=highposition),color="black",size=1.5,linetype=1)+
    #geom_point(data=QTL,aes(x=Chromosome+0.2,y=QTL$QTL.peak))+ #,position = position_dodge(width = 0.35)
    scale_x_continuous(breaks = seq(1,9,1),labels = paste("Chr",seq(1,9,1),sep = ""))+
    scale_y_continuous(breaks = seq(0,150,50),labels = seq(0,150,50),limits = c(0,150))+
    scale_shape_manual(values = c(
      "Inland"=21,
      "Coastal"=23
    ))+
    scale_fill_manual(values = c(
      Coastal="#0072B2",
      Inland="#D55E00"
    ))+
    #scale_color_manual(values = viridis(6))+
    scale_color_manual(values = c(
      "AGB Constitutive"="#009E73",
      "BGB Constitutive"= "#E6AB02",
      "RB Constitutive"= "gray50"
      # "LRWC_Constitutive"= "#440154FF",
      # "LC_Constitutive" = "#009E73",
      # "LC_Responsive"="gray60"
    ))+
    scale_size_manual(values = c(
      "2 - 6" = 1,
      "6 - 10"= 2,
      "10 - 15" = 3
    ))+
    # annotate(geom="text", x=1, y=150, label="A.",size=8,
    #          color="black")+
    labs(x="",y="Genetic Distance (cM)",color="Phenotypes",size="LOD Score",fill="Allele")+
    guides(color = F,shape=F,fill=F)+
    theme(
      panel.grid  = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      axis.title = element_text(size=14,face="bold"),
      axis.text = element_text(face="bold",size=12),
      axis.line.y = element_line(size=0.5),
      axis.line.x= element_blank(),
      axis.ticks =element_blank(),
      legend.text  = element_text(face="bold",size=11),
      legend.title = element_text(face = "bold",size = 12),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.key.size = unit(0.5, "cm"),
      plot.margin = unit(c(1,3,1,1),"cm"),
      legend.position=c(1.05,0.12))
  
)

leg1<-as.tibble(Map) %>% group_by(chr) %>% summarise(length=max(pos,na.rm = T)) %>% 
  ggplot(aes(x=chr))+geom_segment(aes(x=chr, xend=chr,y=0,yend=length),size=0.5)+
  geom_segment(data  = Map,aes(x=(chr-0.1),xend=(chr+0.1),y=pos,yend=pos),size=0.5)+
  geom_linerange(data=QTLCons,aes(x=Chromosome+0.2,ymin=UPOS,ymax=DPOS,color=Phenotypes,size=lod),
                 alpha=0.95,position = position_dodge(width = 0.35),size=3)+
  scale_x_continuous(breaks = seq(1,9,1),labels = paste("Chr",seq(1,9,1),sep = ""))+
  scale_y_continuous(breaks = seq(0,150,50),labels = seq(0,150,50),limits = c(0,150))+
  #scale_color_manual(values = viridis(6))+
  scale_color_manual(values = c(
    "AGB Constitutive"="#009E73",
    "BGB Constitutive"= "orange",
    "RB Constitutive"= "gray60"
    # "LRWC_Constitutive"= "#440154FF",
    # "LC_Constitutive" = "#009E73",
    # "LC_Responsive"="gray60"
  ))+
  labs(x="Chromosme",y="distance (cM)",color="QTL",size="LOD Score")+
  guides(size = FALSE)+
  theme(
    panel.grid  = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(size=14,face="bold"),
    axis.text = element_text(face="bold",size=12),
    axis.line.y = element_line(size=0.5),
    axis.line.x= element_blank(),
    legend.text  = element_text(face="bold",size=11),
    legend.title = element_text(face = "bold",size = 12),
    legend.key = element_rect(colour = NA, fill = NA),
    legend.key.size = unit(0.8, "cm"))

######Responsive and Solitary

QTLRes<-read.xlsx("Results/PHRIL_QTLEffect_ALL.xlsx",sheetIndex = 2)
colnames(QTLRes)
QTLRes$Phenotypes<-paste(QTLRes$Trait,"Responsive",sep = " ")

QTLSol<-read.xlsx("Results/PHRIL_QTLEffect_ALL.xlsx",sheetIndex = 3)
colnames(QTLSol)
#QTLSol$Phenotypes<-paste(QTLSol$Trait,"Ionic",sep = " ")
QTLSol$Phenotypes<-QTLSol$Trait

QTLSalinity<-rbind(QTLRes,QTLSol)

QTLSalinity$Chromosome<-as.numeric(QTLSalinity$Chromosome)
QTLSalinity$lod <- cut(round(QTLSalinity$LOD),breaks=c(2,6,10,15),labels=c("2 - 6","6 - 10",
                                                                          "10 - 15"),right = F)
QTLSalinity$lod <- factor(QTLSalinity$lod,levels=c("2 - 6","6 - 10",
                                                   "10 - 15"))
QTLSalinity$PAllele<-QTLSalinity$PositiveAllele
QTLSalinity$NAllele<-ifelse( QTLSalinity$PositiveAllele == "Coastal", "Inland" , "Coastal")



(p2<-as.tibble(Map) %>% group_by(chr) %>% summarise(length=max(pos,na.rm = T)) %>% 
    ggplot(aes(x=chr))+geom_segment(aes(x=chr, xend=chr,y=0,yend=length),size=0.5)+
    geom_segment(data  = Map,aes(x=(chr-0.075),xend=(chr+0.075),y=pos,yend=pos),size=0.5)+
    geom_linerange(data=QTLSalinity,aes(x=ChrProxy,y=QTL_peak,ymin=UPOS,ymax=DPOS,color=Phenotypes,size=lod),
                   alpha=0.95)+
    geom_linerange(data=QTLSalinity,aes(x=ChrProxy,ymin=QTL_peak,ymax=QTL_peak+1),
                   alpha=1,color="black",size=4)+
    geom_point(data=QTLSalinity,aes(x=ChrProxy,y=DPOS+2,fill=PAllele),
               size=3,alpha=0.95,shape=24)+
    # geom_point(data=QTLSalinity,aes(x=ChrProxy,y=DPOS+2.0,fill=NAllele),
    #            size=3,alpha=0.95,shape=25,show.legend = F)+
    #geom_point(data=QTL,aes(x=Chromosome+0.3,y=QTL_peak,shape=Allele,fill=Allele),
    #              size=4,alpha=0.95)+
    #geom_point(data = Effect,aes(x=Chromosome+0.2+FIL,y=Position),size=0.5)+
    #geom_linerange(data=CNLD,aes(x=chr+0.1,y=pos,ymin=lowposition, ymax=highposition),color="black",size=1.5,linetype=1)+
    #geom_point(data=QTL,aes(x=Chromosome+0.2,y=QTL$QTL.peak))+ #,position = position_dodge(width = 0.35)
    scale_x_continuous(breaks = seq(1,9,1),labels = paste("Chr",seq(1,9,1),sep = ""))+
    scale_y_continuous(breaks = seq(0,150,50),labels = seq(0,150,50),limits = c(0,150))+
    scale_shape_manual(values = c(
      "Inland"=21,
      "Coastal"=23
    ))+
    scale_fill_manual(values = c(
      Coastal="#0072B2",
      Inland="#D55E00"
    ))+
    #scale_color_manual(values = viridis(6))+
    scale_color_manual(values = c(
      "BGB Responsive"=  "orange" ,
      "K"= "dodgerblue",
      "Na"= "firebrick2",
      "Na/K" = "#7570B3"
    ))+
    scale_size_manual(values = c(
      "2 - 6" = 1,
      "6 - 10"= 2,
      "10 - 15" = 3
    ))+ 
    # annotate(geom="text", x=1, y=150, label="B.",size=8,
    #              color="black")+
    labs(x="",y="Genetic Distance (cM)",color="Phenotypes",size="LOD Score",fill="Allele")+
    guides(color = F)+
    theme(
      panel.grid  = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      axis.title = element_text(size=14,face="bold"),
      axis.text = element_text(face="bold",size=12),
      axis.line.y = element_line(size=0.5),
      axis.line.x= element_blank(),
      axis.ticks =element_blank(),
      legend.text  = element_text(face="bold",size=11),
      legend.title = element_text(face = "bold",size = 12),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.key.size = unit(0.5, "cm"),
      plot.margin = unit(c(1,3,1,1),"cm"),
      legend.position=c(1.02,0.8))
  
)

leg2<-as.tibble(Map) %>% group_by(chr) %>% summarise(length=max(pos,na.rm = T)) %>% 
  ggplot(aes(x=chr))+geom_segment(aes(x=chr, xend=chr,y=0,yend=length),size=0.5)+
  geom_segment(data  = Map,aes(x=(chr-0.1),xend=(chr+0.1),y=pos,yend=pos),size=0.5)+
  geom_linerange(data=QTLSalinity,aes(x=Chromosome+0.2,ymin=UPOS,ymax=DPOS,color=Phenotypes,size=lod),
                 alpha=0.95,position = position_dodge(width = 0.35),size=3)+
  scale_x_continuous(breaks = seq(1,9,1),labels = paste("Chr",seq(1,9,1),sep = ""))+
  scale_y_continuous(breaks = seq(0,150,50),labels = seq(0,150,50),limits = c(0,150))+
  scale_color_manual(values = c(
    "BGB Responsive"=  "orange" ,
    "K"= "dodgerblue",
    "Na"= "firebrick2",
    "Na/K" = "#7570B3"
  ))+
  labs(x="Chromosme",y="distance (cM)",color="QTL",size="LOD Score")+
  guides(size = FALSE)+
  theme(
    panel.grid  = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(size=14,face="bold"),
    axis.text = element_text(face="bold",size=12),
    axis.line.y = element_line(size=0.5),
    axis.line.x= element_blank(),
    legend.text  = element_text(face="bold",size=11),
    legend.title = element_text(face = "bold",size = 12),
    legend.key = element_rect(colour = NA, fill = NA),
    legend.key.size = unit(0.8, "cm"))


# get_legend<-function(myggplot){
#   tmp <- ggplot_gtable(ggplot_build(myggplot))
#   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#   legend <- tmp$grobs[[leg]]
#   return(legend)
# }


lg1 <- get_legend(leg1)
lg2 <- get_legend(leg2)

library(cowplot)  
# plot.this <-
#   ggdraw() +
#   draw_plot(p1,x = 0, y = 0,width = 0.55,height = 0.9) +
#   draw_plot(lg1, x = 0.12, y = .75, width = .3, height = .3)+
#   draw_plot(p2,x = 0.47, y = 0,width = 0.55,height = 0.9) +
#   draw_plot(lg2, x = 0.6, y = .73, width = .3, height = .3)

plot.this <-
  ggdraw() +
  draw_plot(p1,x = 0, y = 0.49,width = 1,height = 0.55) +
  draw_label(x=0.03,y=0.97,label = "A)", color = "black", size = 12, fontface = "bold")+
  draw_plot(lg1, x = 0.49, y = .75, width = .3, height = .3)+
  draw_plot(p2,x = 0, y = -0.01,width = 1.05,height = 0.55) +
  draw_label(x=0.03,y=0.48,label = "B)", color = "black", size = 12, fontface = "bold")+
  draw_plot(lg2, x = 0.5, y = .25, width = .3, height = .3)

tiff("Plots/Figure_ALLQTL_HR.tiff",width=8,height=10,units="in",res=300)
plot.this
dev.off()

tiff("Plots/Figure_ALLQTL_LR.tiff",width=8,height=10,units="in",res=100)
plot.this
dev.off()


