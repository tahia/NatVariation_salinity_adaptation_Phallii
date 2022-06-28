setwd("/home/taslima/data/JuengerLab/Research_Article_Preps/NatVariation_salinity_adaptation_Phallii")

library(xlsx)
library(tidyverse)
library(ggplot2)
library(qtl)
library(DescTools)

source("Scripts/RuntopGo.R")


QTLdat<-read.csv("Results/PHRIL_LeafIon_QTLStats.csv")

GFFAnnot<-read.delim(file="DBs/PhalliiHAL_496_v2.1.Annotation.tab",
                     sep = "\t")

#### FOR GO
genedb<-read.delim("DBs/PhalliiHAL_496_v2.1.Onlygene.mod.gff3",
                   header = F,sep = "\t")
geneID2GO <- "DBs/PH_HAL_v2_GO.tab"
info<-read.delim("DBs/PhalliiHAL_496_v2.1.annotation_info.txt",na.strings = "")
head(info)
#Select only columns for geneID and GO
info<-info[,c(2,10,16)]
#Get only genes that have GO annotation
length(which(!is.na(info$GO)))
info<-info[which(!is.na(info$GO)),]


QTLShort<-QTLdat[-which(is.na(QTLdat$qtlnames)),
                 c('qtlnames','chr','lowmarker','highmarker','lowposition','highposition','Trait')]
QTLShort$chr<-paste0("Chr0",QTLShort$chr)
QTLShort$Trait<-as.character(QTLShort$Trait)
#Change trait name of LNa/K as LNabyK to write as sheetname
QTLShort$Trait[which(QTLShort$Trait=="LNa/K")]<-"LNabyK"

GeneinQTLInterval<-function(GFFAnnot=NULL,qtl.name=NULL,chr=NULL, 
                            lowpos=NULL, highpos=NULL,trait=NULL) {
  if(any(is.null(c(chr,lowpos,highpos))))
    stop("Chromosome and QTL interval is missing\n")
  if(is.null(qtl.name))
    stop("QTL name is missing\n")
  
  QAnnot<-GFFAnnot[with(GFFAnnot,
                        Chr == chr &
                          Start>=lowpos &
                          End<=highpos),]
  return(QAnnot)
}

for (i in 1:nrow(QTLShort)){
  Gendat<-GeneinQTLInterval(GFFAnnot = GFFAnnot,
                            qtl.name = QTLShort$qtlnames[i],
                            chr = QTLShort$chr[i], 
                            lowpos = QTLShort$lowmarker[i],
                            highpos = QTLShort$highmarker[i],
                            trait = QTLShort$Trait[i])
  name<-paste(QTLShort$Trait[i],QTLShort$qtlnames[i],sep = "_")
  
  geneTab<-as.data.frame(cbind(GeneID=Gendat$GeneID))
  geneTabAnn<-merge(geneTab, info,by.x="GeneID",by.y="locusName",all.x=T,all.y=F)
  SigGO<-RunTopGO(geneID2GO,info,Gendat$GeneID)
  if (nrow(SigGO)==0) {
    SigGO[1,]<-rep(NA,8)
  }
  if (i==1) {
    write.xlsx(Gendat,file="Results/PHRIL_LeafIon_CandidateGenes.xlsx",sheetName = name,append = F)
    write.xlsx(SigGO,file="Results/PHRIL_LeafIon_CandidateGenesGOEnrichment.xlsx",sheetName = name,append = F)
  } else {
    write.xlsx(Gendat,file="Results/PHRIL_LeafIon_CandidateGenes.xlsx",sheetName = name,append = T)
    write.xlsx(SigGO,file="Results/PHRIL_LeafIon_CandidateGenesGOEnrichment.xlsx",sheetName = name,append = T)
  }
}

############## TARGET GENE Enrichment
ChrLength<-as.data.frame(cbind(Chr=paste0("0",seq(1:9)),
                               Length=c(109.8,115.3,100.5,
                                        83.3,110.4,72.5,
                                        62.3,72.2,137.7)
))
  
# cross<-read.cross(dir = "Data",file = "PHSALT_RQTL_Constitutive_Mean.csv",format = "csv")
# cross<-convert2riself(cross)

load("Data/PHSALT_RQTL_Constitutive_Mean.RData")

### HKT
HKT<-c("PhHAL.1G054100","PhHAL.4G025100","PhHAL.5G395500","PhHAL.7G275900","PhHAL.7G276000")
GFFAnnotHKT<-GFFAnnot[which(GFFAnnot$GeneID %in% HKT),]

head(QTLShort)

HKTPval<-c()

for (i in 1:nrow(QTLShort)){
  cat("Now running QTL:", QTLShort$chr[i], "for trait",QTLShort$Trait[i], format(Sys.time(), "%a %b %d %X %Y"),"......\n" )
  Gendat<-GeneinQTLInterval(GFFAnnot = GFFAnnotHKT,
                            qtl.name = QTLShort$qtlnames[i],
                            chr = QTLShort$chr[i], 
                            lowpos = QTLShort$lowmarker[i],
                            highpos = QTLShort$highmarker[i],
                            trait = QTLShort$Trait[i])
  name<-paste(QTLShort$Trait[i],QTLShort$qtlnames[i],sep = "_")
  Linterval<-QTLShort$highposition[i]-QTLShort$lowposition[i]
  if (nrow(Gendat)>=1) {
    Boot<-c()
    set.seed(123456)
    for (j in 1:1000){
      Chr<-as.character(sample(ChrLength$Chr,1))
      CLength<-as.numeric(as.character(ChrLength$Length[which(ChrLength$Chr==Chr)]))
      Start<-sample(seq(1:(CLength-Linterval )),1)
      End<-Start+Linterval
      
      if ( (paste0("Chr",Chr)==QTLShort$chr[i]) & 
           (Overlap(c(QTLShort$lowmarker[i]: QTLShort$highmarker[i]),
                    c(Start:End))>0 ) 
           ) {
        j=j-1
        next
      } else {
      
      lowpos<-as.numeric(strsplit(find.marker(cross,Chr ,
                                              Start), "_")[[1]][2])
      highpos<-as.numeric(strsplit(find.marker(cross,Chr ,
                                               End), "_")[[1]][2])
      
      #cat("j=", j, Chr,lowpos, highpos,"\n")
      BCount<-GeneinQTLInterval(GFFAnnot = GFFAnnotHKT,
                        qtl.name = QTLShort$qtlnames[i],
                        chr = paste0("Chr",Chr) ,
                        lowpos = lowpos,
                        highpos = highpos,
                        trait = QTLShort$Trait[i])
      Boot<-append(Boot,nrow(BCount))
      }
    }
    pval<-length(which(Boot>= nrow(Gendat)))/1000
    HKTPval<-append(HKTPval,pval)
  } else {
    HKTPval<-append(HKTPval,NA)
  }
  
}

### SOS
SOS<-c("PhHAL.3G490000","PhHAL.4G089400",
       "PhHAL.1G043900","PhHAL.1G130000","PhHAL.2G474900",
       "PhHAL.2G475000","PhHAL.3G008400","PhHAL.3G190900",
       "PhHAL.3G231300","PhHAL.3G333600","PhHAL.3G335800",
       "PhHAL.5G110400","PhHAL.5G256100","PhHAL.6G205900",
       "PhHAL.8G009500","PhHAL.8G120300","PhHAL.9G484500")
GFFAnnotSOS<-GFFAnnot[which(GFFAnnot$GeneID %in% SOS),]

head(QTLShort)

SOSPval<-c()

for (i in 1:nrow(QTLShort)){
  cat("Now running QTL:", QTLShort$chr[i], "for trait",QTLShort$Trait[i], format(Sys.time(), "%a %b %d %X %Y"),"......\n" )
  Gendat<-GeneinQTLInterval(GFFAnnot = GFFAnnotSOS,
                            qtl.name = QTLShort$qtlnames[i],
                            chr = QTLShort$chr[i], 
                            lowpos = QTLShort$lowmarker[i],
                            highpos = QTLShort$highmarker[i],
                            trait = QTLShort$Trait[i])
  name<-paste(QTLShort$Trait[i],QTLShort$qtlnames[i],sep = "_")
  Linterval<-QTLShort$highposition[i]-QTLShort$lowposition[i]
  if (nrow(Gendat)>=1) {
    Boot<-c()
    set.seed(123456)
    for (j in 1:1000){
      Chr<-as.character(sample(ChrLength$Chr,1))
      CLength<-as.numeric(as.character(ChrLength$Length[which(ChrLength$Chr==Chr)]))
      Start<-sample(seq(1:(CLength-Linterval )),1)
      End<-Start+Linterval
      
      if ( (paste0("Chr",Chr)==QTLShort$chr[i]) & 
           (Overlap(c(QTLShort$lowmarker[i]: QTLShort$highmarker[i]),
                    c(Start:End))>0 ) 
      ) {
        j=j-1
        next
      } else {
      
      lowpos<-as.numeric(strsplit(find.marker(cross,Chr ,
                                              Start), "_")[[1]][2])
      highpos<-as.numeric(strsplit(find.marker(cross,Chr ,
                                               End), "_")[[1]][2])
      
      #cat("j=", j, Chr,lowpos, highpos,"\n")
      BCount<-GeneinQTLInterval(GFFAnnot = GFFAnnotSOS,
                                qtl.name = QTLShort$qtlnames[i],
                                chr = paste0("Chr",Chr) ,
                                lowpos = lowpos,
                                highpos = highpos,
                                trait = QTLShort$Trait[i])
      Boot<-append(Boot,nrow(BCount))
    }
    }
    pval<-length(which(Boot>= nrow(Gendat)))/1000
    SOSPval<-append(SOSPval,pval)
  } else {
    SOSPval<-append(SOSPval,NA)
  }
  
}

### HAK/KT/KUP
HAK<-c("PhHAL.1G242900", "PhHAL.2G001500","PhHAL.2G119400",
       "PhHAL.2G274700", "PhHAL.2G373800", "PhHAL.2G471600",
       "PhHAL.2G475100", "PhHAL.4G051700" ,"PhHAL.4G171100",
        "PhHAL.5G026500" , "PhHAL.5G029000 ", "PhHAL.5G030400",
       "PhHAL.5G339300", "PhHAL.6G064700", "PhHAL.6G249100",
        "PhHAL.7G124900", "PhHAL.7G125100", "PhHAL.7G125200",
        "PhHAL.7G278700", "PhHAL.7G281600",  "PhHAL.9G174200",
        "PhHAL.9G174300" , "PhHAL.9G174400" , "PhHAL.9G174500",
       "PhHAL.9G484900")
GFFAnnotHAK<-GFFAnnot[which(GFFAnnot$GeneID %in% HAK),]

head(QTLShort)

HAKPval<-c()

for (i in 1:nrow(QTLShort)){
  cat("Now running QTL:", QTLShort$chr[i], "for trait",QTLShort$Trait[i], format(Sys.time(), "%a %b %d %X %Y"),"......\n" )
  Gendat<-GeneinQTLInterval(GFFAnnot = GFFAnnotHAK,
                            qtl.name = QTLShort$qtlnames[i],
                            chr = QTLShort$chr[i], 
                            lowpos = QTLShort$lowmarker[i],
                            highpos = QTLShort$highmarker[i],
                            trait = QTLShort$Trait[i])
  name<-paste(QTLShort$Trait[i],QTLShort$qtlnames[i],sep = "_")
  Linterval<-QTLShort$highposition[i]-QTLShort$lowposition[i]
  if (nrow(Gendat)>=1) {
    Boot<-c()
    set.seed(123456)
    for (j in 1:1000){
      Chr<-as.character(sample(ChrLength$Chr,1))
      CLength<-as.numeric(as.character(ChrLength$Length[which(ChrLength$Chr==Chr)]))
      Start<-sample(seq(1:(CLength-Linterval )),1)
      End<-Start+Linterval
      
      if ( (paste0("Chr",Chr)==QTLShort$chr[i]) & 
           (Overlap(c(QTLShort$lowmarker[i]: QTLShort$highmarker[i]),
                    c(Start:End))>0 ) 
      ) {
        j=j-1
        next
      } else {
      
      lowpos<-as.numeric(strsplit(find.marker(cross,Chr ,
                                              Start), "_")[[1]][2])
      highpos<-as.numeric(strsplit(find.marker(cross,Chr ,
                                               End), "_")[[1]][2])
      
      #cat("j=", j, Chr,lowpos, highpos,"\n")
      BCount<-GeneinQTLInterval(GFFAnnot = GFFAnnotHAK,
                                qtl.name = QTLShort$qtlnames[i],
                                chr = paste0("Chr",Chr) ,
                                lowpos = lowpos,
                                highpos = highpos,
                                trait = QTLShort$Trait[i])
      Boot<-append(Boot,nrow(BCount))
      }
    }
    pval<-length(which(Boot>= nrow(Gendat)))/1000
    HAKPval<-append(HAKPval,pval)
  } else {
    HAKPval<-append(HAKPval,NA)
  }
  
}


QTLShort$PvalHKT<-HKTPval
QTLShort$PvalSOS<-SOSPval
QTLShort$PvalHAK<-HAKPval

write.csv(QTLShort,"Results/PHRIL_LeafIon_Reference_CandidateGenes_Bootstrapped.csv",row.names = F)
