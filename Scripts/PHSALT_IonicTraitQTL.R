##### RUN QTL ANALYSIS for IONIC trait

library(qtl)
library(qtl2)
library(qtlTools)

cross<-read.cross(format = "csvs",genfile = "Data/PHSALT_ION_RQTL_gen.csv",
           phefile = "Data/PHSALT_ION_RQTL_phe.csv",genotypes = c("AA","BB"))

cross<-convert2riself(cross)
summary(cross)

cross <- calc.genoprob(cross, step=2,error.prob = 0.002)
cross<-sim.geno(cross)

##### This permutation in computationally heavy therefore required multiple cores
cores=8
perm=1000

set.seed(54955149)
operm1.sc <- scanone(cross, method="hk",pheno.col = c(2,5,6), n.perm=perm,
                     n.cluster = cores)

sc1 <- scanone(cross, method="hk",pheno.col = c(2,5,6),
               n.cluster = cores)

summary(sc1,perms=operm1.sc, alpha=0.05)

cat(paste(c("Start building null models at",format(Sys.time(), "%a %b %d %X %Y"),"\n"),sep = " "))
set.seed(54955149)
operm2.null.Ion <- scantwo(cross, method="hk",pheno.col = c(2,5,6), n.perm=perm,
                     n.cluster = cores)

sc1 <- scanone(cross, method="hk",pheno.col = c(2,5,6),
               n.cluster = cores)

summary(sc1,perms=operm2.null.Ion, alpha=0.05)

############# LK
print(penLK <- calc.penalties(operm2.null.Ion[,1],alpha = 0.1))

outsw1LK <- stepwiseqtl(cross, max.qtl=6,method="hk",pheno.col =2,
                        penalties=penLK,verbose=FALSE,keeplodprofile=TRUE,
                        keeptrace=TRUE)
plotLodProfile(outsw1LK)

LKqtl<-makeqtl(cross,chr = outsw1LK$chr,
               pos = outsw1LK$pos, qtl.name = outsw1LK$name ,what = "prob")


modLK<-fitqtl(cross , pheno.col = 2,qtl=LKqtl,
              formula=attributes(outsw1LK)$formula)

summary(modLK)
qqPlot(as.vector(attr(modLK, "residuals")))
shapiro.test(as.vector(attr(modLK, "residuals")))

LKdat<-qtlStats(cross,pheno.col = 2,outsw1LK,calcConfint = T,expandtomarkers=T)
LKdat$Trait<-rep("LK",nrow(LKdat))


############# TLNa
print(penTLNa <- calc.penalties(operm2.null.Ion[,2],alpha = 0.1))

outsw1TLNa <- stepwiseqtl(cross, max.qtl=6,method="hk",pheno.col =5,
                          penalties=penTLNa,verbose=FALSE,keeplodprofile=TRUE,
                          keeptrace=TRUE)
plotLodProfile(outsw1TLNa)

TLNaqtl<-makeqtl(cross,chr = outsw1TLNa$chr,
                 pos = outsw1TLNa$pos, qtl.name = outsw1TLNa$name,
                 what = "prob")

summary(fitqtl(cross , pheno.col = 5,qtl=TLNaqtl,
               formula=attributes(outsw1TLNa)$formula,get.ests = T),pvalues=T)

modTLNa<-fitqtl(cross , pheno.col = 5,qtl=TLNaqtl,
                formula=attributes(outsw1TLNa)$formula)

summary(modTLNa)
qqPlot(as.vector(attr(modTLNa, "residuals")))
shapiro.test(as.vector(attr(modTLNa, "residuals")))

LNadat<-qtlStats(cross,pheno.col = 5,outsw1TLNa,calcConfint = T,expandtomarkers=T)
LNadat$Trait<-rep("LNa",nrow(LNadat))


########## TLNabyK
print(penLNabyK <- calc.penalties(operm2.null.Ion[,3],alpha = 0.1))

outsw1TLNabyK <- stepwiseqtl(cross, max.qtl=6,method="hk",pheno.col = 6,
                             penalties=penLNabyK,verbose=FALSE,keeplodprofile=TRUE,
                             keeptrace=TRUE)

plotLodProfile(outsw1TLNabyK)
summary(outsw1TLNabyK)

TLNabyKqtl<-makeqtl(cross,chr = outsw1TLNabyK$chr,
                    pos = outsw1TLNabyK$pos, qtl.name = outsw1TLNabyK$name,
                    what = "prob")


modTLNabyK<-fitqtl(cross , pheno.col = 6,qtl=TLNabyKqtl,
                   formula=attributes(outsw1TLNabyK2)$formula,get.ests = T)

summary(modTLNabyK)
qqPlot(as.vector(attr(modTLNabyK, "residuals")))
shapiro.test(as.vector(attr(modTLNabyK, "residuals")))

LNabyKdat<-qtlStats(cross,pheno.col = 6,outsw1TLNabyK2,calcConfint = T,expandtomarkers=T)
LNabyKdat$Trait<-rep("LNa/K",nrow(LNabyKdat))

QTLdat<-rbind(LKdat,LNadat,LNabyKdat)
QTLdat$highmarker<-as.numeric(as.character(gsub("Chr0[1-9]_","",QTLdat$highmarker)))
QTLdat$lowmarker<-as.numeric(as.character(gsub("Chr0[1-9]_","",QTLdat$lowmarker)))

write.csv(QTLdat,"Results/PHRIL_LeafIon_QTLStats.csv",row.names = F)