library(DESeq2)
library(ggplot2)
library("RColorBrewer")
library(devtools)
library("BiocParallel")
library("vsn")
library(ade4)
library(factoextra)


RunPrettyPCA<-function(Design=NULL,dat=NULL,type=NULL,collist=NULL,title=NULL) {
  if(is.null(collist)) {collist<- c("#0073C2FF",   "#CD534CFF")}
  if(is.null(Design)) {stop("A Design Matrix is required\n")}
  if(is.null(dat)) {stop("A data matrix is required\n")}
  if(is.null(type)) {type="Regular"}
  if(is.null(title)) {title=""}

  #Get the variances of each gene to select top 1000 genes for PCA
  var<-as.data.frame(apply(dat,1,var))
  colnames(var)<-"GeneVariance"
  #TopVarGene<-rownames(var)[order(var$GeneVariance,decreasing = T)][1:1000]
  TopVarGene<-rownames(var)[order(var$GeneVariance,decreasing = T)][1:nrow(var)]
  dat_sel<-dat[which(rownames(dat) %in% TopVarGene  ),]


  library(FactoMineR)
  ExpdatPrim<-as.data.frame(na.omit(t(dat_sel)))
  ExpdatPrim$ID<-rownames(ExpdatPrim)
  DesignTrim<-Design[,c('Genotype', 'Treatment')]
  DesignTrim$ID<-rownames(DesignTrim)
  #colnames(DesignTrim)<-c("Ecotype","Treatment","ID")
  Expdat<-merge(DesignTrim,ExpdatPrim,by="ID")
  Expdat<-Expdat[,-c(1)]
  
  axis1=1
  axis2=2
  shape=2
  col.ind=1
  axis<-c(axis1,axis2)
  
  seldat<-Expdat[,-c(shape,col.ind)]
  res.pca <- PCA(seldat, graph = FALSE,ncp = 5)
  perVar <- res.pca$eig[,2]
  ind <- data.frame(res.pca$ind$coord[, axis, drop = FALSE])
  var<-as.data.frame(res.pca$var$coord)
  
  ind$shape<-as.vector(Expdat[,shape])
  ind$col<-as.vector(Expdat[,col.ind])
  xlab=paste("Axis",axis1," (",round(as.numeric(perVar[axis1]),2),"%)",sep = "")
  ylab=paste("Axis",axis2," (",round(as.numeric(perVar[axis2]),2),"%)",sep = "")
  
  ### Set the scale manually
  fivenum(ind$Dim.1)
  fivenum(ind$Dim.2)
  
  
  ## Need to work on shape a bit
  shapelist<-c(21,24)


########## Type:Regular

  if (type=="Regular") {
    #title<-paste("PCA on Axis ",axis1, " and ",axis2,sep = "")
  
    plot<-ggplot(ind)+geom_point(aes(x=ind[,1],y=ind[,2],shape=ind$shape,fill=ind$col),size=3,stroke=1.05, alpha=0.95,show.legend = F)+
      geom_hline(yintercept = 0,linetype=2)+geom_vline(xintercept = 0,linetype=2)+
      theme_bw() + labs(x =xlab, y=ylab,title = title,shape=as.name(colnames(Expdat)[shape]),fill=as.name(colnames(Expdat)[col.ind]))+ 
      scale_shape_manual(values=shapelist)+
      scale_fill_manual(values=collist)+
      theme(plot.title = element_text(hjust = 0.5,size=12, face = "bold"),
            legend.title = element_text(size=10, face="bold"),
            legend.text = element_text(size=10),
            legend.position=c(1.05, 0.75),
            axis.title = element_text(size=10,face="bold"),
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.text=element_text(size=10,vjust=1,face = "bold",colour = "black"),
            plot.margin = margin(0.15, 0.15, 0.15, 0.15, "cm") # Change margin here to add/remove extra margin
            #axis.ticks.x = element_blank()
      )
  }

########## Type: Shape

  if (type=="Shape") {
    title<-paste("PCA on Axis ",axis1, " and ",axis2,sep = "")
    plot<-ggplot(ind)+geom_point(aes(x=ind[,1],y=ind[,2],shape=ind$shape),size=3, stroke=1.25,alpha=1,show.legend = T)+
      geom_hline(yintercept = 0,linetype=2)+geom_vline(xintercept = 0,linetype=2)+
      theme_bw() + labs(x =xlab, y=ylab,title = title,shape=as.name(colnames(Expdat)[shape]))+ 
      scale_shape_manual(values=shapelist)+
      theme(plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
            legend.title = element_text(size=14, face="bold"),
            legend.text = element_text(size=12),
            legend.position=c(1.05, 0.75),
            axis.title = element_text(size=12,face="bold"),
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.text=element_text(size=12,vjust=1),
            plot.margin = margin(0.25, 5, 0.25, 0.25, "cm")
            #axis.ticks.x = element_blank()
      )
    
  }

########## Type:Color

  if (type=="Color") {
    title<-paste("PCA on Axis ",axis1, " and ",axis2,sep = "")
    plot<-ggplot(ind)+geom_point(aes(x=ind[,1],y=ind[,2],color=ind$col),size=5,alpha=0.95,show.legend = T)+
      geom_hline(yintercept = 0,linetype=2)+geom_vline(xintercept = 0,linetype=2)+
      theme_bw() + labs(x =xlab, y=ylab,title = title,color=as.name(colnames(Expdat)[col.ind]))+ 
      scale_color_manual(values=collist)+
      theme(plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
            legend.title = element_text(size=14, face="bold"),
            legend.text = element_text(size=12),
            legend.position=c(1.05, 0.75),
            axis.title = element_text(size=12,face="bold"),
            axis.text = element_text(size=12,face="bold",color = "black"),
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"),
            #axis.text=element_text(size=12,vjust=1),
            plot.margin = margin(0.25, 5, 0.25, 0.25, "cm")
            #axis.ticks.x = element_blank()
      )
  
    }
return(plot)
}

