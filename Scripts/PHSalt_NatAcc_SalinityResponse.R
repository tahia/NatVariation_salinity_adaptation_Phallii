setwd("/home/taslima/data/JuengerLab/Research_Article_Preps/NatVariation_salinity_adaptation_Phallii")
library(tidyverse)
library(gtools)
library(ggplot2)
library("xlsx")
library(tidyverse)
library(lme4)
library(raster)
library(maps)
library(rgdal)

dat<-read.csv("Data/PHNatAcc_SALT_EXP_Phenotypes_Final.csv")


library(lmerTest)
model1<-lmer(RWC ~ Treat*POPCLASS + (1 | Treat:Tray)  + (1 | POPCLASS:POPID), data = dat,REML = F)
summary(model1)
anova(model1,type="I")

model2<-lmer(BGB ~ Treat*POPCLASS + (1 | Treat:Tray) + (1 | POPCLASS:POPID), data = dat,REML = F)
summary(model2)
anova(model2,type="I")

model3<-lmer(AGB ~ Treat*POPCLASS + (1 | Treat:Tray)  + (1 | POPCLASS:POPID), data = dat,REML = F)
summary(model3)
anova(model3,type="I")

as_tibble(dat) %>%  dplyr::select(POPID, POPCLASS,POPGROUP,Treat,AGB) %>% 
  group_by(Treat, POPCLASS) %>% drop_na(AGB) %>% summarise(MeanAGB= mean(AGB,na.rm = T), 
                                                                         SD=sd(AGB,na.rm=T), 
                                                                         n = n(), 
                                                                         SE = SD/sqrt(n))

as_tibble(dat) %>%  dplyr::select(POPID, POPCLASS,POPGROUP,Treat,BGB) %>% 
  group_by(Treat, POPCLASS) %>% drop_na(BGB) %>% summarise(MeanBGB= mean(BGB,na.rm = T), 
                                                          SD=sd(BGB,na.rm=T), 
                                                          n = n(), 
                                                          SE = SD/sqrt(n))


######### REACTION NORM PLOT
#########LRWC
RWC_Genotypewise<-as_tibble(dat) %>%  dplyr::select(POPID, POPCLASS,POPGROUP,Treat,RWC) %>% 
  group_by(Treat, POPID,POPCLASS) %>% drop_na(RWC) %>% summarise(MeanRWC= mean(RWC,na.rm = T), 
                                                                  SD=sd(RWC,na.rm=T), 
                                                                  n = n(), 
                                                                  SE = SD/sqrt(n)) %>% 
  mutate(lower = MeanRWC - SE, upper = MeanRWC +SE)

RWC<-dat %>%   dplyr::select(Treat, POPGROUP,POPCLASS, RWC)  %>% 
  group_by(POPCLASS,Treat) %>% drop_na(RWC) %>% summarize(mean = mean(RWC,na.rm=T), 
                                                           std_err= sd(RWC, na.rm=T) / sqrt(n()) )%>% 
  ungroup() %>% 
  mutate(lower = mean - std_err, upper = mean + std_err) %>% 
  ggplot(aes(x = Treat, y = mean, color = POPCLASS,group=POPCLASS, ymin = lower, ymax = upper)) + 
  geom_line(size=0.75) + geom_errorbar(size= .5,width=.1,linetype=1)+geom_point()+
  geom_line(data = RWC_Genotypewise,
            aes(x=Treat,y=MeanRWC,color=POPCLASS,group=POPID, ymin = lower, ymax = upper),
            size=0.3,linetype=2)+
  geom_text(aes(x=1,y=100,label="C)"),color="black",size=5)+
  scale_x_discrete(expand = c(0.05, 0.05),labels=c("Control","Salinity")) +
  scale_y_continuous(limits = c(90,100),breaks = c(90,95,100)) +
  labs(x="",y = "Relative water content\n (%)",title="",color="Ecotype")+
  #guides(colour=F)+
  scale_color_manual(values = c("#0072B2","#D55E00"))+
  theme_classic(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
        axis.title = element_text(size=14,face="bold"),
        axis.text=element_text(size=14,vjust=1,colour = "black"),
        legend.position=c(0.5,0.9))

#lg2 <- get_legend(LRWC)

#LRWC<-LRWC+guides(colour=F)

#########AGB
AGB_Genotypewise<-as_tibble(dat) %>%  dplyr::select(POPID, POPCLASS,POPGROUP,Treat,AGB) %>% 
  group_by(Treat, POPID,POPCLASS) %>% drop_na(AGB) %>% summarise(MeanAGB= mean(AGB,na.rm = T), 
                                                                SD=sd(AGB,na.rm=T), 
                                                                n = n(), 
                                                                SE = SD/sqrt(n)) %>% 
  mutate(lower = MeanAGB - SE, upper = MeanAGB +SE)

AGB_Genotypewise_STI<-AGB_Genotypewise %>% dplyr::select(Treat,POPID,POPCLASS,MeanAGB) %>% 
  spread(Treat,MeanAGB) %>% 
  mutate(STI_AGB=Stress/Control) %>% group_by(POPCLASS)

AGB_Genotypewise_STI %>% dplyr::select(POPCLASS,STI_AGB) %>%  group_by(POPCLASS) %>% 
  summarise(Mean_STI_AGB= mean(STI_AGB,na.rm = T))

modelSWSTI<-lm(STI_AGB ~ POPCLASS, data = AGB_Genotypewise_STI)
summary(modelSWSTI)

AGB_Genotypewise %>% dplyr::select(Treat,POPID,POPCLASS,MeanAGB) %>% spread(Treat,MeanAGB) %>% 
  mutate(STI_AGB=Stress/Control) %>% group_by(POPCLASS) %>% summarise(MeanSTI_AGB=mean(STI_AGB))

AGB<-dat %>%   dplyr::select(Treat, POPGROUP,POPCLASS, AGB)  %>% 
  group_by(POPCLASS,Treat) %>% drop_na(AGB) %>% summarize(mean = mean(AGB,na.rm=T), std_err= sd(AGB, na.rm=T) / sqrt(n()) )%>% 
  ungroup() %>% 
  mutate(lower = mean - std_err, upper = mean + std_err) %>% 
  ggplot(aes(x = Treat, y = mean, color = POPCLASS,group=POPCLASS, ymin = lower, ymax = upper)) + 
  geom_line(size=0.75) + geom_errorbar(size= .5,width=.1,linetype=1)+geom_point()+
  geom_line(data = AGB_Genotypewise,
            aes(x=Treat,y=MeanAGB,color=POPCLASS,group=POPID, ymin = lower, ymax = upper),
            size=0.3,linetype=2)+
  scale_x_discrete(expand = c(0.05, 0.05),labels=c("Control","Salinity")) +
  scale_y_continuous(limits = c(0.04,0.16),breaks = c(0.05,0.1,0.15))+
  geom_text(aes(x=1,y=0.16,label="D)"),color="black",size=5)+
  labs(x="",y = "Above ground biomass\n (g)",title="",color="Ecotype")+
  guides(colour=F)+
  scale_color_manual(values = c("#0072B2","#D55E00"))+
  theme_classic(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
        axis.title = element_text(size=14,face="bold"),
        axis.text=element_text(size=14,vjust=1,colour = "black"))

######### BGB

BGB_Genotypewise<-as_tibble(dat) %>%  dplyr::select(POPID, POPCLASS,POPGROUP,Treat,BGB) %>% 
  group_by(Treat, POPID,POPCLASS) %>% drop_na(BGB) %>% summarise(MeanBGB= mean(BGB,na.rm = T), 
                                                          SD=sd(BGB,na.rm=T), 
                                                          n = n(), 
                                                          SE = SD/sqrt(n)) %>% 
  mutate(lower = MeanBGB - SE, upper = MeanBGB +SE)
  
BGB_Genotypewise_STI<-BGB_Genotypewise %>% dplyr::select(Treat,POPID,POPCLASS,MeanBGB) %>% 
  spread(Treat,MeanBGB) %>% 
  mutate(STI_BGB=Stress/Control) %>% group_by(POPCLASS)

BGB_Genotypewise_STI %>% dplyr::select(POPCLASS,STI_BGB) %>%  group_by(POPCLASS) %>% 
  summarise(Mean_STI_BGB= mean(STI_BGB,na.rm = T))


modelBGBSTI<-lm(STI_BGB ~ POPCLASS, data = BGB_Genotypewise_STI)
summary(modelBGBSTI)


BGB_Genotypewise %>% dplyr::select(Treat,POPID,POPCLASS,MeanBGB) %>% spread(Treat,MeanBGB) %>% 
  mutate(STI_BGB=Stress/Control) %>% group_by(POPCLASS) %>% summarise(MeanSTI_BGB=mean(STI_BGB))

BGB<-dat %>% dplyr::select(Treat, POPGROUP,POPCLASS, BGB)  %>% 
  group_by(POPCLASS,Treat) %>% drop_na(BGB) %>% summarize(mean = mean(BGB,na.rm=T), std_err= sd(BGB, na.rm=T) / sqrt(n()) )%>% 
  ungroup() %>% 
  mutate(lower = mean - std_err, upper = mean + std_err) %>% 
  ggplot(aes(x = Treat, y = mean, color = POPCLASS,group=POPCLASS, ymin = lower, ymax = upper)) + 
  geom_line(size=0.75) + geom_errorbar(size= .5,width=.1,linetype=1)+geom_point()+
  geom_line(data = BGB_Genotypewise,
            aes(x=Treat,y=MeanBGB,color=POPCLASS,group=POPID, ymin = lower, ymax = upper),
            size=0.3,linetype=2)+
  # geom_errorbar(data = RW_Genotypewise,
  #               aes(x=Treat,y=MeanRW,color=POPCLASS,group=POPID, ymin = lower, ymax = upper),
  #               size= .15,width=.05,linetype=1)+
  #geom_point()+
  scale_x_discrete(expand = c(0.05, 0.05),labels=c("Control","Salinity")) +
  scale_y_continuous(limits = c(0.01,0.04),breaks = seq(0.01,0.04,0.01))+
  geom_text(aes(x=1,y=0.04,label="E)"),color="black",size=5)+
  labs(x="",y = "Below ground biomass\n (g)",title="",color="Ecotype")+
  guides(colour=F)+
  scale_color_manual(values = c("#0072B2","#D55E00"))+
  theme_classic(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
        axis.title = element_text(size=14,face="bold"),
        axis.text=element_text(size=14,vjust=1,colour = "black"))


################# RASTER PLOT
Struc<-read.csv("Data/PH_Pop_GeoLocation_Subpop.csv")
ext<-extent(-115, -95, 24, 38)
Struc$Size<-NA
Struc$Shape<-NA
Struc$Size<-factor(Struc$Size)

Na2018<-raster("Data/Na_conc_2018/Na_conc_2018/Na_conc_2018.tif")
#Na2018<-raster("Data/na_tw-2018/na_tw-2018.tif")
crs(Na2018)
Na.gcs <- projectRaster(Na2018, crs="+proj=longlat +datum=WGS84")
PHSalt<-crop(Na.gcs,ext)

dat<-as.data.frame(PHSalt,xy=T)
colnames(dat)[3]<-"Na_conc_2018"
dat<-dat[-which(is.na(dat$Na_conc_2018)),]

mapdata <- read.csv("Data/state-medal-count.csv", header=TRUE, stringsAsFactors=FALSE)
states <- map_data("state")
substates<-states[which(states$region %in% c("texas","new mexico","arizona","oklahoma")),]

p<-ggplot() +
  geom_point(data =  dat, aes(x = x, y = y, color = Na_conc_2018),alpha=0.5) +
  scale_color_gradientn(colours =  terrain.colors(12))+
  #scale_colour_gradient(low = "green", high = "yellow", na.value = NA)+
  geom_polygon(data=substates,aes(x = long, y = lat,group = group),alpha=0, color = "white") + 
  coord_quickmap()+
  labs(x="",y="",color="Sodium deposition \nconcentration (ug/L)",fill="Ecotype")+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10,face="bold"))



q<-p+ geom_point(data = Struc,aes(x=Longitude, y=Latitude,fill=VAR,size=Size),
                 shape=21,alpha=1,colour="black")+
  scale_fill_manual(values = c("#0072B2","#D55E00","#CC79A7", "#0072B2","#000000"),
                    labels=c("FIL","HAL"))+
  scale_size_manual(values = c(2,4),labels=c("1","2"))+
  coord_quickmap()

q<-q+geom_curve(aes(yend = 27.64956, xend = -97.40386, y = 27, x = -97.5), 
                colour = "black", size=1, curvature = -0.2,
                arrow = arrow(length = unit(0.015, "npc"))) + 
  geom_label(aes(x = -97.5, y = 27, label = "FIL2"), 
             hjust = 0, vjust = 0.5, colour = "black", 
             fill = "transparent", label.size = NA,
             family="Helvetica", size = 3,fontface="bold")+
  geom_curve(aes(yend = 30.18502, xend = -97.87398, y = 29.5, x = -96.5), 
             colour = "black", size=1, curvature = -0.2,
             arrow = arrow(length = unit(0.015, "npc"))) + 
  geom_label(aes(x = -96.5, y = 29.5, label = "HAL2"), 
             hjust = 0, vjust = 0.5, colour = "black", 
             fill = "transparent", label.size = NA,
             family="Helvetica", size = 3,fontface="bold")+
  geom_text(aes(x=-116,y=36,label="A)"),color="black",size=5)+
  guides(size=F)

q1<-q+guides(fill=F)
q2<-q+guides(color=F)
r<-q+guides(colour="none",fill="none")

rastL1 <- get_legend(q1)
rastL2 <- get_legend(q2)

library(cowplot)

soil_salinity<-as.data.frame(cbind(Habitat=c(rep("Coastal",3),rep("Inland",2)),
                                   SoilSodium=c(154,254,206,29,32)
                                   ))
soil_salinity$SoilSodium<-as.numeric(soil_salinity$SoilSodium)

soil_salinity_mean<-as_tibble(soil_salinity) %>% group_by(Habitat) %>% 
  summarize(Sodium = mean(SoilSodium,na.rm=T), 
            std_err= sd(SoilSodium, na.rm=T) / sqrt(n()) ) %>% 
  ggplot((aes(x=Habitat,y=Sodium)))+geom_col(aes(fill=Habitat))+
  geom_errorbar(aes(ymin = Sodium-std_err, ymax = Sodium+std_err),
                size= .5,width=.1,linetype=1)+
  labs(x="",y = "Soil salinity (ppm)",title="")+
  guides(fill=F)+
  scale_fill_manual(values = c("#0072B2","#D55E00"))+
  scale_y_continuous(limits = c(0,250),breaks = c(0,100,200))+
  theme_classic(base_size = 12)+
  geom_text(aes(x=0.6,y=240,label="B)"),color="black",size=5)+
  theme(plot.title = element_text(hjust = 0.5,size=8, face = "bold"),
        axis.title = element_text(size=12,face="bold"),
        axis.text=element_text(size=12,vjust=1,colour = "black"),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.background = element_rect(fill = "transparent", colour ="transparent"))



plot.with.inset <-
  ggdraw() +
  draw_plot(r,x = 0.0, y= 0.60,width = 0.5,height = 0.5)+
  draw_plot(rastL1,x=0.38,y=0.6,width= 0.01, height=0.01)+
  draw_plot(rastL2,x=0.45,y=0.6,width= 0.01, height=0.01)+
  draw_plot(soil_salinity_mean,x=0,y=0.5,width= 0.3, height=0.3)+
  draw_plot(RWC,x=0.5,y=0.5,width = 0.5,height = 0.52)+
  draw_plot(AGB,x = 0, y = 0,width = 0.5,height = 0.52) +
  draw_plot(BGB,x = 0.5, y = 0,width = 0.5,height = 0.52) 

tiff("Plots/Figure_PHNatAcc_SalinityResponse.tiff",width=9,height=9,units="in",res=300)
plot.with.inset
dev.off()

