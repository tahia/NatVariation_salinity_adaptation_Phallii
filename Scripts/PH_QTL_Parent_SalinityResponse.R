setwd("/home/taslima/data/JuengerLab/Research_Article_Preps/NatVariation_salinity_adaptation_Phallii")
library(tidyverse)
library(gtools)
library(ggplot2)
library("xlsx")
library(tidyverse)

dat<-read.csv("Data/PH_QTLEXP_Parent_Pheno_FINAL.csv")

library(lmerTest)
model1<-lmer(RWC ~ Treatment*Genotype + (1|Cohort) + (1|Tray)+ (1|Treatment:Tray), data = dat,REML = F)
summary(model1)
anova(model1,type="I")

model2<-lmer(BGB ~ Treatment*Genotype + (1|Cohort) + (1|Tray)+ (1|Treatment:Tray) , data = dat,REML = F)
summary(model2)
anova(model2,type="I")

model3<-lmer(AGB ~ Treatment*Genotype + (1|Cohort) + (1|Tray) + (1|Treatment:Tray) , data = dat,REML = F)
summary(model3)
anova(model3,type = "I")

model4<-lmer(RB ~ Treatment*Genotype + (1|Cohort) + (1|Tray) + (1|Treatment:Tray) , data = dat,REML = F)
summary(model4)
anova(model4,type = "I")

dat2<-dat[which(dat$Cohort==1),]


model5<-lmer(K ~ Treatment*Genotype + (1|Tray), data = dat2)
summary(model5)
anova(model5,type = "I")

model6<-lmer(Na ~ Treatment*Genotype + (1|Tray), data = dat2)
summary(model6)
anova(model6,type = "I")

model7<-lmer(NabyK ~ Treatment*Genotype + (1|Tray), data = dat2)
summary(model7)
anova(model7,type = "I")


dat2 %>%   dplyr::select(Treatment, Genotype,K)  %>% 
  group_by(Treatment, Genotype) %>% drop_na(K) %>% summarize(mean = mean(K,na.rm=T), 
                                                          std_err= sd(K, na.rm=T) / sqrt(n()) )%>% 
  ungroup() %>% 
  mutate(lower = mean - std_err, upper = mean + std_err) %>% 
  ggplot(aes(x = Treatment, y = mean, color = Genotype,group=Genotype, ymin = lower, ymax = upper)) + 
  geom_line(size=0.75) + geom_errorbar(size= .5,width=.1,linetype=1)+geom_point()+
  scale_x_discrete(expand = c(0.05, 0.05),labels=c("Control","Salinity")) +
  #scale_y_continuous(limits = c(90,100),breaks = c(90,95,100)) +
  labs(x="",y = "K",title="",color="Genotype")+
  #guides(colour=F)+
  scale_color_manual(values = c("#0072B2","#D55E00"))+
  theme_classic(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
        axis.title = element_text(size=14,face="bold"),
        axis.text=element_text(size=14,vjust=1,colour = "black"),
        legend.position=c(0.5,0.9))


dat2 %>%   dplyr::select(Treatment, Genotype,Na)  %>% 
  group_by(Treatment, Genotype) %>% drop_na(Na) %>% summarize(mean = mean(Na,na.rm=T), 
                                                             std_err= sd(Na, na.rm=T) / sqrt(n()) )%>% 
  ungroup() %>% 
  mutate(lower = mean - std_err, upper = mean + std_err) %>% 
  ggplot(aes(x = Treatment, y = mean, color = Genotype,group=Genotype, ymin = lower, ymax = upper)) + 
  geom_line(size=0.75) + geom_errorbar(size= .5,width=.1,linetype=1)+geom_point()+
  scale_x_discrete(expand = c(0.05, 0.05),labels=c("Control","Salinity")) +
  #scale_y_continuous(limits = c(90,100),breaks = c(90,95,100)) +
  labs(x="",y = "Na",title="",color="Genotype")+
  #guides(colour=F)+
  scale_color_manual(values = c("#0072B2","#D55E00"))+
  theme_classic(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
        axis.title = element_text(size=14,face="bold"),
        axis.text=element_text(size=14,vjust=1,colour = "black"),
        legend.position=c(0.5,0.9)) 

dat2 %>%   dplyr::select(Treatment, Genotype,NabyK)  %>% 
  group_by(Treatment, Genotype) %>% drop_na(NabyK) %>% summarize(mean = mean(NabyK,na.rm=T), 
                                                                     std_err= sd(NabyK, na.rm=T) / sqrt(n()) )%>% 
  ungroup() %>% 
  mutate(lower = mean - std_err, upper = mean + std_err) %>% 
  ggplot(aes(x = Treatment, y = mean, color = Genotype,group=Genotype, ymin = lower, ymax = upper)) + 
  geom_line(size=0.75) + geom_errorbar(size= .5,width=.1,linetype=1)+geom_point()+
  scale_x_discrete(expand = c(0.05, 0.05),labels=c("Control","Salinity")) +
  #scale_y_continuous(limits = c(90,100),breaks = c(90,95,100)) +
  labs(x="",y = "NabyK",title="",color="Genotype")+
  #guides(colour=F)+
  scale_color_manual(values = c("#0072B2","#D55E00"))+
  theme_classic(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
        axis.title = element_text(size=14,face="bold"),
        axis.text=element_text(size=14,vjust=1,colour = "black"),
        legend.position=c(0.5,0.9)) 

dat %>%   dplyr::select(Treatment, Genotype,AGB)  %>% 
  group_by(Treatment, Genotype) %>% drop_na(AGB) %>% summarize(mean = mean(AGB,na.rm=T), 
                                                             std_err= sd(AGB, na.rm=T) / sqrt(n()) )%>% 
  ungroup() %>% 
  mutate(lower = mean - std_err, upper = mean + std_err) %>% 
  ggplot(aes(x = Treatment, y = mean, color = Genotype,group=Genotype, ymin = lower, ymax = upper)) + 
  geom_line(size=0.75) + geom_errorbar(size= .5,width=.1,linetype=1)+geom_point()+
  scale_x_discrete(expand = c(0.05, 0.05),labels=c("Control","Salinity")) +
  #scale_y_continuous(limits = c(90,100),breaks = c(90,95,100)) +
  labs(x="",y = "AGB",title="",color="Genotype")+
  #guides(colour=F)+
  scale_color_manual(values = c("#0072B2","#D55E00"))+
  theme_classic(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
        axis.title = element_text(size=14,face="bold"),
        axis.text=element_text(size=14,vjust=1,colour = "black"),
        legend.position=c(0.5,0.9))

dat %>%   dplyr::select(Treatment, Genotype,BGB)  %>% 
  group_by(Treatment, Genotype) %>% drop_na(BGB) %>% summarize(mean = mean(BGB,na.rm=T), 
                                                               std_err= sd(BGB, na.rm=T) / sqrt(n()) )%>% 
  ungroup() %>% 
  mutate(lower = mean - std_err, upper = mean + std_err) %>% 
  ggplot(aes(x = Treatment, y = mean, color = Genotype,group=Genotype, ymin = lower, ymax = upper)) + 
  geom_line(size=0.75) + geom_errorbar(size= .5,width=.1,linetype=1)+geom_point()+
  scale_x_discrete(expand = c(0.05, 0.05),labels=c("Control","Salinity")) +
  #scale_y_continuous(limits = c(90,100),breaks = c(90,95,100)) +
  labs(x="",y = "BGB",title="",color="Genotype")+
  #guides(colour=F)+
  scale_color_manual(values = c("#0072B2","#D55E00"))+
  theme_classic(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
        axis.title = element_text(size=14,face="bold"),
        axis.text=element_text(size=14,vjust=1,colour = "black"),
        legend.position=c(0.5,0.9))

dat %>%   dplyr::select(Treatment, Genotype,RWC)  %>% 
  group_by(Treatment, Genotype) %>% drop_na(RWC) %>% summarize(mean = mean(RWC,na.rm=T), 
                                                               std_err= sd(RWC, na.rm=T) / sqrt(n()) )%>% 
  ungroup() %>% 
  mutate(lower = mean - std_err, upper = mean + std_err) %>% 
  ggplot(aes(x = Treatment, y = mean, color = Genotype,group=Genotype, ymin = lower, ymax = upper)) + 
  geom_line(size=0.75) + geom_errorbar(size= .5,width=.1,linetype=1)+geom_point()+
  scale_x_discrete(expand = c(0.05, 0.05),labels=c("Control","Salinity")) +
  #scale_y_continuous(limits = c(90,100),breaks = c(90,95,100)) +
  labs(x="",y = "RWC",title="",color="Genotype")+
  #guides(colour=F)+
  scale_color_manual(values = c("#0072B2","#D55E00"))+
  theme_classic(base_size = 14)+
  theme(plot.title = element_text(hjust = 0.5,size=14, face = "bold"),
        axis.title = element_text(size=14,face="bold"),
        axis.text=element_text(size=14,vjust=1,colour = "black"),
        legend.position=c(0.5,0.9))
