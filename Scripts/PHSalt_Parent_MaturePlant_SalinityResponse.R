library(ggplot2)
library(tidyverse)
library(lmerTest)
library(lme4)
library(car)

#Unit for Fresh weight = gm
#Unit for WP = 1 atm
#atm to Mpa = -WP*0.101325 
#Unit for OA (Osmolarity) = mmol/kg
#mmol/kg to Mpa = -OA (mmol/kg) x2.58X0.001 Van't Hoff equation
#Unit for Na or K = umol/gm of dry leaf

dat<-read.csv("Data/PH_PARENT_SALT_FINAL.csv")


############ SW

modAGFB<-lmer(AGFB ~ Ecotype*Treatment + (1|Unit:Treatment) + (1|Unit), data = dat,REML = F)
summary(modAGFB)
anova(modAGFB,type="I")

# Comment: both G and E 


############## RW
modBGFB<-lmer(BGFB ~ Ecotype*Treatment + (1|Unit:Treatment) + (1|Unit), data = dat,REML = F)
summary(modBGFB)
anova(modBGFB,type="I")
#Anova(modBGFB,type="II",test.statistic="Chisq")

BGFBHAL<-dat %>% dplyr::select(Ecotype, Treatment, BGFB) %>% filter(Ecotype=="HAL")
modBGFBHAL<-aov(BGFB ~ Treatment, data = BGFBHAL)
summary(modBGFBHAL)

BGFBFIL<-dat %>% dplyr::select(Ecotype, Treatment, BGFB) %>% filter(Ecotype=="FIL")
modBGFBFIL<-aov(BGFB ~ Treatment, data = BGFBFIL)
summary(modBGFBFIL)

dat %>% dplyr::select(Ecotype, Treatment, BGFB) %>% filter(Ecotype=="HAL") %>% 
  group_by(Treatment) %>% summarise(BGFB= mean(BGFB,na.rm = T))

# Comment: G and GXE effect

############## RtoS
modRFB<-lmer(RFB ~ Ecotype*Treatment + (1|Unit:Treatment) + (1|Unit), data = dat,REML = F)
summary(modRFB)
anova(modRFB,type="I")
#Anova(modRFB,type="II",test.statistic="Chisq")


# Comment: G and GXE effect

############## PsiW
modPsiW<-lmer(PsiW ~ Ecotype*Treatment + (1|Unit:Treatment) + (1|Unit), data = dat,REML = F)

#Model failed to converge: degenerate  Hessian with 1 negative eigenvalues
summary(modPsiW)
anova(modPsiW,type="I")

PsiWSalinity<-dat %>% dplyr::select(Ecotype, Treatment, PsiW) %>% filter(Treatment=="Salinity")
modPsiWSal<-aov(PsiW ~ Ecotype, data = PsiWSalinity)
summary(modPsiWSal)

dat %>% dplyr::select(Ecotype, Treatment, PsiW) %>% filter(Treatment=="Salinity") %>% 
  group_by(Ecotype) %>% summarise(PsiW= mean(PsiW,na.rm = T))

# Comment: E and GXE effect

############## RWC
modRWC<-lmer(RWC ~ Ecotype*Treatment + (1|Unit:Treatment) + (1|Unit), data = dat,REML = F)

summary(modRWC)
anova(modRWC,type="I")

# Comment: E effect

###### Will keep mmol/kg osmolarity unit

modLO<-lmer(LO ~ Ecotype*Treatment + (1|Unit:Treatment) + (1|Unit), data = dat,REML = F)

summary(modLO)
anova(modLO,type="I")

####### PsiS
# modPsiS<-lmer(PsiS ~ Ecotype*Treatment + (1|Unit:Treatment) + (1|Unit), data = dat,REML = T)
# 
# summary(modPsiS)
# anova(modPsiS)

# Comment: E effect

####### LNa
modLNa<-lmer(LNa ~ Ecotype*Treatment + (1|Unit:Treatment) + (1|Unit), data = dat,REML = F)

summary(modLNa)
anova(modLNa,type="I")

# Comment: E and GXE effect

LNaSalinity<-dat %>% dplyr::select(Ecotype, Treatment, LNa) %>% filter(Treatment=="Salinity")
modLNaSal<-aov(LNa ~ Ecotype, data = LNaSalinity)
summary(modLNaSal)

dat %>% dplyr::select(Ecotype, Treatment, LNa) %>% filter(Treatment=="Salinity") %>% 
  group_by(Ecotype) %>% summarise(LNa= mean(LNa,na.rm = T))

####### LK
modLK<-lmer(LK ~ Ecotype*Treatment + (1|Unit:Treatment) + (1|Unit), data = dat,REML = F)
#Model failed to converge: degenerate  Hessian with 1 negative eigenvalues
summary(modLK)
anova(modLK,type="I")

# Comment: G and E effect

####### NabyK
modNabyK<-lmer(NabyK ~ Ecotype*Treatment + (1|Unit:Treatment) + (1|Unit), data = dat,REML = F)


summary(modNabyK)
anova(modNabyK,type="I")

# Comment: E and GXE effect

NabyKSalinity<-dat %>% dplyr::select(Ecotype, Treatment, NabyK) %>% filter(Treatment=="Salinity")
modNabyKSal<-aov(NabyK ~ Ecotype, data = NabyKSalinity)
summary(modNabyKSal)

dat %>% dplyr::select(Ecotype, Treatment, NabyK) %>% filter(Treatment=="Salinity") %>% 
  group_by(Ecotype) %>% summarise(NabyK= mean(NabyK,na.rm = T))


######## REACTION NORM PLOT

React_theme<-theme_classic(base_size = 10)+
  theme(plot.title = element_text(hjust = 0.5,size=10, face = "bold"),
        axis.title = element_text(size=10,face="bold"),
        axis.text=element_text(size=10,vjust=1,colour = "black",face = "bold"))

####### WATER STATUS RELATED TRAIT
## WP
#y1<-expression(shi)

WP<-dat %>% dplyr::select(Ecotype, Treatment, PsiW) %>% 
  group_by(Ecotype,Treatment) %>% summarize(mean = mean(PsiW,na.rm=T), std_err= sd(PsiW, na.rm=T) / sqrt(n()) )%>% 
  ungroup() %>% 
  mutate(lower = mean - std_err, upper = mean + std_err) %>%  
  ggplot(aes(x = Treatment, y = mean, color = Ecotype,group=Ecotype, ymin = lower, ymax = upper)) + 
  geom_line(size=0.75) + geom_errorbar(size= .5,width=.1,linetype=1)+geom_point()+
  scale_x_discrete(expand = c(0.05, 0.05),labels=c("Control","Salinity")) +
  scale_y_continuous(limits = c(-2.25,-0.75))+
  labs(x="",y = "Leaf water potential\n(Mpa)" ,title="",color="Ecotype")+
  #ylab(expression(psi["W"]*"  MPa"))+
  guides(colour=F)+
  geom_text(aes(x=1,y=-0.75,label="A)"),color="black",size=4)+
  scale_color_manual(values = c("#0072B2","#D55E00"))+
  React_theme

LO<-dat %>% dplyr::select(Ecotype, Treatment, LO) %>% 
  group_by(Ecotype,Treatment) %>% summarize(mean = mean(LO,na.rm=T), std_err= sd(LO, na.rm=T) / sqrt(n()) )%>% 
  ungroup() %>% 
  mutate(lower = mean - std_err, upper = mean + std_err) %>%  
  ggplot(aes(x = Treatment, y = mean, color = Ecotype,group=Ecotype, ymin = lower, ymax = upper)) + 
  geom_line(size=0.75) + geom_errorbar(size= .5,width=.1,linetype=1)+geom_point()+
  scale_x_discrete(expand = c(0.05, 0.05),labels=c("Control","Salinity")) +
  scale_y_continuous(limits = c(350,850))+
  labs(x="",y = "Osmotic Potential\n(mmol/kg)" ,title="",color="Ecotype")+
  guides(colour=F)+
  geom_text(aes(x=1,y=850,label="B)"),color="black",size=4)+
  scale_color_manual(values = c("#0072B2","#D55E00"))+
  React_theme
  
LRWC<-dat %>% dplyr::select(Ecotype, Treatment, RWC) %>% 
  group_by(Ecotype,Treatment) %>% summarize(mean = mean(RWC,na.rm=T), std_err= sd(RWC, na.rm=T) / sqrt(n()) )%>% 
  ungroup() %>% 
  mutate(lower = mean - std_err, upper = mean + std_err) %>%  
  ggplot(aes(x = Treatment, y = mean, color = Ecotype,group=Ecotype, ymin = lower, ymax = upper)) + 
  geom_line(size=0.75) + geom_errorbar(size= .5,width=.1,linetype=1)+geom_point()+
  scale_x_discrete(expand = c(0.05, 0.05),labels=c("Control","Salinity")) +
  scale_y_continuous(limits = c(80,100))+
  labs(x="",y = "Relative water content\n(%)" ,title="",color="      Genotype")+
  #guides(colour=F)+
  geom_text(aes(x=1,y=100,label="C)"),color="black",size=4)+
  scale_color_manual(values = c("#0072B2","#D55E00"),labels=c("FIL2 (Coastal)","HAL2 (Inland)"))+
  React_theme+theme(legend.position = c(.73, .9),
                    legend.title = element_text(size=12,face = "bold"), 
                    legend.text = element_text(size=10,face="bold") ,
                    legend.background = element_rect(fill = "transparent"))

####### GROWTH RELATED TRAITS

####AGFB
AGFB<-dat %>% dplyr::select(Ecotype, Treatment, AGFB) %>% 
  group_by(Ecotype,Treatment) %>% summarize(mean = mean(AGFB,na.rm=T), std_err= sd(AGFB, na.rm=T) / sqrt(n()) )%>% 
  ungroup() %>% 
  mutate(lower = mean - std_err, upper = mean + std_err) %>%  
  ggplot(aes(x = Treatment, y = mean, color = Ecotype,group=Ecotype, ymin = lower, ymax = upper)) + 
  geom_line(size=0.75) + geom_errorbar(size= .5,width=.1,linetype=1)+geom_point()+
  scale_x_discrete(expand = c(0.05, 0.05),labels=c("Control","Salinity")) +
  scale_y_continuous(limits = c(0.9,3.5))+
  labs(x="",y = "Above ground fresh biomass\n(g)" ,title="",color="Ecotype")+
  guides(colour=F)+
  geom_text(aes(x=1,y=3.5,label="D)"),color="black",size=4)+
  scale_color_manual(values = c("#0072B2","#D55E00"))+
  React_theme

####BGFB
BGFB<-dat %>% dplyr::select(Ecotype, Treatment, BGFB) %>% 
  group_by(Ecotype,Treatment) %>% summarize(mean = mean(BGFB,na.rm=T), std_err= sd(BGFB, na.rm=T) / sqrt(n()) )%>% 
  ungroup() %>% 
  mutate(lower = mean - std_err, upper = mean + std_err) %>%  
  ggplot(aes(x = Treatment, y = mean, color = Ecotype,group=Ecotype, ymin = lower, ymax = upper)) + 
  geom_line(size=0.75) + geom_errorbar(size= .5,width=.1,linetype=1)+geom_point()+
  scale_x_discrete(expand = c(0.05, 0.05),labels=c("Control","Salinity")) +
  scale_y_continuous(limits = c(0.25,0.7))+
  labs(x="",y = "Below ground fresh biomass\n(g)" ,title="",color="Ecotype")+
  guides(colour=F)+
  geom_text(aes(x=1,y=0.7,label="E)"),color="black",size=4)+
  scale_color_manual(values = c("#0072B2","#D55E00"))+
  React_theme

RFB<-dat %>% dplyr::select(Ecotype, Treatment, RFB) %>% 
  group_by(Ecotype,Treatment) %>% summarize(mean = mean(RFB,na.rm=T), std_err= sd(RFB, na.rm=T) / sqrt(n()) )%>% 
  ungroup() %>% 
  mutate(lower = mean - std_err, upper = mean + std_err) %>%  
  ggplot(aes(x = Treatment, y = mean, color = Ecotype,group=Ecotype, ymin = lower, ymax = upper)) + 
  geom_line(size=0.75) + geom_errorbar(size= .5,width=.1,linetype=1)+geom_point()+
  scale_x_discrete(expand = c(0.05, 0.05),labels=c("Control","Salinity")) +
  scale_y_continuous(limits = c(2.8,6))+
  labs(x="",y = "Ratio of biomass\n(above ground/below ground)" ,title="",color="Ecotype")+
  guides(colour=F)+
  geom_text(aes(x=1,y=6,label="F)"),color="black",size=4)+
  scale_color_manual(values = c("#0072B2","#D55E00"))+
  React_theme

###### Leaf Ion
LNa<-dat %>% dplyr::select(Ecotype, Treatment, LNa) %>% 
  group_by(Ecotype,Treatment) %>% summarize(mean = mean(LNa,na.rm=T), std_err= sd(LNa, na.rm=T) / sqrt(n()) )%>% 
  ungroup() %>% 
  mutate(lower = mean - std_err, upper = mean + std_err) %>%  
  ggplot(aes(x = Treatment, y = mean, color = Ecotype,group=Ecotype, ymin = lower, ymax = upper)) + 
  geom_line(size=0.75) + geom_errorbar(size= .5,width=.1,linetype=1)+geom_point()+
  scale_x_discrete(expand = c(0.05, 0.05),labels=c("Control","Salinity")) +
  scale_y_continuous(limits = c(0,1250))+
  labs(x="",y = "Na\n(umol/g of dry leaf)" ,title="",color="Ecotype")+
  guides(colour=F)+
  geom_text(aes(x=1,y=1250,label="G)"),color="black",size=4)+
  scale_color_manual(values = c("#0072B2","#D55E00"))+
  React_theme


LK<-dat %>% dplyr::select(Ecotype, Treatment, LK) %>% 
  group_by(Ecotype,Treatment) %>% summarize(mean = mean(LK,na.rm=T), std_err= sd(LK, na.rm=T) / sqrt(n()) )%>% 
  ungroup() %>% 
  mutate(lower = mean - std_err, upper = mean + std_err) %>%  
  ggplot(aes(x = Treatment, y = mean, color = Ecotype,group=Ecotype, ymin = lower, ymax = upper)) + 
  geom_line(size=0.75) + geom_errorbar(size= .5,width=.1,linetype=1)+geom_point()+
  scale_x_discrete(expand = c(0.05, 0.05),labels=c("Control","Salinity")) +
  scale_y_continuous(limits = c(520,750))+
  labs(x="",y = "K\n(umol/g of dry leaf)" ,title="",color="Ecotype")+
  guides(colour=F)+
  geom_text(aes(x=1,y=750,label="H)"),color="black",size=4)+
  scale_color_manual(values = c("#0072B2","#D55E00"))+
  React_theme

LNaK<-dat %>% dplyr::select(Ecotype, Treatment, NabyK) %>% 
  group_by(Ecotype,Treatment) %>% summarize(mean = mean(NabyK,na.rm=T), std_err= sd(NabyK, na.rm=T) / sqrt(n()) )%>% 
  ungroup() %>% 
  mutate(lower = mean - std_err, upper = mean + std_err) %>%  
  ggplot(aes(x = Treatment, y = mean, color = Ecotype,group=Ecotype, ymin = lower, ymax = upper)) + 
  geom_line(size=0.75) + geom_errorbar(size= .5,width=.1,linetype=1)+geom_point()+
  scale_x_discrete(expand = c(0.05, 0.05),labels=c("Control","Salinity")) +
  scale_y_continuous(limits = c(0,2.3))+
  labs(x="",y = "Na/K\n(ratio)" ,title="",color="Ecotype")+
  guides(colour=F)+
  geom_text(aes(x=1,y=2.3,label="I)"),color="black",size=4)+
  scale_color_manual(values = c("#0072B2","#D55E00"))+
  React_theme


library(cowplot)

plot.with.inset <-
  ggdraw() +
  draw_plot(WP,x = 0.0, y= 0.66,width = 0.32,height = 0.35)+
  draw_plot(LO,x = 0.33, y= 0.66,width = 0.32,height = 0.35)+
  draw_plot(LRWC,x = 0.66, y= 0.66,width = 0.32,height = 0.35)+
  draw_plot(AGFB,x = 0.0, y= 0.33,width = 0.32,height = 0.35)+
  draw_plot(BGFB,x = 0.33, y= 0.33,width = 0.32,height = 0.35)+
  draw_plot(RFB,x = 0.66, y= 0.33,width = 0.32,height = 0.35)+
  draw_plot(LNa,x = 0.0, y= 0.0,width = 0.32,height = 0.35)+
  draw_plot(LK,x = 0.33, y= 0.0,width = 0.32,height = 0.35)+
  draw_plot(LNaK,x = 0.66, y= 0.0,width = 0.32,height = 0.35)
  

tiff("Plots/Figure_MaturePlant_SalinityResponse.tiff",width=7.5,height=9,units="in",res=300)
plot.with.inset
dev.off()

