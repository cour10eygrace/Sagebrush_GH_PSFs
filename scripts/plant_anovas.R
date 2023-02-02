rm(list=ls(all=TRUE))
library(multcomp)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(lme4)
library(vegan)
library(dplyr)
library(tidyr)

#load scaled data from model_setup.R script
load("Data/data_anovas.RData")
load(file= "Data/GH_lmer.RData")

#OBJECTIVE 1----
#microbe anovas
#all microbe responses normal after transformation
microbe$elev.treat<-as.factor(microbe$elev.treat)
microbe$Soilelev<-as.factor(microbe$Soilelev)

NAGmod<-aov(NAG~Soiltreat*Soilelev, microbe)
#NAGmod<-aov(NAG~elev.treat, microbe)
summary(NAGmod)
TukeyHSD(NAGmod)

CBHmod<-aov(CBH~Soiltreat*Soilelev, microbe)
#CBHmod<-aov(CBH~elev.treat, microbe)
summary(CBHmod)
TukeyHSD(CBHmod)

chaomod<-aov(chao~Soiltreat*as.factor(Soilelev), microbe)
#chaomod<-aov(chao~elev.treat, microbe)
summary(chaomod) 
TukeyHSD(chaomod)

mpmod<-aov(mprich~Soiltreat*as.factor(Soilelev), microbe)
#mpmod<-aov(mprich~elev.treat, microbe)
summary(mpmod) 
TukeyHSD(mpmod)

mutmod<-aov(mutrich~Soiltreat*as.factor(Soilelev), microbe)
#mpmod<-aov(mprich~elev.treat, microbe)
summary(mutmod) 
TukeyHSD(mutmod)

pathmod<-aov(pathrich~Soiltreat*Soilelev, microbe)
#mpmod<-aov(mprich~elev.treat, microbe)
summary(pathmod) 
TukeyHSD(pathmod)

sapmod<-aov(saprich~Soiltreat*Soilelev, microbe)
#sapmod<-aov(saprich~elev.treat, microbe)
summary(sapmod) 
TukeyHSD(sapmod)

#OBJECTIVE 2----
#plant models-ANOVAS
PSFmic$Seed.elev<-as.factor(PSFmic$Seed.elev)
PSFmic$elev.treat<-as.factor(PSFmic$elev.treat)
PSFmic$soil.elev<-as.factor(PSFmic$soil.elev)

psfmod<-aov(PSF_u~soil.treat*soil.elev*Seed.elev, PSFmic)
summary(psfmod)
TukeyHSD(psfmod)

heightmod<-aov(height~soil.treat*soil.elev*Seed.elev, PSFmic)
summary(heightmod)
TukeyHSD(heightmod)

rootmod<-aov(root_shootz~soil.treat*soil.elev*Seed.elev, PSFmic)
summary(rootmod)
TukeyHSD(rootmod)

SLAmod<-aov(SLA~soil.treat*soil.elev*Seed.elev, PSFmic)
summary(SLAmod)
TukeyHSD(SLAmod)

LDMCmod<-aov(LDMC~soil.treat*soil.elev*Seed.elev, PSFmic)
summary(LDMCmod)
TukeyHSD(LDMCmod)

CNmod<-aov(leafCN~soil.treat*soil.elev*Seed.elev, PSFmic)
summary(CNmod)
TukeyHSD(CNmod)

Cmod<-aov(leafC~soil.treat*soil.elev*Seed.elev, PSFmic)
summary(Cmod)
TukeyHSD(Cmod)

Nmod<-aov(leafN~soil.treat*soil.elev*Seed.elev, PSFmic)
summary(Nmod)
TukeyHSD(Nmod)

N15mod<-aov(leaf15N~soil.treat*soil.elev*Seed.elev, PSFmic)
summary(N15mod)
TukeyHSD(N15mod)

C13mod<-aov(leaf13C~soil.treat*soil.elev*Seed.elev, PSFmic)
summary(C13mod)
TukeyHSD(C13mod)

#OBJECTIVE 2b----
#mixed effects models
PSFmod<-lmer(PSF_u~ chao + saprich + pathrich+ pcoa2+ (1|soil.elev:soil.treat), data = PSFmic,  
             control = lmerControl(optimizer= "optimx", optCtrl  = list(method="nlminb")))
summary(PSFmod)

heightmod<-lmer(height~ chao + saprich + pathrich+ pcoa2+ (1|soil.elev:soil.treat), data = PSFmic)
summary(heightmod)

rsmod<-lmer(root_shootz~ chao + saprich + pathrich+ pcoa2+ (1|soil.elev:soil.treat), data = PSFmic)
summary(rsmod)

C13mod<-lmer(leaf13C~ chao + saprich + pathrich+ pcoa2+ (1|soil.elev:soil.treat), data = PSFmic)
summary(C13mod)

N15mod<-lmer(leaf15N~ chao + saprich + pathrich+ pcoa2+ (1|soil.elev:soil.treat), data = PSFmic)
summary(N15mod)
coef(N15mod)

#PLOTS----

#microbe plots----
#community comp
microbeplot<-dplyr::select(microbe, NAG, CBH, chao,  mutrich, pathrich, saprich, elev.treat, Soiltreat, 
           Soilelev)%>%#,Dim1,Dim2)
pivot_longer(cols=c(NAG, CBH, chao,  mutrich, pathrich, saprich), #Dim1, Dim2),
                       names_to = 'micpar')
microbeplot$micpar = factor(microbeplot$micpar,
                            levels = c( "chao",  "pathrich", "saprich", "mutrich", "CBH",  "NAG")) #, "Dim1", "Dim2" ))

#Fig 2
ggplot(microbeplot,
       aes(x=as.factor(Soilelev),y=value, fill=Soiltreat)) + 
  geom_boxplot() + facet_wrap(~micpar, scale="free_y")+
  scale_fill_brewer(palette = "Paired", direction = -1)+
  theme_classic() 

#plant growth plots----
plantplot<-dplyr::select(PSFmic,PSF_u, height,root_shootz, SLA, LDMC, leafCN,leafC, leafN, leaf13C, leaf15N, 
                         Seed.elev, soil.elev, soil.treat) %>%
  pivot_longer(cols=c(PSF_u, height,root_shootz, SLA, LDMC, leafCN,leafC, leafN, leaf13C, leaf15N), 
               names_to = 'plantpar')
plantplot$plantpar = factor(plantplot$plantpar,
                            levels = c("PSF_u", "height", "root_shootz", "SLA", "LDMC", "leafC", "leafN", 'leafCN', 'leaf13C', 'leaf15N'))

#FIGURE 5
ggplot(subset(plantplot,plantpar!="SLA"&plantpar!="LDMC"&plantpar!="leafCN"&plantpar!="leafC"&plantpar!="leafN"),
       aes(x=as.factor(Seed.elev),y=value, fill=elev.treat)) + 
  geom_boxplot() + facet_wrap(~plantpar, scale="free_y")+
  scale_fill_brewer(palette = "Paired", direction = -1)+
  theme_classic() 

#Figure S1
ggplot(plantplot,
       aes(x=as.factor(Seed.elev),y=value, fill=elev.treat)) + 
  geom_boxplot() + facet_wrap(~plantpar, scale="free_y")+
  scale_fill_brewer(palette = "Paired", direction = -1)+
  theme_classic() 
