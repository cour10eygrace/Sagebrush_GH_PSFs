rm(list=ls(all=TRUE))
library(lme4)
library(lmerTest)
library(brms)
library(dplyr)
library(optimx)
library(tidyr)
library(ggplot2)

#load and organize data----
GH1<-read.csv("Data/Green House experiment plus leaf traits.csv")
GH1$soil.treat=factor(GH1$soil.treat)#run this to fix soil.treat from 5 levels to 4

#Live vs. sterile inoculum 
summary(lm(GH1$Total.biomass.final..g.~as.factor(Live.inoc), GH1))
hist(GH1$Total.biomass.final..g.)
hist(log(GH1$Total.biomass.final..g.))
summary(lm(log(GH1$Total.biomass.final..g.)~as.factor(Live.inoc), GH1))
t.test(GH1$Total.biomass.final..g.~as.factor(Live.inoc), GH1)

#Alpha diversity 
alphadiv<-read.csv("Data/alphadiv_averaged.csv")
GH1<-unite(GH1, elev.treat, soil.elev, soil.treat, remove=F)
join<-dplyr::select(GH1, elev.treat, Pot)

alphadiv<-left_join(alphadiv, join)

#Funguild 
#from funguild_calculate.R
funguild<-read.csv("Data/funguild_abund_rich_bysample.csv")%>%rename(ID=sample)%>%dplyr::select(-X)
metadataGH<-read.csv("Data/map_SampleSeq_info.csv")%>%rename(ID=X.SampleID)
funguild_meta<-left_join(funguild, metadataGH)

funguild<-left_join(funguild_meta, join)

#community comp
#from betadiv.R script
PCoA_bray_coords<-read.csv("Data/PCoA_bray_coords.csv")%>%rename(ID=X)
PCoA_meta<-left_join(PCoA_bray_coords, metadataGH)

PCoA_bray_coords<-left_join(PCoA_meta, join)
PSF<-subset(GH1, Live.inoc>0, na.rm=T)

PSFmic<-left_join(PSF,  alphadiv, by="Pot")%>%distinct(.)
PSFmic<-left_join(PSFmic, funguild,  by="Pot")%>%distinct(.)
PSFmic<-left_join(PSFmic, PCoA_bray_coords, by="Pot")%>%distinct(.)
PSFmic<-select(PSFmic, -elev.treat.y)%>%rename(elev.treat=elev.treat.x)#fix dup col

PSFmic<-mutate(PSFmic,root_shoot=Final.BG.biomass..g./Final.AG.biomass..g., 
               RGR=Total.biomass.final..g./No_days)

#enzymes-remove outliers
PSFmic<-group_by(filter(PSFmic, NAG.activity.inoc<1000), elev.treat)%>%mutate(NAG_avg=mean(NAG.activity.inoc), 
                                              CBH_avg=mean(CBH.activity.inoc))

#ORDINATE LEAF TRAITS----
traits<-dplyr::select(PSFmic, SLA..cm2.g.,LDMC..g.g.,Wt..C., Wt..N.,
                      C.N.ratio ,X.13C, X.15N , Pot)

##markos code updated
PCA_trt = prcomp(traits[, c(2:8)], scale.=TRUE)
biplot(PCA_trt, choices = c(1,2)) #show all 
biplot(PCA_trt, choices = c(3,4))#show SLA, 15N, 13C better 

##this gives loadings
PCA_trt

##this give amount explained by each axis
summary(PCA_trt) 
PCA_traitsdf<-as.data.frame(PCA_trt$x)
PCA_traitsdf<-PCA_traitsdf[,c(1:4)]

#merge back with other data 
PCA_traitsdf<-cbind(PCA_traitsdf, traits)
PSFmic<-left_join(PSFmic, PCA_traitsdf)

#test for normality 
shapiro.test(PSFmic$PC1)#p=0.2 normal
shapiro.test(PSFmic$PC2)#p=0.0.6 normal
shapiro.test(PSFmic$PC3)#p=0.5 normal
shapiro.test(PSFmic$PC4)#p=0.3 normal
shapiro.test(log(PSFmic$Final.height..cm.))#p=0.3 log normal 
shapiro.test(sqrt(PSFmic$SLA..cm2.g.))#p=0.09 sqrt normal 
shapiro.test(log(PSFmic$LDMC..g.g.))#p=0.6 log normal 
shapiro.test(log(PSFmic$C.N.ratio))#p=0.2 log normal 
shapiro.test(PSFmic$Wt..C.)#p=0.4 normal 
shapiro.test(PSFmic$Wt..N.)#p=0.2 normal 
shapiro.test(PSFmic$X.15N)#p=0.8 normal 
shapiro.test(sqrt(PSFmic$X.13C+36))#p=0.5 sqrt normal 
shapiro.test(PSFmic$PSFu)#p<0.05 
shapiro.test(log(PSFmic$PSFu+1))#p<0.05 
shapiro.test(sqrt(PSFmic$PSFu+1))#p=0.07 

#try transform tukey 
library(rcompanion)
transformTukey(PSFmic$PSFu+1)#lambda 0.25
shapiro.test((PSFmic$PSFu+1)^0.25) #p=0.27

shapiro.test(PSFmic$root_shoot)#p<0.05 
shapiro.test(log(PSFmic$root_shoot))#p<0.05
shapiro.test(sqrt(PSFmic$root_shoot))#p<0.05
#try transform tukey 
transformTukey(PSFmic$root_shoot+2)#lambda -1.85
shapiro.test(-1*(PSFmic$root_shoot+2)^-1.85)#p=0.1

#normalize and mean center response data 
stnorm <- function(x) (x-mean(x,na.rm=TRUE))/ (sqrt(sum((x-mean(x,na.rm=TRUE))^2,na.rm=TRUE)/(length(x)-1)))
PSFmic$height<-stnorm(log(PSFmic$Final.height..cm.))
PSFmic$PSF_u<-stnorm((PSFmic$PSFu+1)^0.25)
PSFmic$root_shootz<-stnorm(-1*(PSFmic$root_shoot+2)^-1.85)
PSFmic$SLA<-stnorm(sqrt(PSFmic$SLA..cm2.g.))
PSFmic$LDMC<-stnorm(log(PSFmic$LDMC..g.g.))
PSFmic$leafCN<-stnorm(log(PSFmic$C.N.ratio))
PSFmic$leafC<-stnorm(PSFmic$Wt..C.)
PSFmic$leafN<-stnorm(PSFmic$Wt..N.)
PSFmic$leaf15N<-stnorm(PSFmic$X.15N)
PSFmic$leaf13C<-stnorm(sqrt(PSFmic$X.13C+36))
#PC 1-4 already mean centered and normal 

#mean center predictors 
PSFmic$saprich = stnorm(PSFmic$saprich)
PSFmic$mprich = stnorm(PSFmic$mprich) 
PSFmic$mutrich = stnorm(PSFmic$mutrich) 
PSFmic$pathrich = stnorm(PSFmic$pathrich) 
PSFmic$pcoa1<-PSFmic$Dim1
PSFmic$pcoa2<-PSFmic$Dim2
PSFmic$chao = stnorm(PSFmic$chao1)
PSFmic$CBH = stnorm(PSFmic$CBH.activity.inoc)
PSFmic$NAG = stnorm(PSFmic$NAG.activity.inoc)

PSFmic<-dplyr::select(PSFmic, chao, CBH, NAG, pcoa1,pcoa2,mprich, saprich,mutrich, pathrich,
  height, PSF_u, root_shootz, PC1, PC2, PC3,PC4, SLA, LDMC, leafCN, leafC, leafN, 
  leaf13C, leaf15N, soil.elev, soil.treat, Seed.elev, elev.treat, Pot)
#save(PSFmic,file= "Data/GH_lmer.RData")

#microbial data----
#raw microbe data for ggplots 
#Alpha diversity 
alphadiv<-read.csv("Data/alphadiv_averaged.csv")
GH1<-unite(GH1, elev.treat, soil.elev, soil.treat, remove=F)
join<-dplyr::select(GH1, elev.treat, Pot)
#calculate averages for each soil x elev type across 5 samples of each 
alphadiv<-left_join(alphadiv, join)

#Funguild 
#from funguild_reclaculate.R
funguild<-read.csv("Data/funguild_abund_rich_bysample.csv")%>%rename(ID=sample)%>%dplyr::select(-X)
metadataGH<-read.csv("Data/map_SampleSeq_info.csv")%>%rename(ID=X.SampleID)
funguild_meta<-left_join(funguild, metadataGH)

funguild<-left_join(funguild_meta, join)

#community comp
PCoA_bray_coords<-read.csv("Data/PCoA_bray_coords.csv")%>%rename(ID=X)
PCoA_meta<-left_join(PCoA_bray_coords, metadataGH)
PCoA_bray_coords<-left_join(PCoA_meta, join)                                 

microbe<-left_join(alphadiv, funguild)
microbe<-left_join(microbe, PCoA_bray_coords)
microbe<-left_join(microbe, dplyr::select(PSF, Pot, NAG.activity.inoc, CBH.activity.inoc))
microbe$X<-NULL

shapiro.test(log(microbe$CBH))#p=0.9 log transform= normal
shapiro.test(log(microbe$NAG))#p=0.9 log transform =normal
shapiro.test(log(microbe$chao1))#p=0.4 log transform =normal
shapiro.test(sqrt(microbe$mprich))#p=0.07 sqrt transform= normal
shapiro.test(microbe$saprich)#p=0.8 no transform= normal
shapiro.test(microbe$mutrich) #p=0.7 normal no transform 
shapiro.test(microbe$pathrich) #p=0.06 normal no transform 

#normalize and mean center 
microbe$CBH = stnorm(log(microbe$CBH))
microbe$NAG = stnorm(log(microbe$NAG))
microbe$chao = stnorm(log(microbe$chao1))
microbe$mprich = stnorm(sqrt(microbe$mprich)) 
microbe$saprich = stnorm(microbe$saprich)
microbe$mutrich = stnorm(microbe$mutrich)
microbe$pathrich = stnorm(log(microbe$pathrich))

#save(microbe, PSFmic,file= "Data/data_anovas.RData")
