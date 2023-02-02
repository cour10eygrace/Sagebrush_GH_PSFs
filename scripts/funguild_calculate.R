#calculate funguild relative abundances by sample for the following criteria 
#Mutaualists only AMF 
#Pathogens only "Plant-Pathogen" guild, none with multiple assignments
#Saprotrophs only "Soil saprotroph" guild or specific mention of soil inhabiting in Notes or Citation 
#both OTU richness and relative abundance 

library(tidyr)
library(dplyr)

#ALL OTUs
#has taxonomy column and funguild assignments 
ITS<-read.csv("Data/Courtney.guilds_indicsp.csv") #2797 OTUs 

#load mapping file, remove 'sample ID' cell before

metadataGH<-read.csv("Data/map_SampleSeq_info.csv")
metadataGH$ID<-metadataGH$X.SampleID4
rownames(metadataGH)<-metadataGH[,1]
map<-metadataGH[,-1]

#expand taxonomy column--ignore warnings
ITS<-separate(ITS, Taxonomy, c("x","x1", "x2", "x3","x4", "x5", "k", "k1", "p", "p1","c", "c1","o", "o1","f", "f1","g", "g1","s", "s1", "s2", "s3"), remove=FALSE)
ITS2<-filter(ITS, x3=="None")%>%separate(Taxonomy, c("x","x1", "x2", "x3", "k", "k1", "p", "p1","c", "c1","o", "o1"), remove=FALSE)
ITS<-filter(ITS, x3!="None")

#filter out anything not fungi #2470
ITS_Fungi<-rbind(ITS2, ITS)%>%filter(k1=="Fungi")%>%dplyr::select(-x,-x1,-x2, -x3,-x4,-x5, -k, -p, -c, -o,-f,-g,-s)
#calculate only funguild assigned OTUs #749 
ITS_Funguild<-subset(ITS_Fungi, Guild!="-")

#calculate path, mututalist, sap richness
path<-filter(ITS_Fungi, Guild=="Plant Pathogen")
path2<-path[,c(1:41)]
path2<-path2[, -1]
path_by_sample<-path2 %>%gather(key = 'sample', value='value')%>%group_by(sample)%>%mutate(rich=ifelse(value>0, 1, 0))%>%
                summarise(pathabund=sum(value), pathrich=sum(rich))

mut<-filter(ITS_Fungi, Guild=="Arbuscular Mycorrhizal")
mut2<-mut[,c(1:41)]
mut2<-mut2[, -1]
mut_by_sample<-mut2 %>%gather(key = 'sample', value='value')%>%group_by(sample)%>%mutate(rich=ifelse(value>0, 1, 0))%>%
  summarise(mutabund=sum(value), mutrich=sum(rich))

sap<-filter(ITS_Fungi, Trophic.Mode=="Saprotroph"& !grepl("Dung", Guild) & !grepl("Parasite", Guild))
sap2<-filter(sap,grepl("soil", Notes)|grepl("Soil", Guild)|grepl("Sterkenburg", Citation.Source))
sap3<-filter(sap, f1=="Hypocreaceae"|f1=="Chaetosphaeriaceae"|f1=="Hyaloscyphaceae"|
               f1=="Herpotrichiellaceae" | f1=="Mycenaceae" |f1=="Marasmiaceae"|f1=="Thelephoraceae"|
               g1=="Wallemia"|o1=="Pezizales"|g1=="Clavulinopsis") #Taxa from Purahong 2019, Tedersoo 2008, Jancic 2015 
sap<-rbind(sap2, sap3)            
sap2<-sap[,c(1:41)]
sap2<-sap2[, -1]
sap_by_sample<-sap2%>%gather(key = 'sample', value='value')%>%group_by(sample)%>%mutate(rich=ifelse(value>0, 1, 0))%>%
  summarise(sapabund=sum(value), saprich=sum(rich))

funguild_abund_rich<-left_join(sap_by_sample, path_by_sample)
funguild_abund_rich<-left_join(funguild_abund_rich, mut_by_sample)
funguild_abund_rich$mpratio<-funguild_abund_rich$mutabund/funguild_abund_rich$pathabund
funguild_abund_rich$mprich<-funguild_abund_rich$mutrich/funguild_abund_rich$pathrich

#write.csv(funguild_abund_rich, "Data/funguild_abund_rich_bysample.csv")
