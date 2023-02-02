## load libraries
library(permute)
library(indicspecies)
library(tidyr)
library(dplyr)
library(ggplot2)

#run inicsp analysis----
#ALL OTUs 
#has taxonomy column and funguild assignments 
ITS<-read.csv("Data/bioinformatics/Courtney.guilds_indicsp.csv") #2797 OTUs 

#load mapping file, remove 'sample ID' cell before
metadataGH<-read.csv("Data/map_SampleSeq_info.csv")
metadataGH$ID<-metadataGH$X.SampleID
rownames(metadataGH)<-metadataGH[,1]
map<-metadataGH[,-1]

#expand taxonomy column 
ITS<-separate(ITS, Taxonomy, c("x","x1", "x2", "x3","x4", "x5", "k", "k1", "p", "p1","c", "c1","o", "o1","f", "f1","g", "g1","s", "s1", "s2", "s3"), remove=FALSE)
ITS2<-filter(ITS, x3=="None")%>%separate(Taxonomy, c("x","x1", "x2", "x3", "k", "k1", "p", "p1","c", "c1","o", "o1"), remove=FALSE)
ITS<-filter(ITS, x3!="None")

#filter out anything not fungi 
ITS_Fungi<-rbind(ITS2, ITS)%>%filter(k1=="Fungi")%>%select(-x,-x1,-x2, -x3,-x4,-x5, -k, -p, -c, -o,-f,-g,-s)


#remove last rows, taxonomy, from otu table
otu<-ITS_Fungi[,c(1:41)]
str(otu)
#make new dataframe for taxa-2470 Fungal OTUs
taxa<-(ITS_Fungi[,c(1,42:51)])
taxa$X.OTU.ID<-as.integer(taxa$X.OTU.ID)



#transpose OTU table so that samples are rows and OTU IDs are columns
OTU.transp<-t(otu)
OTU.transp<-as.data.frame(OTU.transp)
str(OTU.transp)
colnames(OTU.transp)<-(OTU.transp[1,])
OTU.transp=OTU.transp[-1,]

#Discard all samples in one matrix that don't appear in the other
map.sub<-map[row.names(map) %in% row.names(OTU.transp),]
otu.transp.sub<-OTU.transp[row.names(OTU.transp) %in% row.names(map),]

#match order of mapping file with order of samples in OTU table
p<-match(row.names(OTU.transp),row.names(map))
map2<-map[p , ]

#check
nrow(map2)==nrow(OTU.transp)

#create factor name to call for Uniqueclass i.e. cluster in indicspp
elev<-as.factor(map2$Soilelev)
trt<-as.factor(map2$Soiltreat)
map2<-unite(map2, elev_treat, Soilelev, Soiltreat, remove=FALSE)
elev_trt<-as.factor(map2$elev_treat)

#indicator species analysis by unique soil elev x treat identifier
indval_elev_trt<-multipatt(x=as.data.frame(OTU.transp),
                       cluster=elev_trt,
                       func="IndVal.g",
                       duleg=TRUE,
                       control=how(nperm=999))
#Get summary information of OTUs that are indicators
#reorganize output and combine with taxonomy info
summary(indval_elev_trt)
elev_trt_out<-indval_elev_trt$sign
elev_trt_out$X.OTU.ID<-as.integer(row.names(elev_trt_out))
elev_trt_out<-mutate(elev_trt_out, elev_trt=case_when(index==1~"3100_GC", index==2~"3100_SC",index==3~"3100_SR1",
index==4~"3100_SR5", index==5~"3700_GC", index==6~"3700_SC", index==7~"3700_SR1",index==8~"3700_SR5"))%>%
  left_join(.,taxa)%>%filter(p.value<0.05)
elev_trt_out<-elev_trt_out[,-c(1:9)]
#filter out anything not at least assigned to class 
elev_trt_out<-filter(elev_trt_out, !is.na(p1))


#post analysis----
elev_trt_out$X<-NULL
library(FUNGuildR)

unique(elev_trt_out$elev_trt)
str(elev_trt_out)

check<-group_by(elev_trt_out, elev_trt, c1)%>%summarise(n=n())%>%group_by(elev_trt)%>%mutate(tot=sum(n))

group_by(elev_trt_out, elev_trt)%>%summarise(n=n())

post_indic_class<-group_by(elev_trt_out, elev_trt, p1, c1)%>%summarise(count=n())%>%
  filter(!is.na(c1))
unique(post_indic_class$c1)
unique(post_indic_class$p1)

post_indic_ord<-group_by(elev_trt_out, elev_trt, p1, c1, o1)%>%summarise(count=n())%>%
  filter(!is.na(o1))
unique(post_indic_ord$o1)


post_indic_spp<-filter(elev_trt_out, !is.na(s1))

#stacked bar chart 
specColor <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
  "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99","#FFD700", "#C4A484","#B15928",
  "#FF80FF" , "#CE50CA","#787f63", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7",
  "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
  "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D",
  "#8A7C64", "#599861", "#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7",
  "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
  "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D",
  "#8A7C64", "#599861")

ggplot(check, aes(x=elev_trt, y=n,fill=c1)) + geom_bar(stat = 'identity', position="stack")+
  theme_classic()+ ylab("indic spp (n)")+ xlab("soil elev x soil source")+
  scale_fill_manual(values = specColor, na.value = "grey")+
  geom_text(aes(label=n),size=3, position = position_stack(vjust=0.5))


str(ITS)
funguild<-select(ITS, Taxonomy, Taxon, Taxon.Level, Trophic.Mode, Guild,
                 Confidence.Ranking, Growth.Morphology, Trait, Notes, Citation.Source)
#expand taxonomy column 
funguild<-separate(funguild, Taxonomy, c("x","x1", "x2", "x3","x4", "x5", "k", "k1", "p", "p1","c", "c1","o", "o1","f", "f1","g", "g1","s", "s1", "s2", "s3"), remove=FALSE)
funguild2<-filter(funguild, x3=="None")%>%separate(Taxonomy, c("x","x1", "x2", "x3", "k", "k1", "p", "p1","c", "c1","o", "o1"), remove=FALSE)
funguild<-filter(funguild, x3!="None")

#filter out anything not fungi 
funguild<-rbind(funguild2, funguild)%>%filter(k1=="Fungi")%>%select(-x,-x1,-x2, -x3,-x4,-x5, -k, -p, -c, -o,-f,-g,-s)
str(post_indic_spp)
post_indic_spp<-left_join(post_indic_spp, funguild)
post_indic_ordx<-left_join(post_indic_ord, funguild)

elev_trt_out<-left_join(elev_trt_out, funguild)%>%distinct(.)
elev_trt_out$Taxonomy<-NULL
