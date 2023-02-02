#load libraries 
library(ggplot2)
library(dplyr)
library(vegan)
library(ape)
library(tidyr)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

#starting with distance matrix
#fungi 
bray<-read.csv("Data/bray_curtis_dm.csv")
rownames(bray)<-bray$X
bray$X<-NULL
map<-read.csv("Data/map_SampleSeq_info.csv") #sample info
rownames(map)<-map$X.SampleID
map$X.SampleID<-NULL
map$Soiltreat2<-NULL 
pmap<-match(rownames(bray), rownames(map))
map2b<-map[pmap, ]
bray1<-as.dist(as(bray,"matrix"))

#run NMDS
#fungi
NMDS_bray<-metaMDS(bray, k=3, trymax=200)#stress=0.124
NMDS_bray_coords<-as.data.frame(NMDS_bray$points)
NMDS_bray_coords<-merge(NMDS_bray_coords, map2b, by="row.names")

#run adonis in Vegan 
#fungi 
adonis(bray1~  Soiltreat,data=map2b,permutations=999, strata =  map2b$Soilelev)
adonis(bray1~  Soiltreat*Soilelev ,data=map2b,permutations=999)

#pairwise
pairwise.adonis2(bray1 ~ as.factor(Soiltreat),data=map2b, strata = 'Soilelev', nperm=999)
pairwise.adonis2(bray1 ~ as.factor(Soilelev),data=map2b, strata = 'Soiltreat', nperm=999)

#combine soil elevation and soil treat 
map2b$Soilelev<-as.character(as.factor(map2b$Soilelev))
map2b<-unite(map2b, 'elevtrt', c(Soilelev, Soiltreat), remove = F)

pairwise.adonis2(bray1 ~ elevtrt,data=map2b,  nperm=999)


#ggplot NMDS-FIGURE 3
#library(ggplot2)
#source_url("https://raw.github.com/low-decarie/FAAV/master/r/stat-ellipse.R")   
ggplot(data=NMDS_bray_coords, aes(x=MDS1,y=MDS2,colour=Soiltreat))+  
  geom_point(aes(x=MDS1,y=MDS2,colour=Soiltreat,  shape=as.factor(Soilelev)), size=2.5)+ 
  stat_ellipse(lty=2)+ theme_classic()  + scale_color_brewer(palette = "Paired", direction = -1)

#betadisper
bray.disp.trt<-betadisper(bray1,map2b$Soiltreat)
anova(bray.disp.trt)
plot(bray.disp.trt)
permutest(bray.disp.trt, pairwise = TRUE, permutations = 999)
 
#PCoA 
PCoA_bray<-ape::pcoa(bray)#run pcoa ape
PCoA_bray$values$Relative_eig[1:3] #percent explained by each axis 
biplot.pcoa(PCoA_bray)

PCoA_fungi<-cmdscale(bray) #run pcoa vegan
#PCoA_fungi<-merge(PCoA_fungi, map2b, by="row.names")
#row.names(PCoA_fungi)<-PCoA_fungi$Row.names
#PCoA_fungi$Row.names<-NULL

#ggplot pcoa Fig S2
#with loading Axes
efit<-envfit(PCoA_fungi, map2b, permutations = 999, strata = NULL, 
             display = "sites", na.rm = FALSE)
efit<-as.data.frame(efit$factors$centroids)
efit$group<-rownames(efit)
efit<-subset(efit, grepl("Soiltreat", group))
efit<-efit%>%mutate(group=stringr::str_remove(group, pattern = "Soiltreat"))




df.plot<-data.frame(PCoA_bray$vectors)

x_label<-round(PCoA_bray$values$Relative_eig[1]*100,2) #percent explained by each axis 

y_label<-round(PCoA_bray$values$Relative_eig[2]*100,2) #percent explained by each axis 


df.plot$group<-map2b$Soiltreat

ggplot(data=df.plot,aes(x=Axis.1,y=Axis.2,
                        
                        
                        color=group))+
  
  
  geom_point(size=5)+
  
  
  theme_bw()+
  

 geom_segment(data=efit, aes(x = 0, y = 0, xend = Dim1, yend = Dim2, colour=group), 
              arrow = arrow())+
 
  geom_text(data = efit, aes(x =Dim1, y = Dim2, label = group),size = 3) + 
  
  
  theme(panel.grid = element_blank())+
  
  
  geom_vline(xintercept = 0,lty="dashed")+
  
  
  geom_hline(yintercept = 0,lty="dashed")+
  
  
  labs(x=paste0("PCoA1 ",x_label,"%"),
       
       
       y=paste0("PCoA2 ",y_label,"%"))  +
  
  
  stat_ellipse(data=df.plot,
               
               
             geom = "polygon",
              
               
            aes(fill=group),
               
               
               alpha=0.1)

