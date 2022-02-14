##############################################
###### Data format, 20210202, Jie hu ####
###### pool p1, p2 and corridor

rm(list = ls()) #remove the data in the environment
ngs.r<-read.table("whitecloverfungi.txt", header=T, row.names=1, sep="\t")  ## otu abundance rawdata
taxo<-read.table("taxo.txt", header=T, row.names=1, sep="\t")  ## otu taxonomic rawdata
taxo.f<-as.data.frame(t(taxo))
env<-read.table("sample_information.txt", sep="\t", header=T,row.names=1)        ## sample properties, we need this
ngs.s<-ngs.r[13:833,-1]  ##selected samples, remove "TOTAL"

rownames(ngs.s)<-rownames(env)
env<-env[order(rownames(env)),]
ngs.s<-ngs.s[order(rownames(ngs.s)),]

ngs.f<-as.data.frame(t(ngs.s))

#### overview, calculate value for method description ####
colSums(ngs.f)
mean(colSums(ngs.f))
min(colSums(ngs.f))
max(colSums(ngs.f))
sum(ngs.f)

## loading packages 
library(vegan) # only the ones that you need
###rarefaction curve for selecting the minimum reads data  ####
#### Figure S1 for PNAS paper ####
ngs.ra<-ngs.r[14:833,]
rarecurve(ngs.ra,step=100,xlim=c(0,6500),ylim=c(0,120),xlab="Sequence depth", ylab="Number of OTUs",col="blue",label=FALSE)
abline(v=4292)
#rareslope(ngs.ra, sample=4292)
#####

#ngs.rarefied<-rarefy(ngs.s, sample=3310, se = FALSE, MARGIN = 1)
#write.table(ngs.rarefied,file="rarefied_matrix.txt",sep="\t")

#ngs.rarefied<-as.data.frame(rarefy(ngs.s,sample=3310))
#ngs.rarefied<-drarefy(ngs.s, sample=3310)

#### calculate relative abundance of each cluster
ngs.relat<-0*ngs.f                                                      ## create a empty matrix
for(i in 1:dim(ngs.relat)[2]) ngs.relat[,i]<-ngs.f[,i]/colSums(ngs.f)[i] ## loop to calculate relative abundance
colSums(ngs.relat)                                                       ## check if calculate correctly

#### Standardization of abundance = rarefaction ####
ngs.b<-round(ngs.relat*min(colSums(ngs.f)),0)   ## recover otu number based relative otu abundance, rarefy
ngs.rt<-as.data.frame(cbind(ngs.b,taxo))              ## combine otu and taxonomy
ngs.rarefied<-ngs.b[,-c(822,823)]               ## subset based on dimensionality of env.txt

mean(colSums(ngs.rarefied))
min(colSums(ngs.rarefied))
max(colSums(ngs.rarefied))
ngs.rarefied.t<-t(ngs.rarefied)

write.table(ngs.rarefied.t,file="rarefied_matrix.txt",sep="\t")

ngs.b$zero.num<-rowSums(ngs.b==0)             ## calculate the zero amount of each row
ngs.b$sum<-rowSums(ngs.b)                     ## calculate sum of each row, it makes some mistakes!!!!

colb<-data.frame(colSums(ngs.b!=0))           ## OTU richness of each sample
ngs.01<-1*(ngs.rarefied>0)                    ## transfer to 0/1 matrix 
ngs.01.t<-1*(ngs.rarefied.t>0)   

#filter<-rowSums(ngs.b.filter)>sum(ngs.b.filter)/10000
#ngs.c<-ngs.b.filter[filter,]           ## select core Clusters in microbiome
#ngs.c01<-data.frame(1*(ngs.c>0))
#tax.c<-taxo.i[filter,]  
#### this script will calculate diversity index for different community 

ngs.cluster <- as.data.frame(t(ngs.rarefied))      ## transposes
ngs.cluster.01 <- as.data.frame(t(ngs.01))         ## transposes

#### Figure S2 for white clover paper ####
ngs.all<-cbind(ngs.cluster,env)
ngs.brachy<-subset(ngs.cluster,ngs.all$plant=="Brachypodium")
ngs.trifo<-subset(ngs.cluster,ngs.all$plant=="Trifolium")
ngs.brachy1<-as.data.frame(t(ngs.brachy))
ngs.trifo1<-as.data.frame(t(ngs.trifo))

ngs.bra.phy<-aggregate(ngs.brachy1, list(taxo$blast_perc_identity), FUN = sum)
ngs.tri.phy<-aggregate(ngs.trifo1, list(taxo$blast_perc_identity), FUN = sum)

sum(rowSums(ngs.tri.phy[,-1]))
sum(rowSums(ngs.bra.phy[,-1]))

library(ggplot2)
library(gridExtra)
ngs.tri.phy$seq<-rowSums(ngs.tri.phy[,-1])
ngs.tri.phy$prop<-rowSums(ngs.tri.phy[,-1])/sum(rowSums(ngs.tri.phy[,-1]))
ngs.tri.phy$phylum<-c("Ascomycota","Basidiomycota","Chytridiomycota","Glomeromycota","Multi-affiliation","GZygomycota")
seq.tri<-ggplot(ngs.tri.phy, aes(x=phylum, y=prop))+ 
  geom_bar(stat="identity", width=0.55,fill="blue")+
  ylim(0,0.6)+
  labs(y='Proportion of sequences in T. repens root',x='')+
  ggprism::theme_prism()+theme(axis.text.x=element_text(angle = 60))

ngs.bra.phy$seq<-rowSums(ngs.bra.phy[,-1])
ngs.bra.phy$prop<-rowSums(ngs.bra.phy[,-1])/sum(rowSums(ngs.bra.phy[,-1]))
ngs.bra.phy$phylum<-c("Ascomycota","Basidiomycota","Chytridiomycota","Glomeromycota","Multi-affiliation","GZygomycota")
seq.mb<-ggplot(ngs.bra.phy, aes(x=phylum, y=prop))+ 
  geom_bar(stat="identity", width=0.55,fill="blue")+
  ylim(0,0.6)+
  labs(y='Proportion of sequences in B. pinnatum root',x='')+
  ggprism::theme_prism()+theme(axis.text.x=element_text(angle = 60))

pdf("Figure S2.pdf",height=6,width=12)
grid.arrange(seq.tri,seq.mb,ncol=2)
dev.off() 

seq.tri1<-ggplot(ngs.tri.phy, aes(x=phylum, y=seq))+ 
  geom_bar(stat="identity", width=0.55,fill="blue")+
  labs(y='Number of sequences in T. repens root',x='')+
  ggprism::theme_prism()+theme(axis.text.x=element_text(angle = 60))

seq.mb1<-ggplot(ngs.bra.phy, aes(x=phylum, y=seq))+ 
  geom_bar(stat="identity", width=0.55,fill="blue")+
  labs(y='Number of sequences in B. pinnatum root',x='')+
  ggprism::theme_prism()+theme(axis.text.x=element_text(angle = 60))

pdf("Figure S2.1.pdf",height=6,width=12)
grid.arrange(seq.tri1,seq.mb1,ncol=2)
dev.off() 
#####

#### calculate diversity index ####

otu.shannon<- diversity(ngs.cluster, index = "shannon") 
otu.richness<-rowSums(ngs.cluster>0)
otu.evenness<- otu.shannon/log(otu.richness)                  ## Pielou's evenness
otu.alpha<-fisher.alpha(ngs.rarefied,MARGIN=2)

n<-rowSums(ngs.cluster==1)                                    ## number of singletons
m<-rowSums(ngs.cluster==2)                                    ## number of doubletons
otu.Chao1<-otu.richness+n*(n-1)/(2*(m+1))                     ## calculate chao1 estimate 

## gathering the parameters together ####
richness<-cbind(env,otu.richness,otu.shannon,otu.alpha,otu.evenness)
write.table(richness,file="Global richness at individual level.txt",sep="\t")

xyplot(otu.richness~plant,richness)
xyplot(otu.shannon~year,richness)

library(ggplot2)
pdf("Richness of plant species .pdf",height=6,width=12)
ggplot(data=richness,aes(x=plant,y=otu.richness))+geom_point()+geom_boxplot()+
  facet_wrap(corridor~code,ncol=5)
dev.off()

lm1<-lm(otu.richness~code+corridor+plant,data=richness)
anova(lm1)
summary(lm1)

###
richness.may<-rbind(richness.2018.may,richness.2019)
richness.october<-rbind(richness.2018.oct,richness.2017)
richness.indip<-subset(richness,richness$position!="cor" &richness$position!="mb")

pdf("Community Richness and Evenness of individual.pdf",height=5,width=11)

indi.r<-ggplot(data=richness.indip,aes(x=code,y=otu.richness))+
  geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  xlab("Different campaign") + 
  ylab("Individual Sample OTU richness")

indi.e<-ggplot(data=richness.indip,aes(x=code,y=otu.evenness))+
  geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  xlab("Different position") + 
  ylab("Individual Sample Pielou's evenness")

grid.arrange(indi.r,indi.e,ncol=2)
dev.off() 

#### figures about fungal community diversity ####
library(ggplot2)
pdf("2017.pdf",height=4,width=6)
ggplot(data=richness.2017,aes(x=position,y=otu.evenness))+geom_point()+
  facet_wrap(~corridor,ncol=2)

ggplot(data=richness.2017,aes(x=position,y=otu.richness))+geom_point()+
  facet_wrap(~corridor,ncol=2)
dev.off()

library(gridExtra)
pdf("Community Richness and Evenness.pdf",height=9,width=28)

R2017oct<-ggplot(data=richness.2017,aes(x=position,y=otu.richness))+
  geom_boxplot(aes(fill=as.factor(corridor)))+
  ylim(40,105)+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  xlab("Different position") + ylab("2017 October OTU richness")
R2018may<-ggplot(data=richness.2018.may,aes(x=position,y=otu.richness))+
  geom_boxplot(aes(fill=as.factor(corridor)))+
  ylim(40,105)+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  xlab("Different position") + ylab("2018 May  OTU richness")
R2018june<-ggplot(data=richness.2018.june,aes(x=position,y=otu.richness))+
  geom_boxplot(aes(fill=as.factor(corridor)))+
  ylim(40,105)+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  xlab("Different position") + ylab("2018 June OTU richness")
R2018oct<-ggplot(data=richness.2018.oct,aes(x=position,y=otu.richness))+
  geom_boxplot(aes(fill=as.factor(corridor)))+
  ylim(40,105)+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  xlab("Different position") + ylab("2018 October  OTU richness")
R2019may<-ggplot(data=richness.2019,aes(x=position,y=otu.richness))+
  geom_boxplot(aes(fill=as.factor(corridor)))+
  ylim(40,105)+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  xlab("Different position") + ylab("2019 May  OTU richness")

R2017oct1<-ggplot(data=richness.2017,aes(x=position,y=otu.evenness))+
  geom_boxplot(aes(fill=as.factor(corridor)))+
  ylim(0.2,0.7)+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  xlab("Different position") + ylab("2017 October OTU evenness")
R2018may1<-ggplot(data=richness.2018.may,aes(x=position,y=otu.evenness))+
  geom_boxplot(aes(fill=as.factor(corridor)))+
  ylim(0.2,0.7)+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  xlab("Different position") + ylab("2018 May  OTU evenness")
R2018june1<-ggplot(data=richness.2018.june,aes(x=position,y=otu.evenness))+
  geom_boxplot(aes(fill=as.factor(corridor)))+
  ylim(0.2,0.7)+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  xlab("Different position") + ylab("2018 June OTU evenness")
R2018oct1<-ggplot(data=richness.2018.oct,aes(x=position,y=otu.evenness))+
  geom_boxplot(aes(fill=as.factor(corridor)))+
  ylim(0.2,0.7)+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  xlab("Different position") + ylab("2018 October  OTU evenness")
R2019may1<-ggplot(data=richness.2019,aes(x=position,y=otu.evenness))+
  geom_boxplot(aes(fill=as.factor(corridor)))+
  ylim(0.2,0.7)+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  xlab("Different position") + ylab("2019 May  OTU evenness")

grid.arrange(R2017oct,R2018may,R2018june,R2018oct,R2019may,R2017oct1,R2018may1,R2018june1,R2018oct1,R2019may1,ncol=5)

dev.off()

ggplot(data=richness.2018.may,aes(x=position,y=otu.evenness))+geom_point()+
  facet_wrap(~corridor,ncol=2)

ggplot(data=richness.2018.oct,aes(x=position,y=otu.evenness))+geom_point()+
  facet_wrap(~corridor,ncol=2)

ggplot(data=richness.2019,aes(x=position,y=otu.evenness))+geom_point()+
  facet_wrap(~corridor,ncol=2)
