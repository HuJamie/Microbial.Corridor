######
##### Pool three replicates together, to check the population level diversity and composition ####

ngs.pool<-aggregate(t(ngs.rarefied), list(env$month,env$pool.code,env$year,env$code,env$position,env$corridor,env$mesocosm), FUN=mean)    ## pool data

env.pool<-ngs.pool[,1:7]
ngs.pool1<-ngs.pool[,-c(1:7)]
ngs.pool2<-as.data.frame(t(ngs.pool1))

library(vegan)
otu.shannon<- diversity(ngs.pool1,index="shannon") 
otu.richness<-rowSums(ngs.pool1>0)
otu.evenness<- otu.shannon/log(otu.richness)                  ## Pielou's evenness
#otu.alpha<-fisher.alpha(ngs.pool2,MARGIN=2)                   ## function accepts only integers (counts)

n<-rowSums(ngs.pool1==1)                                    ## number of singletons
m<-rowSums(ngs.pool1==2)                                    ## number of doubletons
otu.Chao1<-otu.richness+n*(n-1)/(2*(m+1))                   ## calculate chao1 estimate 

richness.pool<-cbind(env.pool,otu.richness,otu.shannon,otu.evenness)

write.table(richness.pool,file="richness.pool.txt",sep="\t")
richness.poolf<-read.table("richness.poolf.txt", sep="\t", header=T,row.names=1)    ## sample properties, we need this

richness.pp<-subset(richness.poolf,richness.poolf$position!="mb"&richness.poolf$position!="cor")
richness.pp<-subset(richness.pool,richness.pool$Group.5!="mb"&richness.pool$Group.5!="cor")
colnames(richness.pp)<-c("month","sample","year","code","position","corridor","meso","otu.richness","otu.shannon","otu.evenness")

library(gridExtra)
library(ggplot2)
pdf("Community Richness and Evenness of population.pdf",height=5,width=11)

pool.r<-ggplot(data=richness.pp,aes(x=Group.4,y=otu.richness))+
  geom_boxplot(aes(fill=as.factor(Group.6)))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.6,color=as.factor(Group.6)))+ 
  xlab("Different campaign") + 
  ylab("Pool Sample OTU richness")

pool.e<-ggplot(data=richness.pp,aes(x=Group.4,y=otu.evenness))+
  geom_boxplot(aes(fill=as.factor(Group.6)))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.6,color=as.factor(Group.6)))+ 
  xlab("Different Campaign") + 
  ylab("Pool Sample OTU evenness")

grid.arrange(pool.r,pool.e,ncol=2)
dev.off()  

richness.p2017<-subset(richness.poolf,richness.poolf$year==2017)
richness.p2018<-subset(richness.poolf,richness.poolf$year==2018)
richness.p2018.may<-subset(richness.poolf,richness.poolf$year==2018 &richness.poolf$month== " may")
richness.p2018.june<-subset(richness.poolf,richness.poolf$year==2018 &richness.poolf$month== "june")
richness.p2018.oct<-subset(richness.poolf,richness.poolf$year==2018 &richness.poolf$month== "october")
richness.p2019<-subset(richness.poolf,richness.poolf$year==2019)

####Bray-Curtis dissimilarity at patch level for all root mycobiota ####
#### Distance among population 1 and 2 ####
ngs.poolf<-cbind(richness.poolf[,1:6],ngs.pool1)
ngs.poolpp<-subset(ngs.poolf,ngs.poolf$position!="cor"&ngs.poolf$position!="mb")

#### 2017 ####
ngs.pool.2017<-subset(ngs.poolpp,ngs.poolpp$year==2017)
ngs.pool.2017.0<-subset(ngs.pool.2017,ngs.pool.2017$corridor==0)
ngs.pool.2017.1<-subset(ngs.pool.2017,ngs.pool.2017$corridor==1)

richness.p2017.0<-subset(richness.p2017,richness.p2017$corridor==0&richness.p2017$position!="mb"&richness.p2017$position!="cor")
richness.p2017.1<-subset(richness.p2017,richness.p2017$corridor==1&richness.p2017$position!="mb"&richness.p2017$position!="cor")

## calculate the distance of two population in the same mesocosm
p1<-as.vector(t(ngs.pool.2017.1[3,-c(1:6)]))
p2<-as.vector(t(ngs.pool.2017.1[4,-c(1:6)]))
dists<-vegdist(rbind(p1,p2),method="bray")

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
dist.2017<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(ngs.pool.2017.1)[1])
{
  p1<-as.vector(t(ngs.pool.2017.1[i,-c(1:6)]))
  p2<-as.vector(t(ngs.pool.2017.1[i+1,-c(1:6)]))
  dist.2017[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}


for (i in 1:dim(ngs.pool.2017.0)[1])
{
  p1<-as.vector(t(ngs.pool.2017.0[i,-c(1:6)]))
  p2<-as.vector(t(ngs.pool.2017.0[i+1,-c(1:6)]))
  dist.2017[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

dist.2017f<-dist.2017[ c(TRUE,FALSE), ]  # rows
colnames(dist.2017f)<-c("2017 October With corridor","2017 October Without corridor")

hist(dist.2017f)

##### 2018 May ####
ngs.pool.2018m<-subset(ngs.poolf,ngs.poolf$year==2018&ngs.poolf$month==" may")
ngs.pool.2018m.0<-subset(ngs.pool.2018m,ngs.pool.2018m$corridor==0&ngs.pool.2018m$position!="mb"&ngs.pool.2018m$position!="cor")
ngs.pool.2018m.1<-subset(ngs.pool.2018m,ngs.pool.2018m$corridor==1&ngs.pool.2018m$position!="mb"&ngs.pool.2018m$position!="cor")

richness.p2018m.0<-subset(richness.p2018.may,richness.p2018.may$corridor==0&richness.p2018.may$position!="mb"&richness.p2018.may$position!="cor")
richness.p2018m.1<-subset(richness.p2018.may,richness.p2018.may$corridor==1&richness.p2018.may$position!="mb"&richness.p2018.may$position!="cor")

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
dist.2018m<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(ngs.pool.2018m.1)[1])
{
  p1<-as.vector(t(ngs.pool.2018m.1[i,-c(1:6)]))
  p2<-as.vector(t(ngs.pool.2018m.1[i+1,-c(1:6)]))
  dist.2018m[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}


for (i in 1:dim(ngs.pool.2018m.0)[1])
{
  p1<-as.vector(t(ngs.pool.2018m.0[i,-c(1:6)]))
  p2<-as.vector(t(ngs.pool.2018m.0[i+1,-c(1:6)]))
  dist.2018m[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

dist.2018mf<-dist.2018m[ c(TRUE,FALSE), ]  # rows
colnames(dist.2018mf)<-c("2018 May With corridor","2018 May Without corridor")

hist(dist.2018mf)

#### 2018 June ####
ngs.pool.2018j<-subset(ngs.poolf,ngs.poolf$year==2018&ngs.poolf$month=="june")
ngs.pool.2018j.0<-subset(ngs.pool.2018j,ngs.pool.2018j$corridor==0&ngs.pool.2018j$position!="mb"&ngs.pool.2018j$position!="cor")
ngs.pool.2018j.1<-subset(ngs.pool.2018j,ngs.pool.2018j$corridor==1&ngs.pool.2018j$position!="mb"&ngs.pool.2018j$position!="cor")

richness.p2018j.0<-subset(richness.p2018.june,richness.p2018.june$corridor==0&richness.p2018.june$position!="mb"&richness.p2018.june$position!="cor")
richness.p2018j.1<-subset(richness.p2018.june,richness.p2018.june$corridor==1&richness.p2018.june$position!="mb"&richness.p2018.june$position!="cor")


richness.p2018j.pp<-subset(richness.p2018.june,richness.p2018.june$position!="mb"&richness.p2018.june$position!="cor")

ngs.pool.2018j.pp<-subset(ngs.pool.2018j,ngs.pool.2018j$position!="cor"&ngs.pool.2018j$position!="mb")
rda.2018j<-rda(log2(ngs.pool.2018j.pp[,-c(1:5)]+1)~corridor,ngs.pool.2018j.pp) 

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
dist.2018j<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(ngs.pool.2018j.1)[1])
{
  p1<-as.vector(t(ngs.pool.2018j.1[i,-c(1:6)]))
  p2<-as.vector(t(ngs.pool.2018j.1[i+1,-c(1:6)]))
  dist.2018j[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

for (i in 1:dim(ngs.pool.2018j.0)[1])
{
  p1<-as.vector(t(ngs.pool.2018j.0[i,-c(1:6)]))
  p2<-as.vector(t(ngs.pool.2018j.0[i+1,-c(1:6)]))
  dist.2018j[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

dist.2018jf<-dist.2018j[c(TRUE,FALSE), ]  # rows
colnames(dist.2018jf)<-c("2018 June With corridor","2018 June Without corridor")

hist(dist.2018jf)

#### 2018 October ####
ngs.pool.2018o<-subset(ngs.poolf,ngs.poolf$year==2018&ngs.poolf$month=="october")
ngs.pool.2018o.0<-subset(ngs.pool.2018o,ngs.pool.2018o$corridor==0&ngs.pool.2018o$position!="mb"&ngs.pool.2018o$position!="cor")
ngs.pool.2018o.1<-subset(ngs.pool.2018o,ngs.pool.2018o$corridor==1&ngs.pool.2018o$position!="mb"&ngs.pool.2018o$position!="cor")

richness.p2018o.0<-subset(richness.p2018.oct,richness.p2018.oct$corridor==0&richness.p2018.oct$position!="mb"&richness.p2018.oct$position!="cor")
richness.p2018o.1<-subset(richness.p2018.oct,richness.p2018.oct$corridor==1&richness.p2018.oct$position!="mb"&richness.p2018.oct$position!="cor")


richness.p2018o.pp<-subset(richness.p2018.oct,richness.p2018.oct$position!="mb"&richness.p2018.oct$position!="cor")

ngs.pool.2018o.pp<-subset(ngs.pool.2018o,ngs.pool.2018o$position!="cor"&ngs.pool.2018o$position!="mb")
rda.2018o<-rda(log2(ngs.pool.2018o.pp[,-c(1:5)]+1)) 

library(BiodiversityR)
plot.otu.rda<-plot(rda.2018o,scaling=2, label=FALSE)
ordisymbol(plot.otu.rda,richness.p2018o.pp,"corridor",col=1,pchs=F,legend=T,legend.x="topright",legend.ncol=1)

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
dist.2018o<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(ngs.pool.2018o.1)[1])
{
  p1<-as.vector(t(ngs.pool.2018o.1[i,-c(1:6)]))
  p2<-as.vector(t(ngs.pool.2018o.1[i+1,-c(1:6)]))
  dist.2018o[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

for (i in 1:dim(ngs.pool.2018o.0)[1])
{
  p1<-as.vector(t(ngs.pool.2018o.0[i,-c(1:6)]))
  p2<-as.vector(t(ngs.pool.2018o.0[i+1,-c(1:6)]))
  dist.2018o[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

dist.2018of<-dist.2018o[ c(TRUE,FALSE), ]  # rows
colnames(dist.2018of)<-c("2018 October With corridor","2018 October Without corridor")

hist(dist.2018of)

#####2019####
ngs.pool.2019<-subset(ngs.poolf,ngs.poolf$year==2019)
ngs.pool.2019.0<-subset(ngs.pool.2019,ngs.pool.2019$corridor==0&richness.p2019$position!="mb"&richness.p2019$position!="cor")
ngs.pool.2019.1<-subset(ngs.pool.2019,ngs.pool.2019$corridor==1&richness.p2019$position!="mb"&richness.p2019$position!="cor")

richness.p2019.0<-subset(richness.p2019,richness.p2019$corridor==0&richness.p2019$position!="mb"&richness.p2019$position!="cor")
richness.p2019.1<-subset(richness.p2019,richness.p2019$corridor==1&richness.p2019$position!="mb"&richness.p2019$position!="cor")
richness.p2019.pp<-subset(richness.p2019,richness.p2019$position!="mb"&richness.p2019$position!="cor")

ngs.pool.2019.pp<-subset(ngs.pool.2019,ngs.pool.2019$position!="cor"&ngs.pool.2019$position!="mb")

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
dist.2019<-as.data.frame(matrix(NA,nrow=20,ncol=2))

for (i in 1:dim(ngs.pool.2019.1)[1])
{
  p1<-as.vector(t(ngs.pool.2019.1[i,-c(1:6)]))
  p2<-as.vector(t(ngs.pool.2019.1[i+1,-c(1:6)]))
  dist.2019[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

for (i in 1:dim(ngs.pool.2019.0)[1])
{
  p1<-as.vector(t(ngs.pool.2019.0[i,-c(1:6)]))
  p2<-as.vector(t(ngs.pool.2019.0[i+1,-c(1:6)]))
  dist.2019[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

dist.2019f<-dist.2019[ c(TRUE,FALSE), ]  # rows
colnames(dist.2019f)<-c("2019 May Without corridor","2019 May With corridor")
hist(dist.2019f)

dist.all<-cbind(dist.2017f,dist.2018mf,dist.2018jf,dist.2018of,dist.2019f)

dist.allf<-stack(dist.all)
boxplot(values~ind,dist.allf)

write.table(dist.allf,file="dist.all.txt",sep="\t")

dist.alle<-read.table("dist.all20210421.txt", sep="\t", header=T)        ## sample properties, we need this
dist.allc<-read.table("dist.alle.txt", sep="\t", header=T)        ## sample properties, we need this

distance.e<-ggplot(data=dist.alle,aes(x=code,y=distance))+ylim(0,1)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(corridor)))+ 
  geom_boxplot(aes(fill=as.factor(corridor)))+
  xlab("Different campagin") + ylab("Diatance between P1 and P2")

distance.c<-ggplot(data=dist.allc,aes(x=code,y=distance))+ylim(0,1)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(corridor)))+ 
  geom_boxplot(aes(fill=as.factor(corridor)))+
  xlab("Different campagin") + ylab("Diatance between P1 and P2")

grid.arrange(distance.e,distance.c,ncol=2)

ggplot(data=dist.alle,aes(x=code,y=distance))+geom_jitter(width=0.2)+geom_boxplot(aes(fill=as.factor(corridor)))+
  facet_wrap(~corridor,ncol=2)

lmd<-lm(distance~code+corridor,dist.allc)
anova(lmd)
summary(lmd)

dist.alle.1<-subset(dist.alle,dist.alle$corridor=="With")
dist.alle.0<-subset(dist.alle,dist.alle$corridor=="Without")

ggplot(data=dist.alle.1,aes(x=code,y=distance))+ylim(0,1)+
  geom_jitter(width=0.18)+
  geom_smooth()+
  xlab("Different campaigns") + ylab("Bray-Curtis dissimilarity between P1 and P2 with corridor")

ggplot(data=dist.alle.0,aes(x=code,y=distance))+ylim(0,1)+
  geom_jitter(width=0.18)+
  xlab("Different campaigns") + ylab("Bray-Curtis dissimilarity between P1 and P2 without corridor")

lmd1<-lm(distance~code,dist.alle.1)
anova(lmd1)
summary(lmd1)

lmd0<-lm(distance~code,dist.alle.0)
anova(lmd0)
summary(lmd0)


boxplot(distance~campagin,dist.allf.1)
boxplot(distance~campagin,dist.allf.0)
ggplot(data=dist.allf,aes(x=campagin,y=distance))+geom_point()+
  facet_wrap(~corridor,ncol=2)

distance<-ggplot(data=dist.allflog,aes(x=campagin,y=distance))+ylim(5,10)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(corridor)))+ 
  geom_boxplot(aes(fill=as.factor(corridor),color=as.factor(campagin)))+
   xlab("Different campagin") + ylab("Diatance between P1 and P2")

distance.0<-ggplot(data=dist.allflog.0,aes(x=campagin,y=distance))+ylim(5,10)+
  geom_point(position=position_dodge(width=0.75))+ 
  geom_boxplot(aes(fill=as.factor(campagin)))+
  xlab("Different campagin") + ylab("Diatance between P1 and P2 without corridor")

distance.1<-ggplot(data=dist.allflog.1,aes(x=campagin,y=distance))+ylim(5,10)+
  geom_point(position=position_dodge(width=0.75))+ 
  geom_boxplot(aes(fill=as.factor(campagin)))+
  xlab("Different campagin") + ylab("Diatance between P1 and P2 with corridor")


library(gridExtra)
pdf("Distance between population 1 and population 2.pdf",height=4,width=10)
grid.arrange(distance.0,distance.1,ncol=2)
dev.off()

pdf("Distance between populations together.pdf",height=4,width=10)
distance
dev.off()


lm1<-lm(distance~campagin+corridor,dist.allflog)
anova(lm1)
summary(lm1)

lm2<-lm(distance~campagin,dist.allflog.0)
anova(lm2)

lm3<-lm(distance~campagin,dist.allflog.1)
anova(lm3)

##### pool population 1 and population 2 for meta-population level ####
ngs.t<-as.data.frame(t(ngs.rarefied))
ngs.te<-cbind(env,ngs.t)

ngs.pp<-subset(ngs.te,ngs.te$position!="mb"&ngs.te$position!="cor")
ngs.ppool<-aggregate(ngs.pp[,18:283], list(ngs.pp$mesocosm,ngs.pp$code,ngs.pp$corridor), FUN=mean)    ## pool data

env.ppool<-ngs.ppool[,1:3]
colnames(env.ppool)<-c("macrocosm","campaign","corridor")
ngs.ppool1<-ngs.ppool[,-c(1:3)]
ngs.ppool2<-as.data.frame(t(ngs.ppool1))

library(vegan)
otu.shannon.pp<- diversity(ngs.ppool1,index="shannon") 
otu.richness.pp<-rowSums(ngs.ppool1>0)
otu.evenness.pp<- otu.shannon.pp/log(otu.richness.pp)                  ## Pielou's evenness

richness.ppool<-cbind(env.ppool,otu.richness.pp,otu.shannon.pp,otu.evenness.pp)

richness.pp2017<-subset(richness.ppool,richness.ppool$campaign=="17O")
richness.pp2018.may<-subset(richness.ppool,richness.ppool$campaign=="18 M")
richness.pp2018.june<-subset(richness.ppool,richness.ppool$campaign=="18J")
richness.pp2018.oct<-subset(richness.ppool,richness.ppool$campaign=="18O")
richness.pp2019<-subset(richness.ppool,richness.ppool$campaign=="19M")


library(gridExtra)
pdf("Community Richness and Evenness of meta-population.pdf",height=5,width=11)

pp.richness<-ggplot(data=richness.ppool,aes(x=campaign,y=otu.richness.pp))+
  geom_boxplot(aes(fill=as.factor(corridor)))+
  ylim(120,200)+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylab("OTU Richness")

pp.evenness<-ggplot(data=richness.ppool,aes(x=campaign,y=otu.evenness.pp))+
  geom_boxplot(aes(fill=as.factor(corridor)))+
  ylim(0.48,0.7)+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylab("OTU Pielou Evenness")

grid.arrange(pp.richness,pp.evenness,ncol=2)

dev.off()


pp.shannon<-ggplot(data=richness.ppool,aes(x=campaign,y=otu.shannon.pp))+
  geom_boxplot(aes(fill=as.factor(corridor)))+
  ylim(2.5,3.5)+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylab("OTU Shannon Diversity")
