#####PCoA analysis 20210422

library(ecodist)
library(vegan)
library(ggplot2)

#### PCoA at individual level for all plants ####
ngs.log<-log2(ngs.rarefied.t+1)

ngs.bray<-vegdist(ngs.log, method = "bray") # dissimilarity matrix using bray-curtis distance indices on the dataset native to vegan
pcoaVS<-pco(ngs.bray,negvals="zero", dround = 0) # if negvals = 0 sets all negative eigenvalues to zero; if = "rm" corrects for negative eigenvalues using method 1 of Legendre and Anderson 1999

plot(pcoaVS$vectors[,2]~pcoaVS$vectors[,1],xlab = "PCoA1", ylab="PCoA2",pch=19,col=env$plant,
     axes = TRUE, main = "PCoA (ecodist) on White clover data")

pcoa<-cbind(pcoaVS$vectors[,2],pcoaVS$vectors[,1],env)

## code from Marine
#p1=p + stat_ellipse(geom = "polygon", type="norm",  level = 0.95,alpha=0.2, aes(fill=Pratique))

#### Figure. 2A all plants####
pdf("Figure. 2A PCoA of all samples 20210426 plant norm.pdf",height=6,width=8)
ggplot(data=pcoa,aes(x=pcoaVS$vectors[, 1],y=pcoaVS$vectors[, 2]))+
  xlim(-0.035,0.02)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(plant)),size=3)+ 
  stat_ellipse(geom="polygon",type="t",level=0.95,alpha=0.2,aes(fill=plant))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens ans Brachypodium")
dev.off()

adonis.plant<-adonis(ngs.log~plant,env,by=NULL, method="bray") 
adonis.plant

pcoaVS$values # eigenvalue for each component. This is a measure of the variance explained by each dimension
pcoaVS$vectors # eigenvectors. Each column contains the scores for that dimension.

#####Different sampling campaigns with both plants####
pcoa.17o<-subset(pcoa,pcoa$code=="17O")
ngs.17o<-subset(ngs.rarefied.t,env$code=="17O")
p.17o<-ggplot(data=pcoa.17o,aes(x=pcoa.17o[,2],y=pcoa.17o[,1]))+
  xlim(-0.04,0.02)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(plant)),size=3)+ 
  stat_ellipse(geom="polygon",type="t",level=0.95,alpha=0.2,aes(fill=plant))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens ans Brachypodium in 17O")

adonis.17o<-adonis(ngs.17o~plant,pcoa.17o,by=NULL, method="bray") 
adonis.17o

pcoa.18m<-subset(pcoa,pcoa$code=="18 M")
ngs.18m<-subset(ngs.log,env$code=="18 M")
p.18m<-ggplot(data=pcoa.18m,aes(x=pcoa.18m[,2],y=pcoa.18m[,1]))+
  xlim(-0.04,0.02)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(plant)),size=3)+ 
  stat_ellipse(geom="polygon",type="t",level=0.95,alpha=0.2,aes(fill=plant))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens ans Brachypodium in 18M")
adonis.18m<-adonis(ngs.18m~plant,pcoa.18m,by=NULL, method="bray") 
adonis.18m

pcoa.18j<-subset(pcoa,pcoa$code=="18J")
ngs.18j<-subset(ngs.log,env$code=="18J")
p.18j<-ggplot(data=pcoa.18j,aes(x=pcoa.18j[,2],y=pcoa.18j[,1]))+
  xlim(-0.04,0.02)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(plant)),size=3)+ 
  stat_ellipse(geom="polygon",type="t",level=0.95,alpha=0.2,aes(fill=plant))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens ans Brachypodium in 18J")
adonis.18j<-adonis(ngs.18j~plant,pcoa.18j,by=NULL, method="bray") 
adonis.18j

pcoa.18o<-subset(pcoa,pcoa$code=="18O")
ngs.18o<-subset(ngs.log,env$code=="18O")
p.18o<-ggplot(data=pcoa.18o,aes(x=pcoa.18o[,2],y=pcoa.18o[,1]))+
  xlim(-0.04,0.02)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(plant)),size=3)+ 
  stat_ellipse(geom="polygon",type="t",level=0.95,alpha=0.2,aes(fill=plant))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens ans Brachypodium in 18O")
adonis.18o<-adonis(ngs.18o~plant,pcoa.18o,by=NULL, method="bray") 
adonis.18o

pcoa.19m<-subset(pcoa,pcoa$code=="19M")
ngs.19m<-subset(ngs.log,env$code=="19M")
p.19m<-ggplot(data=pcoa.19m,aes(x=pcoa.19m[,2],y=pcoa.19m[,1]))+
  xlim(-0.04,0.02)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(plant)),size=3)+ 
  stat_ellipse(geom="polygon",type="t",level=0.95,alpha=0.2,aes(fill=plant))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens ans Brachypodium in 19M")
adonis.19m<-adonis(ngs.19m~plant,pcoa.19m,by=NULL, method="bray") 
adonis.19m

#######Figure S3####
pdf("Figure S3. all plants at different sampling campaigns.pdf",height=5,width=33)
grid.arrange(p.17o,p.18m,p.18j,p.18o,p.19m,ncol=5)
dev.off() 

#####T. Repens only in patches ####
ngs.nomb<-subset(as.data.frame(t(ngs.rarefied)),env$position!="mb"&env$position!="cor")
env.nomb<-subset(env,env$position!="mb"&env$position!="cor")
nomb.log<-log2(ngs.nomb+1)

ngs.bray.nm<-vegdist(nomb.log, method = "bray") # dissimilarity matrix using bray-curtis distance indices on the dataset native to vegan
pcoaVS.nm<-pco(ngs.bray.nm,negvals="zero", dround = 0) # if negvals = 0 sets all negative eigenvalues to zero; if = "rm" corrects for negative eigenvalues using method 1 of Legendre and Anderson 1999

library(ape)
Var<-pcoa(ngs.bray)
Var$values  #### The percentage of "variance"

plot(pcoaVS.nm$vectors[,2]~pcoaVS.nm$vectors[,1],xlab = "PCoA1", ylab="PCoA2",pch=19,col=env.nomb$code,
     axes = TRUE, main = "PCoA (ecodist) on White clover data")

pcoa.nm<-cbind(pcoaVS.nm$vectors[,2],pcoaVS.nm$vectors[,1],env.nomb)

####Figure 2B####
pdf("Figure. 2B PCoA of all T. Repens in patches 202111.pdf",height=6,width=8)
ggplot(data=pcoa.nm,aes(x=pcoaVS.nm$vectors[, 1],y=pcoaVS.nm$vectors[, 2]))+
  xlim(-0.05,0.05)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(code)),size=3)+ 
  stat_ellipse(geom="polygon",type="t",level=0.95,alpha=0.2,aes(fill=code))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens at individual level")
dev.off()

adonis.code<-adonis(nomb.log~code,env.nomb,by=NULL, method="bray") 
adonis.code

########Try this with random effect #####
adonis(nomb.log~code*corridor+(1|pool.code),env.nomb, permutations = 999, method = "bray")

adonis.cc<-adonis(nomb.log~code*corridor,env.nomb,by=NULL, method="bray") 
adonis.cc

##########

##### all T. repens including patches and corridors ####
ngs.allt<-subset(as.data.frame(t(ngs.rarefied)),env$position!="mb")
env.allt<-subset(env,env$position!="mb")
allt.log<-log2(ngs.allt+1)

ngs.bray.allt<-vegdist(allt.log, method = "bray") # dissimilarity matrix using bray-curtis distance indices on the dataset native to vegan
pcoaVS.allt<-pco(ngs.bray.allt,negvals="zero", dround = 0) # if negvals = 0 sets all negative eigenvalues to zero; if = "rm" corrects for negative eigenvalues using method 1 of Legendre and Anderson 1999
pcoa.allt<-cbind(pcoaVS.allt$vectors[,2],pcoaVS.allt$vectors[,1],env.allt)

pdf("PCoA of all T. repens samples p+c 20210831.pdf",height=5,width=7)
ggplot(data=pcoa.allt,aes(x=pcoaVS.allt$vectors[,1],y=pcoaVS.allt$vectors[,2]))+
  xlim(-0.052,0.05)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(position.code)),size=3)+ 
  stat_ellipse(geom="polygon",type="norm",level=0.95,alpha=0.2,aes(fill=position.code))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on all T. repens")
dev.off()

pcoa.nm1<-subset(pcoa.nm,pcoa.nm$corridor==1)
pcoa.nm0<-subset(pcoa.nm,pcoa.nm$corridor==0)

pdf("PCoA of T. repens samples with corridor 20210428.pdf",height=5,width=7)
ggplot(data=pcoa.nm1,aes(x=pcoa.nm1[,2],y=pcoa.nm1[,1]))+
  xlim(-0.06,0.06)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(code)),size=3)+ 
  stat_ellipse(geom="polygon",type="t",level=0.95,alpha=0.2,aes(fill=code))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens with corridor")
dev.off()

pdf("PCoA of T. repens samples without corridor 20210428.pdf",height=5,width=7)
ggplot(data=pcoa.nm0,aes(x=pcoa.nm0[,2],y=pcoa.nm0[,1]))+
  xlim(-0.06,0.06)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(code)),size=3)+ 
  stat_ellipse(geom="polygon",type="t",level=0.95,alpha=0.2,aes(fill=code))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens without corridor")
dev.off()

pdf("PCoA of all T. repens samples 20210428 patches.pdf",height=5,width=7)
ggplot(data=pcoa.nm,aes(x=pcoaVS.nm$vectors[,1],y=pcoaVS.nm$vectors[,2]))+
  xlim(-0.05,0.06)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(corridor)),size=3)+
  #stat_ellipse(geom="polygon",type="t",level=0.95,alpha=0.2,aes(fill=corridor))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens")
dev.off()

###### 2017 October
pcoa.nm.2017<-subset(pcoa.nm,pcoa.nm$code=="17O")
pcoa.nm1.2017<-subset(pcoa.nm1,pcoa.nm1$code=="17O")
pcoa.nm0.2017<-subset(pcoa.nm0,pcoa.nm0$code=="17O")
pcoa.allt.2017<-subset(pcoa.allt,pcoa.allt$code=="17O")

pdf("PCoA of 2017 all T. repens samples 20210831.pdf",height=5,width=7)
allt.17o<-ggplot(data=pcoa.allt.2017,aes(x=pcoa.allt.2017[,2],y=pcoa.allt.2017[,1]))+
  xlim(-0.06,0.06)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(position.code)),size=3)+
  stat_ellipse(geom="polygon",type="norm",level=0.95,alpha=0.2,aes(fill=position.code))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA on all T. repens 2017 october")
allt.17o
dev.off()

pdf("PCoA of 2017 T. repens patch samples 20210428.pdf",height=5,width=7)
ggplot(data=pcoa.nm.2017,aes(x=pcoa.nm.2017[,2],y=pcoa.nm.2017[,1]))+
  xlim(-0.06,0.06)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(position)),size=3)+
  stat_ellipse(geom="polygon",type="norm",level=0.95,alpha=0.2,aes(fill=position))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens")
dev.off()

pdf("PCoA of 2017 with corridor T. repens samples 20210428.pdf",height=5,width=7)
ggplot(data=pcoa.nm1.2017,aes(x=pcoa.nm1.2017[,2],y=pcoa.nm1.2017[,1]))+
  xlim(-0.06,0.06)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(position)),size=3)+
  stat_ellipse(geom="polygon",type="norm",level=0.95,alpha=0.2,aes(fill=position))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens")
dev.off()

pdf("PCoA of 2017 without corridor T. repens samples 20210428.pdf",height=5,width=7)
ggplot(data=pcoa.nm0.2017,aes(x=pcoa.nm0.2017[,2],y=pcoa.nm0.2017[,1]))+
  xlim(-0.06,0.06)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(position)),size=3)+
  stat_ellipse(geom="polygon",type="norm",level=0.95,alpha=0.2,aes(fill=position))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens")
dev.off()

#### significance test for community composition ####
ngs.2017<-subset(ngs.rarefied.t,env$code=="17O"&env$plant=="Trifolium"&env$position!="cor")
env.2017<-subset(env,env$code=="17O"&env$plant=="Trifolium"&env$position!="cor")
adonis.17o<-adonis(ngs.2017~corridor,env.2017,by=NULL, method="bray") 
adonis.17o

ngs.2018m<-subset(ngs.rarefied.t,env$code=="18 M"&env$plant=="Trifolium"&env$position!="cor")
env.2018m<-subset(env,env$code=="18 M"&env$plant=="Trifolium"&env$position!="cor")
adonis.18m<-adonis(ngs.2018m~corridor,env.2018m,by=NULL, method="bray") 
adonis.18m

pcoa.allt.2018m<-subset(pcoa.allt,pcoa.allt$code=="18 M")
pdf("PCoA of 2018 may all T. repens samples p+c 20210831.pdf",height=5,width=7)
allt.18m<-ggplot(data=pcoa.allt.2018m,aes(x=pcoa.allt.2018m[,2],y=pcoa.allt.2018m[,1]))+
  xlim(-0.06,0.06)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(position.code)),size=3)+
  stat_ellipse(geom="polygon",type="norm",level=0.95,alpha=0.2,aes(fill=position.code))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA on all T. repens 2018 May")
allt.18m
dev.off()

ngs.2018j<-subset(ngs.rarefied.t,env$code=="18J"&env$plant=="Trifolium"&env$position!="cor")
env.2018j<-subset(env,env$code=="18J"&env$plant=="Trifolium"&env$position!="cor")
adonis.18j<-adonis(ngs.2018j~corridor,env.2018j,by=NULL, method="bray") 
adonis.18j

pcoa.allt.2018j<-subset(pcoa.allt,pcoa.allt$code=="18J")
pdf("PCoA of 2018 june all T. repens samples p+c 20210831.pdf",height=5,width=7)
allt.18j<-ggplot(data=pcoa.allt.2018j,aes(x=pcoa.allt.2018j[,2],y=pcoa.allt.2018j[,1]))+
  xlim(-0.06,0.06)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(position.code)),size=3)+
  stat_ellipse(geom="polygon",type="norm",level=0.95,alpha=0.2,aes(fill=position.code))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA on all T. repens 2018 June")
allt.18j
dev.off()

ngs.2018o<-subset(ngs.rarefied.t,env$code=="18O"&env$plant=="Trifolium"&env$position!="cor")
env.2018o<-subset(env,env$code=="18O"&env$plant=="Trifolium"&env$position!="cor")
adonis.18o<-adonis(ngs.2018o~corridor,env.2018o,by=NULL, method="bray") 
adonis.18o

pcoa.allt.2018o<-subset(pcoa.allt,pcoa.allt$code=="18O")
pdf("PCoA of 2018 october all T. repens samples p+c 20210831.pdf",height=5,width=7)
allt.18o<-ggplot(data=pcoa.allt.2018o,aes(x=pcoa.allt.2018o[,2],y=pcoa.allt.2018o[,1]))+
  xlim(-0.06,0.06)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(position.code)),size=3)+
  stat_ellipse(geom="polygon",type="norm",level=0.95,alpha=0.2,aes(fill=position.code))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA on all T. repens 2018 October")
allt.18o
dev.off()

ngs.2019m<-subset(ngs.rarefied.t,env$code=="19M"&env$plant=="Trifolium"&env$position!="cor")
env.2019m<-subset(env,env$code=="19M"&env$plant=="Trifolium"&env$position!="cor")
adonis.19m<-adonis(ngs.2019m~corridor,env.2019m,by=NULL, method="bray") 
adonis.19m

###### 2019 May
pcoa.nm.2019<-subset(pcoa.nm,pcoa.nm$code=="19M")
pcoa.nm1.2019<-subset(pcoa.nm1,pcoa.nm1$code=="19M")
pcoa.nm0.2019<-subset(pcoa.nm0,pcoa.nm0$code=="19M")
pcoa.allt.2019<-subset(pcoa.allt,pcoa.allt$code=="19M")

pdf("PCoA of 2019 all T. repens samples p+c 20210831.pdf",height=5,width=7)
allt.19m<-ggplot(data=pcoa.allt.2019,aes(x=pcoa.allt.2019[,2],y=pcoa.allt.2019[,1]))+
  xlim(-0.06,0.06)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(position.code)),size=3)+
  stat_ellipse(geom="polygon",type="norm",level=0.95,alpha=0.2,aes(fill=position.code))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA on all T. repens 2019 May")
allt.19m
dev.off()

#####Figure S4####
pdf("Figure S4. all T. repens p+c.pdf",height=5,width=33)
grid.arrange(allt.17o,allt.18m,allt.18j,allt.18o,allt.19m,ncol=5)
dev.off()  
######

pdf("PCoA of 2019 T. repens samples 20210428.pdf",height=5,width=7)
ggplot(data=pcoa.nm.2019,aes(x=pcoa.nm.2019[,2],y=pcoa.nm.2019[,1]))+
  xlim(-0.06,0.06)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(position)),size=3)+
  stat_ellipse(geom="polygon",type="norm",level=0.95,alpha=0.2,aes(fill=position))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens")
dev.off()

pdf("PCoA of 2019 with corridor T. repens samples 20210428.pdf",height=5,width=7)
ggplot(data=pcoa.nm1.2019,aes(x=pcoa.nm1.2019[,2],y=pcoa.nm1.2019[,1]))+
  xlim(-0.08,0.06)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(position)),size=3)+
  stat_ellipse(geom="polygon",type="norm",level=0.95,alpha=0.2,aes(fill=position))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens")
dev.off()

pdf("PCoA of 2019 without corridor T. repens samples 20210428.pdf",height=5,width=7)
ggplot(data=pcoa.nm0.2019,aes(x=pcoa.nm0.2019[,2],y=pcoa.nm0.2019[,1]))+
  xlim(-0.06,0.04)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(position)),size=3)+
  stat_ellipse(geom="polygon",type="norm",level=0.95,alpha=0.2,aes(fill=position))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens")
dev.off()

#######################################################
### PCoA at population level for T. repens ####
#######################################################

ngs.poolt<-ngs.poolpp[,7:272]
env.poolt<-ngs.poolpp[,1:6]
ngs.poolt.log<-log2(ngs.poolt+1)

adonis.ccp<-adonis(ngs.poolt.log~code,env.poolt,by=NULL, method="bray") 
adonis.ccp

ngs.bray.nm<-vegdist(ngs.poolt.log, method = "bray") # dissimilarity matrix using bray-curtis distance indices on the dataset native to vegan
pcoaVS.nm<-pco(ngs.bray.nm,negvals="zero", dround = 0) # if negvals = 0 sets all negative eigenvalues to zero; if = "rm" corrects for negative eigenvalues using method 1 of Legendre and Anderson 1999

plot(pcoaVS.nm$vectors[,2]~pcoaVS.nm$vectors[,1],xlab = "PCoA1", ylab="PCoA2",pch=19,col=env.poolt$code,
     axes = TRUE, main = "PCoA (ecodist) on White clover data")

pcoa.nm.p<-cbind(pcoaVS.nm$vectors[,2],pcoaVS.nm$vectors[,1],env.poolt)
pcoa.nm.p1<-subset(pcoa.nm.p,pcoa.nm.p$corridor==1)
pcoa.nm.p0<-subset(pcoa.nm.p,pcoa.nm.p$corridor==0)

pdf("PCoA of all T. repens samples at population level 20210428 corridor.pdf",height=5,width=7)
ggplot(data=pcoa.nm.p,aes(x=pcoa.nm.p[,2],y=pcoa.nm.p[,1]))+
  xlim(-0.11,0.16)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(corridor)),size=3)+
  #stat_ellipse(geom="polygon",type="t",level=0.95,alpha=0.2,aes(fill=corridor))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens")
dev.off()

#####Figure 2C####
pdf("Figure 2C. PCoA of all T. repens samples at population level 20210428 code.pdf",height=5,width=7)
ggplot(data=pcoa.nm.p,aes(x=pcoa.nm.p[,2],y=pcoa.nm.p[,1]))+
  xlim(-0.14,0.16)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(code)),size=3)+
  stat_ellipse(geom="polygon",type="t",level=0.95,alpha=0.2,aes(fill=code))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens at population level")
dev.off()

#####
pdf("PCoA of T. repens samples with corridor at populaton level 20210428.pdf",height=5,width=7)
ggplot(data=pcoa.nm.p1,aes(x=pcoa.nm.p1[,2],y=pcoa.nm.p1[,1]))+
  xlim(-0.12,0.15)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(code)),size=3)+ 
  stat_ellipse(geom="polygon",type="t",level=0.95,alpha=0.2,aes(fill=code))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens with corridor at population level")
dev.off()

pdf("PCoA of T. repens samples without corridor at population level 20210428.pdf",height=5,width=7)
ggplot(data=pcoa.nm.p0,aes(x=pcoa.nm.p0[,2],y=pcoa.nm.p0[,1]))+
  xlim(-0.13,0.15)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(code)),size=3)+ 
  stat_ellipse(geom="polygon",type="t",level=0.95,alpha=0.2,aes(fill=code))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens without corridor at population level")
dev.off()


### 2017 October 
pcoa.nm.p2017<-subset(pcoa.nm.p,pcoa.nm.p$code=="2017october")
pcoa.nm1.p2017<-subset(pcoa.nm.p1,pcoa.nm.p1$code=="2017october")
pcoa.nm0.p2017<-subset(pcoa.nm.p0,pcoa.nm.p0$code=="2017october")

pdf("PCoA of 2017 T. repens samples at population level 20210428.pdf",height=5,width=7)
ggplot(data=pcoa.nm.p2017,aes(x=pcoa.nm.p2017[,2],y=pcoa.nm.p2017[,1]))+
  xlim(-0.13,0.11)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(corridor)),size=3)+
  #stat_ellipse(geom="polygon",type="norm",level=0.95,alpha=0.2,aes(fill=corridor))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens")
dev.off()

pdf("PCoA of 2017 with corridor T. repens samples at population level 20210428.pdf",height=5,width=7)
ggplot(data=pcoa.nm1.p2017,aes(x=pcoa.nm1.p2017[,2],y=pcoa.nm1.p2017[,1]))+
  xlim(-0.15,0.12)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(position)),size=3)+
  stat_ellipse(geom="polygon",type="norm",level=0.95,alpha=0.2,aes(fill=position))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens")
dev.off()

pdf("PCoA of 2017 without corridor T. repens samples at population level 20210428.pdf",height=5,width=7)
ggplot(data=pcoa.nm0.p2017,aes(x=pcoa.nm0.p2017[,2],y=pcoa.nm0.p2017[,1]))+
  xlim(-0.12,0.12)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(position)),size=3)+
  stat_ellipse(geom="polygon",type="norm",level=0.95,alpha=0.2,aes(fill=position))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens")
dev.off()

#### 2019 May 
pcoa.nm.p2019<-subset(pcoa.nm.p,pcoa.nm.p$code=="2019 may")
pcoa.nm1.p2019<-subset(pcoa.nm.p1,pcoa.nm.p1$code=="2019 may")
pcoa.nm0.p2019<-subset(pcoa.nm.p0,pcoa.nm.p0$code=="2019 may")

pdf("PCoA of 2019 T. repens samples at population level 20210428.pdf",height=5,width=7)
ggplot(data=pcoa.nm.p2019,aes(x=pcoa.nm.p2019[,2],y=pcoa.nm.p2019[,1]))+
  xlim(-0.08,0.14)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(corridor)),size=3)+
  #stat_ellipse(geom="polygon",type="t",level=0.95,alpha=0.2,aes(fill=corridor))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens")
dev.off()

pdf("PCoA of 2019 with corridor T. repens samples at population level 20210428.pdf",height=5,width=7)
ggplot(data=pcoa.nm1.p2019,aes(x=pcoa.nm1.p2019[,2],y=pcoa.nm1.p2019[,1]))+
  xlim(-0.16,0.2)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(position)),size=3)+
  stat_ellipse(geom="polygon",type="t",level=0.95,alpha=0.2,aes(fill=position))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens")
dev.off()

pdf("PCoA of 2019 without corridor T. repens samples at population level 20210428.pdf",height=5,width=7)
ggplot(data=pcoa.nm0.p2019,aes(x=pcoa.nm0.p2019[,2],y=pcoa.nm0.p2019[,1]))+
  xlim(-0.18,0.20)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(position)),size=3)+
  stat_ellipse(geom="polygon",type="norm",level=0.95,alpha=0.2,aes(fill=position))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens")
dev.off()

#### statistics at population level####
ngs.2017p<-subset(ngs.poolt,env.poolt$code=="2017october")
env.2017p<-subset(env.poolt,env.poolt$code=="2017october")
adonis.17op<-adonis(ngs.2017p~corridor,env.2017p,by=NULL, method="bray") 
adonis.17op

ngs.2018mp<-subset(ngs.poolt,env.poolt$code=="2018 may")
env.2018mp<-subset(env.poolt,env.poolt$code=="2018 may")
adonis.18mp<-adonis(ngs.2018mp~corridor,env.2018mp,by=NULL, method="bray") 
adonis.18mp

ngs.2018jp<-subset(ngs.poolt,env.poolt$code=="2018june")
env.2018jp<-subset(env.poolt,env.poolt$code=="2018june")
adonis.18jp<-adonis(ngs.2018jp~corridor,env.2018jp,by=NULL, method="bray") 
adonis.18jp

ngs.2018op<-subset(ngs.poolt,env.poolt$code=="2018october")
env.2018op<-subset(env.poolt,env.poolt$code=="2018october")
adonis.18op<-adonis(ngs.2018op~corridor,env.2018op,by=NULL, method="bray") 
adonis.18op

ngs.2019mp<-subset(ngs.poolt,env.poolt$code=="2019 may")
env.2019mp<-subset(env.poolt,env.poolt$code=="2019 may")
adonis.19mp<-adonis(ngs.2019mp~corridor,env.2019mp,by=NULL, method="bray") 
adonis.19mp

### PCoA at meta-population level for T. repens ####
ngs.mp<-ngs.ppool[,5:270]
env.mp<-ngs.ppool[,1:4]
ngs.mp.log<-log2(ngs.mp+1)

adonis.codemp<-adonis(ngs.mp.log~Group.2,env.mp,by=NULL, method="bray") 
adonis.codemp

adonis.ccmp<-adonis(ngs.mp.log~Group.2*Group.3,env.mp,by=NULL, method="bray") 
adonis.ccmp

ngs.bray.nm<-vegdist(ngs.mp.log, method = "bray") # dissimilarity matrix using bray-curtis distance indices on the dataset native to vegan
pcoaVS.nm<-pco(ngs.bray.nm,negvals="zero", dround = 0) # if negvals = 0 sets all negative eigenvalues to zero; if = "rm" corrects for negative eigenvalues using method 1 of Legendre and Anderson 1999

plot(pcoaVS.nm$vectors[,2]~pcoaVS.nm$vectors[,1],xlab = "PCoA1", ylab="PCoA2",pch=19,col=env.poolt$code,
     axes = TRUE, main = "PCoA (ecodist) on White clover data")

pcoa.nm.mp<-cbind(pcoaVS.nm$vectors[,2],pcoaVS.nm$vectors[,1],env.mp)
pcoa.nm.mp1<-subset(pcoa.nm.mp,pcoa.nm.mp$Group.3==1)
pcoa.nm.mp0<-subset(pcoa.nm.mp,pcoa.nm.mp$Group.3==0)

Var<-pcoa(ngs.bray.nm)
Var$values  #### The percentage of "variance"

pdf("PCoA of all T. repens samples at meta-population level 20210428 corridor.pdf",height=5,width=7)
ggplot(data=pcoa.nm.mp,aes(x=pcoa.nm.mp[,2],y=pcoa.nm.mp[,1]))+
  xlim(-0.40,0.25)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(Group.3)),size=3)+
  #stat_ellipse(geom="polygon",type="t",level=0.95,alpha=0.2,aes(fill=corridor))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens")
dev.off()

#####Figure 2D ####
pdf("Figure. 2D PCoA of all T. repens samples at meta-population level 20210503 code.pdf",height=5,width=7)
ggplot(data=pcoa.nm.mp,aes(x=pcoa.nm.mp[,2],y=pcoa.nm.mp[,1]))+
  xlim(-0.40,0.40)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(Group.2)),size=3)+
  stat_ellipse(geom="polygon",type="norm",level=0.95,alpha=0.2,aes(fill=Group.2))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens at meta-pupolation level")
dev.off()

pdf("PCoA of T. repens samples with corridor at meta-populaton level 20210428.pdf",height=5,width=7)
ggplot(data=pcoa.nm.mp1,aes(x=pcoa.nm.mp1[,2],y=pcoa.nm.mp1[,1]))+
  xlim(-0.40,0.35)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(Group.2)),size=3)+ 
  stat_ellipse(geom="polygon",type="norm",level=0.95,alpha=0.2,aes(fill=Group.2))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens with corridor at meta-pupolation level")
dev.off()

pdf("PCoA of T. repens samples without corridor at meta-population level 20210428.pdf",height=5,width=7)
ggplot(data=pcoa.nm.mp0,aes(x=pcoa.nm.mp0[,2],y=pcoa.nm.mp0[,1]))+
  xlim(-0.45,0.45)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(Group.2)),size=3)+ 
  stat_ellipse(geom="polygon",type="t",level=0.95,alpha=0.2,aes(fill=Group.2))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens without corridor at meta-pupolation level")
dev.off()

### 2017 October 
pcoa.nm.mp2017<-subset(pcoa.nm.mp,pcoa.nm.mp$Group.2=="17O")
pcoa.nm1.mp2017<-subset(pcoa.nm.mp1,pcoa.nm.mp1$Group.2=="17O")
pcoa.nm0.mp2017<-subset(pcoa.nm.mp0,pcoa.nm.mp0$Group.2=="17O")

pdf("PCoA of 2017 T. repens samples at meta-population level 20210428.pdf",height=5,width=7)
ggplot(data=pcoa.nm.mp2017,aes(x=pcoa.nm.mp2017[,2],y=pcoa.nm.mp2017[,1]))+
  xlim(-0.25,0.20)+
  geom_point(position=position_dodge(width=0.75),aes(color=as.factor(Group.3)),size=3)+
  #stat_ellipse(geom="polygon",type="t",level=0.95,alpha=0.2,aes(fill=Group.3))+
  xlab("PCoA 1") + ylab("PCoA 2")+ ggtitle("PCoA (ecodist) on T. repens at meta-population level")
dev.off()

####statistics at meta-population level####
ngs.2017mp<-subset(ngs.mp,env.mp$Group.2=="17O")
env.2017mp<-subset(env.mp,env.mp$Group.2=="17O")
adonis.17omp<-adonis(ngs.2017mp~Group.3,env.2017mp,by=NULL, method="bray") 
adonis.17omp

ngs.2018mmp<-subset(ngs.mp,env.mp$Group.2=="18 M")
env.2018mmp<-subset(env.mp,env.mp$Group.2=="18 M")
adonis.18mmp<-adonis(ngs.2018mmp~Group.3,env.2018mmp,by=NULL, method="bray") 
adonis.18mmp

ngs.2018jmp<-subset(ngs.mp,env.mp$Group.2=="18J")
env.2018jmp<-subset(env.mp,env.mp$Group.2=="18J")
adonis.18jmp<-adonis(ngs.2018jmp~Group.3,env.2018jmp,by=NULL, method="bray") 
adonis.18jmp

ngs.2018omp<-subset(ngs.mp,env.mp$Group.2=="18O")
env.2018omp<-subset(env.mp,env.mp$Group.2=="18O")
adonis.18omp<-adonis(ngs.2018omp~Group.3,env.2018omp,by=NULL, method="bray") 
adonis.18omp

ngs.2019mmp<-subset(ngs.mp,env.mp$Group.2=="19M")
env.2019mmp<-subset(env.mp,env.mp$Group.2=="19M")
adonis.19mmp<-adonis(ngs.2019mmp~Group.3,env.2019mmp,by=NULL, method="bray") 
adonis.19mmp

