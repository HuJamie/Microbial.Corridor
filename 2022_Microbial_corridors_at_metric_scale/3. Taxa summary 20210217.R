####taxonomic summary

####abundance inside each phylum ####
ngs.fungi<-aggregate(ngs.rarefied, list(taxo$blast_perc_identity), FUN = sum)    ## otu reads inside phylum,based on normalized matrix
rownames(ngs.fungi)<-ngs.fungi[,1]  
ngs.fungi.t<-t(ngs.fungi[,-1])

matrix.b<-cbind(richness,ngs.fungi.t)

library(ggplot2)
pdf("Abunadnce within each phylum at individual level.pdf",height=5,width=15)
a.ab<-ggplot(data=matrix.b,aes(x=code,y=Ascomycota))+geom_boxplot()+geom_point()+
  facet_wrap(~corridor,ncol=1)

b.ab<-ggplot(data=matrix.b,aes(x=code,y=Basidiomycota))+geom_boxplot()+geom_point()+
  facet_wrap(~corridor,ncol=1)

c.ab<-ggplot(data=matrix.b,aes(x=code,y=Chytridiomycota))+geom_boxplot()+geom_point()+
  facet_wrap(~corridor,ncol=1)

g.ab<-ggplot(data=matrix.b,aes(x=code,y=Glomeromycota))+geom_boxplot()+geom_point()+
  facet_wrap(~corridor,ncol=1)

z.ab<-ggplot(data=matrix.b,aes(x=code,y=Zygomycota))+geom_boxplot()+geom_point()+
  facet_wrap(~corridor,ncol=1)

grid.arrange(a.ab,b.ab,c.ab,g.ab,z.ab,ncol=5)
dev.off()

## abundance in phylum at population level ####
matrix.tri<-subset(matrix.b,matrix.b$position!="mb"&matrix.b$position!="cor")
phylum.abun.pp<-aggregate(matrix.tri[,22:27], list(matrix.tri$pool.code,matrix.tri$code,matrix.tri$corridor), FUN=mean)    ## pool data

library(ggplot2)
pdf("Abunadnce within each phylum at population level.pdf",height=5,width=15)
a.abp<-ggplot(data=phylum.abun.pp,aes(x=Group.2,y=Ascomycota))+geom_boxplot()+geom_point()+
  facet_wrap(~Group.3,ncol=1)

b.abp<-ggplot(data=phylum.abun.pp,aes(x=Group.2,y=Basidiomycota))+geom_boxplot()+geom_point()+
  facet_wrap(~Group.3,ncol=1)

c.abp<-ggplot(data=phylum.abun.pp,aes(x=Group.2,y=Chytridiomycota))+geom_boxplot()+geom_point()+
  facet_wrap(~Group.3,ncol=1)

g.abp<-ggplot(data=phylum.abun.pp,aes(x=Group.2,y=Glomeromycota))+geom_boxplot()+geom_point()+
  facet_wrap(~Group.3,ncol=1)

z.abp<-ggplot(data=phylum.abun.pp,aes(x=Group.2,y=Zygomycota))+geom_boxplot()+geom_point()+
  facet_wrap(~Group.3,ncol=1)

grid.arrange(a.abp,b.abp,c.abp,g.abp,z.abp,ncol=5)
dev.off()


## abundance in phylum at meta-population level ####
phylum.abun.mp<-aggregate(matrix.tri[,22:27], list(matrix.tri$mesocosm,matrix.tri$code,matrix.tri$corridor), FUN=mean)    ## pool data

pdf("Abunadnce within each phylum at meta-population level.pdf",height=5,width=15)
a.abmp<-ggplot(data=phylum.abun.mp,aes(x=Group.2,y=Ascomycota))+geom_boxplot()+geom_point()+
  facet_wrap(~Group.3,ncol=1)

b.abmp<-ggplot(data=phylum.abun.mp,aes(x=Group.2,y=Basidiomycota))+geom_boxplot()+geom_point()+
  facet_wrap(~Group.3,ncol=1)

c.abmp<-ggplot(data=phylum.abun.mp,aes(x=Group.2,y=Chytridiomycota))+geom_boxplot()+geom_point()+
  facet_wrap(~Group.3,ncol=1)

g.abmp<-ggplot(data=phylum.abun.mp,aes(x=Group.2,y=Glomeromycota))+geom_boxplot()+geom_point()+
  facet_wrap(~Group.3,ncol=1)

z.abmp<-ggplot(data=phylum.abun.mp,aes(x=Group.2,y=Zygomycota))+geom_boxplot()+geom_point()+
  facet_wrap(~Group.3,ncol=1)

grid.arrange(a.abmp,b.abmp,c.abmp,g.abp,z.abmp,ncol=5)
dev.off()


##### calculate diversity index inside each phylum at individual level ####
ngs.phylum<-cbind(taxo$blast_perc_identity,ngs.rarefied)

ngs.Asco<-subset(ngs.phylum,ngs.phylum$`taxo$blast_perc_identity`=="Ascomycota")
ngs.Asco.r<-ngs.Asco[,-1]
library(vegan)
otu.shannon.asco<- diversity(t(ngs.Asco.r),index="shannon") 
otu.richness.asco<-rowSums(t(ngs.Asco.r)>0)
otu.evenness.asco<- otu.shannon.asco/log(otu.richness.asco)  ## Pielou's evenness

########Basidiomycota
ngs.Basidio<-subset(ngs.phylum,ngs.phylum$`taxo$blast_perc_identity`=="Basidiomycota")
ngs.Basidio.r<-ngs.Basidio[,-1]
otu.shannon.Basidio<- diversity(t(ngs.Basidio.r),index="shannon") 
otu.richness.Basidio<-rowSums(t(ngs.Basidio.r)>0)
otu.evenness.Basidio<- otu.shannon.Basidio/log(otu.richness.Basidio)  ## Pielou's evenness

##### Chytridiomycota
ngs.Chytridio<-subset(ngs.phylum,ngs.phylum$`taxo$blast_perc_identity`=="Chytridiomycota")
ngs.Chytridio.r<-ngs.Chytridio[,-1]
otu.shannon.Chytridio<- diversity(t(ngs.Chytridio.r),index="shannon") 
otu.richness.Chytridio<-rowSums(t(ngs.Chytridio.r)>0)
otu.evenness.Chytridio<- otu.shannon.Chytridio/log(otu.richness.Chytridio)  ## Pielou's evenness

######Glomeromycota
ngs.Glomero<-subset(ngs.phylum,ngs.phylum$`taxo$blast_perc_identity`=="Glomeromycota")
ngs.Glomero.r<-ngs.Glomero[,-1]
otu.shannon.Glomero<- diversity(t(ngs.Glomero.r),index="shannon") 
otu.richness.Glomero<-rowSums(t(ngs.Glomero.r)>0)
otu.evenness.Glomero<- otu.shannon.Glomero/log(otu.richness.Glomero)  ## Pielou's evenness

######Zygomycota
ngs.Zygo<-subset(ngs.phylum,ngs.phylum$`taxo$blast_perc_identity`=="Zygomycota")
ngs.Zygo.r<-ngs.Zygo[,-1]
otu.shannon.Zygo<- diversity(t(ngs.Zygo.r),index="shannon") 
otu.richness.Zygo<-rowSums(t(ngs.Zygo.r)>0)
otu.evenness.Zygo<- otu.shannon.Zygo/log(otu.richness.Zygo)  ## Pielou's evenness


richness.with.phylum<-cbind(env,otu.richness.asco,otu.shannon.asco,otu.evenness.asco,
                            otu.richness.Basidio,otu.shannon.Basidio,otu.evenness.Basidio,
                            otu.richness.Chytridio,otu.shannon.Chytridio,otu.evenness.Chytridio,
                            otu.richness.Glomero,otu.shannon.Glomero,otu.evenness.Glomero,
                            otu.richness.Zygo,otu.shannon.Zygo,otu.evenness.Zygo)

richness.nomb<-subset(richness,richness$plant!="Brachypodium")
richness.with.phylum.nomb<-subset(richness.with.phylum,richness.with.phylum$plant!="Brachypodium"&richness.with.phylum$position!="cor")

library(ggplot2)
library(gridExtra)
pdf("Richness in phylum at individual level.pdf",height=6,width=15)
a<-ggplot(data=richness.with.phylum.nomb,aes(x=code,y=otu.richness.asco))+geom_jitter(width=0.2)+geom_boxplot()+
  facet_wrap(~corridor,ncol=1)

b<-ggplot(data=richness.with.phylum.nomb,aes(x=code,y=otu.richness.Basidio))+geom_jitter(width=0.2)+geom_boxplot()+
  facet_wrap(~corridor,ncol=1)

c<-ggplot(data=richness.with.phylum.nomb,aes(x=code,y=otu.richness.Chytridio))+geom_jitter(width=0.2)+geom_boxplot()+
  facet_wrap(~corridor,ncol=1)

g<-ggplot(data=richness.with.phylum.nomb,aes(x=code,y=otu.richness.Glomero))+geom_jitter(width=0.2)+geom_boxplot()+
  facet_wrap(~corridor,ncol=1)

z<-ggplot(data=richness.with.phylum.nomb,aes(x=code,y=otu.richness.Zygo))+geom_jitter(width=0.2)+geom_boxplot()+
  facet_wrap(~corridor,ncol=1)

grid.arrange(a,b,c,g,z,ncol=5)
dev.off()

pdf("Evenness in phylum at individual level.pdf",height=6,width=15)
ae<-ggplot(data=richness.with.phylum,aes(x=code,y=otu.evenness.asco))+geom_jitter(width=0.2)+geom_boxplot()+
  facet_wrap(~corridor,ncol=1)

be<-ggplot(data=richness.with.phylum.nomb,aes(x=code,y=otu.evenness.Basidio))+geom_jitter(width=0.2)+geom_boxplot()+
  facet_wrap(~corridor,ncol=1)

ce<-ggplot(data=richness.with.phylum.nomb,aes(x=code,y=otu.evenness.Chytridio))+geom_jitter(width=0.2)+geom_boxplot()+
  facet_wrap(~corridor,ncol=1)

ge<-ggplot(data=richness.with.phylum.nomb,aes(x=code,y=otu.evenness.Glomero))+geom_jitter(width=0.2)+geom_boxplot()+
  facet_wrap(~corridor,ncol=1)

ze<-ggplot(data=richness.with.phylum.nomb,aes(x=code,y=otu.evenness.Zygo))+geom_jitter(width=0.2)+geom_boxplot()+
  facet_wrap(~corridor,ncol=1)

grid.arrange(ae,be,ce,ge,ze,ncol=5)
dev.off()

### Statistics for diversity inside each phylum at individual level ####
rich.with.phy.indip1<-subset(richness.with.phylum.nomb,richness.with.phylum.nomb$corridor==1)
rich.with.phy.indip0<-subset(richness.with.phylum.nomb,richness.with.phylum.nomb$corridor==0)
rich.2017.phy.indip<-subset(richness.with.phylum.nomb,richness.with.phylum.nomb$code=="17O")
rich.2018.may.phy.indip<-subset(richness.with.phylum.nomb,richness.with.phylum.nomb$code=="18 M")
rich.2018.june.phy.indip<-subset(richness.with.phylum.nomb,richness.with.phylum.nomb$code=="18J")
rich.2018.oct.phy.indip<-subset(richness.with.phylum.nomb,richness.with.phylum.nomb$code=="18O")
rich.2019.phy.indip<-subset(richness.with.phylum.nomb,richness.with.phylum.nomb$code=="19M")

####Corridor effect for Ascomycota####
lm.asc<-lmer(otu.richness.asco~corridor+(1|pool.code),rich.2017.phy.indip)
Anova(lm.asc)
AIC(lm.asc)
r.squaredGLMM(lm.asc)

lm.asc<-lmer(otu.evenness.asco~corridor+(1|pool.code),rich.2017.phy.indip)
Anova(lm.asc)
AIC(lm.asc)
r.squaredGLMM(lm.asc)

lm.asc<-lmer(otu.richness.asco~corridor+(1|pool.code),rich.2018.may.phy.indip)
Anova(lm.asc)
AIC(lm.asc)
r.squaredGLMM(lm.asc)

lm.asc<-lmer(otu.evenness.asco~corridor+(1|pool.code),rich.2018.may.phy.indip)
Anova(lm.asc)
AIC(lm.asc)
r.squaredGLMM(lm.asc)

lm.asc<-lmer(otu.richness.asco~corridor+(1|pool.code),rich.2018.june.phy.indip)
Anova(lm.asc)
AIC(lm.asc)
r.squaredGLMM(lm.asc)

lm.asc<-lmer(otu.evenness.asco~corridor+(1|pool.code),rich.2018.june.phy.indip)
Anova(lm.asc)
AIC(lm.asc)
r.squaredGLMM(lm.asc)

lm.asc<-lmer(otu.richness.asco~corridor+(1|pool.code),rich.2018.oct.phy.indip)
Anova(lm.asc)
AIC(lm.asc)
r.squaredGLMM(lm.asc)

lm.asc<-lmer(otu.evenness.asco~corridor+(1|pool.code),rich.2018.oct.phy.indip)
Anova(lm.asc)
AIC(lm.asc)
r.squaredGLMM(lm.asc)

lm.asc<-lmer(otu.richness.asco~corridor+(1|pool.code),rich.2019.phy.indip)
Anova(lm.asc)
AIC(lm.asc)
r.squaredGLMM(lm.asc)

lm.asc<-lmer(otu.evenness.asco~corridor+(1|pool.code),rich.2019.phy.indip)
Anova(lm.asc)
AIC(lm.asc)
r.squaredGLMM(lm.asc)

#### Corridor effect for Basidiomycota ####
lm.bas<-lmer(otu.richness.Basidio~corridor+(1|pool.code),rich.2017.phy.indip)
Anova(lm.bas)
AIC(lm.bas)
r.squaredGLMM(lm.bas)

lm.bas<-lmer(otu.evenness.Basidio~corridor+(1|pool.code),rich.2017.phy.indip)
Anova(lm.bas)
AIC(lm.bas)
r.squaredGLMM(lm.bas)

lm.bas<-lmer(otu.richness.Basidio~corridor+(1|pool.code),rich.2018.may.phy.indip)
Anova(lm.bas)
AIC(lm.bas)
r.squaredGLMM(lm.bas)

lm.bas<-lm(otu.evenness.Basidio~corridor+(1|pool.code),rich.2018.may.phy.indip)
Anova(lm.bas)
AIC(lm.bas)
r.squaredGLMM(lm.bas)

lm.bas<-lmer(otu.richness.Basidio~corridor+(1|pool.code),rich.2018.june.phy.indip)
Anova(lm.bas)
AIC(lm.bas)
r.squaredGLMM(lm.bas)

lm.bas<-lmer(otu.evenness.Basidio~corridor+(1|pool.code),rich.2018.june.phy.indip)
Anova(lm.bas)
AIC(lm.bas)
r.squaredGLMM(lm.bas)

lm.bas<-lmer(otu.richness.Basidio~corridor+(1|pool.code),rich.2018.oct.phy.indip)
Anova(lm.bas)
AIC(lm.bas)
r.squaredGLMM(lm.bas)

lm.bas<-lmer(otu.evenness.Basidio~corridor+(1|pool.code),rich.2018.oct.phy.indip)
Anova(lm.bas)
AIC(lm.bas)
r.squaredGLMM(lm.bas)

lm.bas<-lmer(otu.richness.Basidio~corridor+(1|pool.code),rich.2019.phy.indip)
Anova(lm.bas)
AIC(lm.bas)
r.squaredGLMM(lm.bas)

lm.bas<-lmer(otu.evenness.Basidio~corridor+(1|pool.code),rich.2019.phy.indip)
Anova(lm.bas)
AIC(lm.bas)
r.squaredGLMM(lm.bas)

#### Corridor effect for Chytridiomycota####
lm.chy<-lmer(otu.richness.Chytridio~corridor+(1|pool.code),rich.2017.phy.indip)
Anova(lm.chy)
AIC(lm.chy)
r.squaredGLMM(lm.chy)

lm.chy<-lmer(otu.evenness.Chytridio~corridor+(1|pool.code),rich.2017.phy.indip)
Anova(lm.chy)
AIC(lm.chy)
r.squaredGLMM(lm.chy)

lm.chy<-lmer(otu.richness.Chytridio~corridor+(1|pool.code),rich.2018.may.phy.indip)
Anova(lm.chy)
AIC(lm.chy)
r.squaredGLMM(lm.chy)

lm.chy<-lmer(otu.evenness.Chytridio~corridor+(1|pool.code),rich.2018.may.phy.indip)
Anova(lm.chy)
AIC(lm.chy)
r.squaredGLMM(lm.chy)

lm.chy<-lmer(otu.richness.Chytridio~corridor+(1|pool.code),rich.2018.june.phy.indip)
Anova(lm.chy)
AIC(lm.chy)
r.squaredGLMM(lm.chy)

lm.chy<-lmer(otu.evenness.Chytridio~corridor+(1|pool.code),rich.2018.june.phy.indip)
Anova(lm.chy)
AIC(lm.chy)
r.squaredGLMM(lm.chy)

lm.chy<-lmer(otu.richness.Chytridio~corridor+(1|pool.code),rich.2018.oct.phy.indip)
Anova(lm.chy)
AIC(lm.chy)
r.squaredGLMM(lm.chy)

lm.chy<-lmer(otu.evenness.Chytridio~corridor+(1|pool.code),rich.2018.oct.phy.indip)
Anova(lm.chy)
AIC(lm.chy)
r.squaredGLMM(lm.chy)

lm.chy<-lmer(otu.richness.Chytridio~corridor+(1|pool.code),rich.2019.phy.indip)
Anova(lm.chy)
AIC(lm.chy)
r.squaredGLMM(lm.chy)

lm.chy<-lmer(otu.evenness.Chytridio~corridor+(1|pool.code),rich.2019.phy.indip)
Anova(lm.chy)
AIC(lm.chy)
r.squaredGLMM(lm.chy)

#### Corridor effect for Glomeromycota

lm.glomero<-lmer(otu.richness.Glomero~corridor+(1|pool.code),rich.2017.phy.indip)
Anova(lm.glomero)
AIC(lm.glomero)
r.squaredGLMM(lm.glomero)

lm.glomero<-lmer(otu.evenness.Glomero~corridor+(1|pool.code),rich.2017.phy.indip)
Anova(lm.glomero)
AIC(lm.glomero)
r.squaredGLMM(lm.glomero)

lm.glomero<-lmer(otu.richness.Glomero~corridor+(1|pool.code),rich.2018.may.phy.indip)
Anova(lm.glomero)
AIC(lm.glomero)
r.squaredGLMM(lm.glomero)
lm.glomero<-lmer(otu.evenness.Glomero~corridor+(1|pool.code),rich.2018.may.phy.indip)
Anova(lm.glomero)
AIC(lm.glomero)
r.squaredGLMM(lm.glomero)

lm.glomero<-lmer(otu.richness.Glomero~corridor+(1|pool.code),rich.2018.june.phy.indip)
Anova(lm.glomero)
AIC(lm.glomero)
r.squaredGLMM(lm.glomero)
lm.glomero<-lmer(otu.evenness.Glomero~corridor+(1|pool.code),rich.2018.june.phy.indip)
Anova(lm.glomero)
AIC(lm.glomero)
r.squaredGLMM(lm.glomero)

lm.glomero<-lmer(otu.richness.Glomero~corridor+(1|pool.code),rich.2018.oct.phy.indip)
Anova(lm.glomero)
AIC(lm.glomero)
r.squaredGLMM(lm.glomero)
lm.glomero<-lmer(otu.evenness.Glomero~corridor+(1|pool.code),rich.2018.oct.phy.indip)
Anova(lm.glomero)
AIC(lm.glomero)
r.squaredGLMM(lm.glomero)

lm.glomero<-lmer(otu.richness.Glomero~corridor+(1|pool.code),rich.2019.phy.indip)
Anova(lm.glomero)
AIC(lm.glomero)
r.squaredGLMM(lm.glomero)

lm.glomero<-lmer(otu.evenness.Glomero~corridor+(1|pool.code),rich.2019.phy.indip)
Anova(lm.glomero)
AIC(lm.glomero)
r.squaredGLMM(lm.glomero)

#### Corridor effect on Zygomycota


phy.div.stat.in<-cbind(adr.i,br.i,cr.i,gr.i,zr.i)
write.table(phy.div.stat.in,file="phy.div.stat.in.txt",sep="\t")

phy.eve.stat.in<-cbind(ade.i,be.i,ce.i,ge.i,ze.i)
write.table(phy.eve.stat.in,file="phy.eve.stat.in.txt",sep="\t")

lm.zygo<-lmer(otu.richness.Zygo~corridor+(1|pool.code),rich.2017.phy.indip)
Anova(lm.zygo)
AIC(lm.zygo)
r.squaredGLMM(lm.zygo)
lm.zygo<-lmer(otu.evenness.Zygo~corridor+(1|pool.code),rich.2017.phy.indip)
Anova(lm.zygo)
AIC(lm.zygo)
r.squaredGLMM(lm.zygo)

lm.zygo<-lmer(otu.richness.Zygo~corridor+(1|pool.code),rich.2018.may.phy.indip)
Anova(lm.zygo)
AIC(lm.zygo)
r.squaredGLMM(lm.zygo)
lm.zygo<-lmer(otu.evenness.Zygo~corridor+(1|pool.code),rich.2018.may.phy.indip)
Anova(lm.zygo)
AIC(lm.zygo)
r.squaredGLMM(lm.zygo)

lm.zygo<-lmer(otu.richness.Zygo~corridor+(1|pool.code),rich.2018.june.phy.indip)
Anova(lm.zygo)
AIC(lm.zygo)
r.squaredGLMM(lm.zygo)
lm.zygo<-lmer(otu.evenness.Zygo~corridor+(1|pool.code),rich.2018.june.phy.indip)
Anova(lm.zygo)
AIC(lm.zygo)
r.squaredGLMM(lm.zygo)

lm.zygo<-lmer(otu.richness.Zygo~corridor+(1|pool.code),rich.2018.oct.phy.indip)
Anova(lm.zygo)
AIC(lm.zygo)
r.squaredGLMM(lm.zygo)
lm.zygo<-lmer(otu.evenness.Zygo~corridor+(1|pool.code),rich.2018.oct.phy.indip)
Anova(lm.zygo)
AIC(lm.zygo)
r.squaredGLMM(lm.zygo)

lm.zygo<-lmer(otu.richness.Zygo~corridor+(1|pool.code),rich.2019.phy.indip)
Anova(lm.zygo)
AIC(lm.zygo)
r.squaredGLMM(lm.zygo)
lm.zygo<-lmer(otu.evenness.Zygo~corridor+(1|pool.code),rich.2019.phy.indip)
Anova(lm.zygo)
AIC(lm.zygo)
r.squaredGLMM(lm.zygo)

#### calculate diversity index inside each phylum at population level ####

env.pool<-ngs.pool[,1:7]

ngs.poolp<-subset(ngs.pool,ngs.pool$Group.5!="cor"&ngs.pool$Group.5!="mb")
ngs.poolp.t<-as.data.frame(t(ngs.poolp[,8:273]))
ngs.poolp.phylum<-cbind(taxo$blast_perc_identity,ngs.poolp.t)

ngs.pool.t<-as.data.frame(t(ngs.pool[,8:273]))
ngs.pool.phylum<-cbind(taxo$blast_perc_identity,ngs.pool.t)

ngs.pool.Asco<-subset(ngs.pool.phylum,ngs.pool.phylum$`taxo$blast_perc_identity`=="Ascomycota")
ngs.pool.Asco.r<-ngs.pool.Asco[,-1]
library(vegan)
otu.shannon.asco.p<- diversity(t(ngs.pool.Asco.r),index="shannon") 
otu.richness.asco.p<-rowSums(t(ngs.pool.Asco.r)>0)
otu.evenness.asco.p<- otu.shannon.asco.p/log(otu.richness.asco.p)  ## Pielou's evenness

########Basidiomycota
ngs.pool.Basidio<-subset(ngs.pool.phylum,ngs.pool.phylum$`taxo$blast_perc_identity`=="Basidiomycota")
ngs.pool.Basidio.r<-ngs.pool.Basidio[,-1]
otu.shannon.Basidio.p<- diversity(t(ngs.pool.Basidio.r),index="shannon") 
otu.richness.Basidio.p<-rowSums(t(ngs.pool.Basidio.r)>0)
otu.evenness.Basidio.p<- otu.shannon.Basidio.p/log(otu.richness.Basidio.p)  ## Pielou's evenness

##### Chytridiomycota
ngs.pool.Chytridio<-subset(ngs.pool.phylum,ngs.pool.phylum$`taxo$blast_perc_identity`=="Chytridiomycota")
ngs.pool.Chytridio.r<-ngs.pool.Chytridio[,-1]
otu.shannon.Chytridio.p<- diversity(t(ngs.pool.Chytridio.r),index="shannon") 
otu.richness.Chytridio.p<-rowSums(t(ngs.pool.Chytridio.r)>0)
otu.evenness.Chytridio.p<- otu.shannon.Chytridio.p/log(otu.richness.Chytridio.p)  ## Pielou's evenness

######Glomeromycota
ngs.pool.Glomero<-subset(ngs.pool.phylum,ngs.pool.phylum$`taxo$blast_perc_identity`=="Glomeromycota")
ngs.pool.Glomero.r<-ngs.pool.Glomero[,-1]
otu.shannon.Glomero.p<- diversity(t(ngs.pool.Glomero.r),index="shannon") 
otu.richness.Glomero.p<-rowSums(t(ngs.pool.Glomero.r)>0)
otu.evenness.Glomero.p<- otu.shannon.Glomero.p/log(otu.richness.Glomero.p)  ## Pielou's evenness

######Zygomycota 
ngs.pool.Zygo<-subset(ngs.pool.phylum,ngs.pool.phylum$`taxo$blast_perc_identity`=="Zygomycota")
ngs.pool.Zygo.r<-ngs.pool.Zygo[,-1]
otu.shannon.Zygo.p<- diversity(t(ngs.pool.Zygo.r),index="shannon") 
otu.richness.Zygo.p<-rowSums(t(ngs.pool.Zygo.r)>0)
otu.evenness.Zygo.p<- otu.shannon.Zygo.p/log(otu.richness.Zygo.p)  ## Pielou's evenness

richness.pool.phylum<-cbind(env.pool,otu.richness.asco.p,otu.shannon.asco.p,otu.evenness.asco.p,
                            otu.richness.Basidio.p,otu.shannon.Basidio.p,otu.evenness.Basidio.p,
                            otu.richness.Chytridio.p,otu.shannon.Chytridio.p,otu.evenness.Chytridio.p,
                            otu.richness.Glomero.p,otu.shannon.Glomero.p,otu.evenness.Glomero.p,
                            otu.richness.Zygo.p,otu.shannon.Zygo.p,otu.evenness.Zygo.p)

richness.pool.phylum.tri<-subset(richness.pool.phylum,richness.pool.phylum$Group.5!="mb"&richness.pool.phylum$Group.5!="cor")
richness.pool.phylum.mb<-subset(richness.pool.phylum,richness.pool.phylum$Group.5=="mb")

library(ggplot2)
library(gridExtra)
pdf("Richness in phylum tri species at population level.pdf",height=6,width=15)
a.p<-ggplot(data=richness.pool.phylum.tri,aes(x=Group.4,y=otu.richness.asco.p))+geom_boxplot()+geom_jitter(width=0.2)+
  ylim(0,55)+
  facet_wrap(~Group.6,ncol=1)

b.p<-ggplot(data=richness.pool.phylum.tri,aes(x=Group.4,y=otu.richness.Basidio.p))+geom_boxplot()+geom_jitter(width=0.2)+
  ylim(0,55)+
  facet_wrap(~Group.6,ncol=1)

c.p<-ggplot(data=richness.pool.phylum.tri,aes(x=Group.4,y=otu.richness.Chytridio.p))+geom_boxplot()+geom_jitter(width=0.2)+
  ylim(0,55)+
  facet_wrap(~Group.6,ncol=1)

g.p<-ggplot(data=richness.pool.phylum.tri,aes(x=Group.4,y=otu.richness.Glomero.p))+geom_boxplot()+geom_jitter(width=0.2)+
  ylim(0,55)+
  facet_wrap(~Group.6,ncol=1)

z.p<-ggplot(data=richness.pool.phylum.tri,aes(x=Group.4,y=otu.richness.Zygo.p))+geom_boxplot()+geom_jitter(width=0.2)+
  ylim(0,55)+
  facet_wrap(~Group.6,ncol=1)

grid.arrange(a.p,b.p,c.p,g.p,z.p,ncol=5)
dev.off()

pdf("Richness in phylum pool mb species.pdf",height=6,width=15)
amb.p<-ggplot(data=richness.pool.phylum.mb,aes(x=Group.4,y=otu.richness.asco.p))+geom_boxplot()+geom_jitter(width=0.2)+
  ylim(0,55)+
  facet_wrap(~Group.6,ncol=1)

bmb.p<-ggplot(data=richness.pool.phylum.mb,aes(x=Group.4,y=otu.richness.asco.p))+geom_boxplot()+geom_jitter(width=0.2)+
  ylim(0,55)+
  facet_wrap(~Group.6,ncol=1)

cmb.p<-ggplot(data=richness.pool.phylum.mb,aes(x=Group.4,y=otu.richness.Chytridio.p))+geom_boxplot()+geom_jitter(width=0.2)+
  ylim(0,55)+
  facet_wrap(~Group.6,ncol=1)

gmb.p<-ggplot(data=richness.pool.phylum.mb,aes(x=Group.4,y=otu.richness.Glomero.p))+geom_boxplot()+geom_jitter(width=0.2)+
  ylim(0,55)+
  facet_wrap(~Group.6,ncol=1)

zmb.p<-ggplot(data=richness.pool.phylum.mb,aes(x=Group.4,y=otu.richness.Zygo.p))+geom_boxplot()+geom_jitter(width=0.2)+
  ylim(0,55)+
  facet_wrap(~Group.6,ncol=1)

grid.arrange(amb.p,bmb.p,cmb.p,gmb.p,zmb.p,ncol=5)
dev.off()

#### statistics for diversity  inside each phylum at population level ####
richness.pool.phylum.tri

rich.2017.phy.p<-subset(richness.pool.phylum.tri,richness.pool.phylum.tri$Group.4=="17O")
rich.2018.may.phy.p<-subset(richness.pool.phylum.tri,richness.pool.phylum.tri$Group.4=="18 M")
rich.2018.june.phy.p<-subset(richness.pool.phylum.tri,richness.pool.phylum.tri$Group.4=="18J")
rich.2018.oct.phy.p<-subset(richness.pool.phylum.tri,richness.pool.phylum.tri$Group.4=="18O")
rich.2019.phy.p<-subset(richness.pool.phylum.tri,richness.pool.phylum.tri$Group.4=="19M")

#### Corridor effect for Ascomycota#### 
lm.asc<-lmer(otu.richness.asco.p~Group.6+(1|Group.7),rich.2017.phy.p)
Anova(lm.asc)
AIC(lm.asc)
r.squaredGLMM(lm.asc)

lm.asc<-lmer(otu.evenness.asco.p~Group.6+(1|Group.7),rich.2017.phy.p)
Anova(lm.asc)
AIC(lm.asc)
r.squaredGLMM(lm.asc)

lm.asc<-lmer(otu.richness.asco.p~Group.6+(1|Group.7),rich.2018.may.phy.p)
Anova(lm.asc)
AIC(lm.asc)
r.squaredGLMM(lm.asc)
lm.asc<-lmer(otu.evenness.asco.p~Group.6+(1|Group.7),rich.2018.may.phy.p)
Anova(lm.asc)
AIC(lm.asc)
r.squaredGLMM(lm.asc)

lm.asc<-lmer(otu.richness.asco.p~Group.6+(1|Group.7),rich.2018.june.phy.p)
Anova(lm.asc)
AIC(lm.asc)
r.squaredGLMM(lm.asc)
lm.asc<-lmer(otu.evenness.asco.p~Group.6+(1|Group.7),rich.2018.june.phy.p)
Anova(lm.asc)
AIC(lm.asc)
r.squaredGLMM(lm.asc)

lm.asc<-lmer(otu.richness.asco.p~Group.6+(1|Group.7),rich.2018.oct.phy.p)
Anova(lm.asc)
AIC(lm.asc)
r.squaredGLMM(lm.asc)
lm.asc<-lmer(otu.evenness.asco.p~Group.6+(1|Group.7),rich.2018.oct.phy.p)
Anova(lm.asc)
AIC(lm.asc)
r.squaredGLMM(lm.asc)

lm.asc<-lmer(otu.richness.asco.p~Group.6+(1|Group.7),rich.2019.phy.p)
Anova(lm.asc)
AIC(lm.asc)
r.squaredGLMM(lm.asc)
lm.asc<-lmer(otu.evenness.asco.p~Group.6+(1|Group.7),rich.2019.phy.p)
Anova(lm.asc)
AIC(lm.asc)
r.squaredGLMM(lm.asc)

#### Corridor effect for Basidiomycota ####
lm.bas<-lmer(otu.richness.Basidio.p~Group.6+(1|Group.7),rich.2017.phy.p)
Anova(lm.bas)
AIC(lm.bas)
r.squaredGLMM(lm.bas)
lm.bas<-lmer(otu.evenness.Basidio.p~Group.6+(1|Group.7),rich.2017.phy.p)
Anova(lm.bas)
AIC(lm.bas)
r.squaredGLMM(lm.bas)

lm.bas<-lmer(otu.richness.Basidio.p~Group.6+(1|Group.7),rich.2018.may.phy.p)
Anova(lm.bas)
AIC(lm.bas)
r.squaredGLMM(lm.bas)
lm.bas<-lmer(otu.evenness.Basidio.p~Group.6+(1|Group.7),rich.2018.may.phy.p)
Anova(lm.bas)
AIC(lm.bas)
r.squaredGLMM(lm.bas)

lm.bas<-lmer(otu.richness.Basidio.p~Group.6+(1|Group.7),rich.2018.june.phy.p)
Anova(lm.bas)
AIC(lm.bas)
r.squaredGLMM(lm.bas)
lm.bas<-lmer(otu.evenness.Basidio.p~Group.6+(1|Group.7),rich.2018.june.phy.p)
Anova(lm.bas)
AIC(lm.bas)
r.squaredGLMM(lm.bas)

lm.bas<-lmer(otu.richness.Basidio.p~Group.6+(1|Group.7),rich.2018.oct.phy.p)
Anova(lm.bas)
AIC(lm.bas)
r.squaredGLMM(lm.bas)
lm.bas<-lmer(otu.evenness.Basidio.p~Group.6+(1|Group.7),rich.2018.oct.phy.p)
Anova(lm.bas)
AIC(lm.bas)
r.squaredGLMM(lm.bas)

lm.bas<-lmer(otu.richness.Basidio.p~Group.6+(1|Group.7),rich.2019.phy.p)
Anova(lm.bas)
AIC(lm.bas)
r.squaredGLMM(lm.bas)
lm.bas<-lmer(otu.evenness.Basidio.p~Group.6+(1|Group.7),rich.2019.phy.p)
Anova(lm.bas)
AIC(lm.bas)
r.squaredGLMM(lm.bas)

#### Corridor effect for Chytridiomycota
lm.chy<-lmer(otu.richness.Chytridio.p~Group.6+(1|Group.7),rich.2017.phy.p)
Anova(lm.chy)
AIC(lm.chy)
r.squaredGLMM(lm.chy)
lm.chy<-lmer(otu.evenness.Chytridio.p~Group.6+(1|Group.7),rich.2017.phy.p)
Anova(lm.chy)
AIC(lm.chy)
r.squaredGLMM(lm.chy)

lm.chy<-lmer(otu.richness.Chytridio.p~Group.6+(1|Group.7),rich.2018.may.phy.p)
Anova(lm.chy)
AIC(lm.chy)
r.squaredGLMM(lm.chy)
lm.chy<-lmer(otu.evenness.Chytridio.p~Group.6+(1|Group.7),rich.2018.may.phy.p)
Anova(lm.chy)
AIC(lm.chy)
r.squaredGLMM(lm.chy)

lm.chy<-lmer(otu.richness.Chytridio.p~Group.6+(1|Group.7),rich.2018.june.phy.p)
Anova(lm.chy)
AIC(lm.chy)
r.squaredGLMM(lm.chy)
lm.chy<-lmer(otu.evenness.Chytridio.p~Group.6+(1|Group.7),rich.2018.june.phy.p)
Anova(lm.chy)
AIC(lm.chy)
r.squaredGLMM(lm.chy)

lm.chy<-lmer(otu.richness.Chytridio.p~Group.6+(1|Group.7),rich.2018.oct.phy.p)
Anova(lm.chy)
AIC(lm.chy)
r.squaredGLMM(lm.chy)
lm.chy<-lmer(otu.evenness.Chytridio.p~Group.6+(1|Group.7),rich.2018.oct.phy.p)
Anova(lm.chy)
AIC(lm.chy)
r.squaredGLMM(lm.chy)

lm.chy<-lmer(otu.richness.Chytridio.p~Group.6+(1|Group.7),rich.2019.phy.p)
Anova(lm.chy)
AIC(lm.chy)
r.squaredGLMM(lm.chy)
lm.chy<-lmer(otu.evenness.Chytridio.p~Group.6+(1|Group.7),rich.2019.phy.p)
Anova(lm.chy)
AIC(lm.chy)
r.squaredGLMM(lm.chy)

#### Corridor effect for Glomeromycota
lm.glo<-lmer(otu.richness.Glomero.p~Group.6+(1|Group.7),rich.2017.phy.p)
Anova(lm.glo)
AIC(lm.glo)
r.squaredGLMM(lm.glo)
lm.glo<-lmer(otu.evenness.Glomero.p~Group.6+(1|Group.7),rich.2017.phy.p)
Anova(lm.glo)
AIC(lm.glo)
r.squaredGLMM(lm.glo)

lm.glo<-lmer(otu.richness.Glomero.p~Group.6+(1|Group.7),rich.2018.may.phy.p)
Anova(lm.glo)
AIC(lm.glo)
r.squaredGLMM(lm.glo)
lm.glo<-lmer(otu.evenness.Glomero.p~Group.6+(1|Group.7),rich.2018.may.phy.p)
Anova(lm.glo)
AIC(lm.glo)
r.squaredGLMM(lm.glo)

lm.glo<-lmer(otu.richness.Glomero.p~Group.6+(1|Group.7),rich.2018.june.phy.p)
Anova(lm.glo)
AIC(lm.glo)
r.squaredGLMM(lm.glo)
lm.glo<-lmer(otu.evenness.Glomero.p~Group.6+(1|Group.7),rich.2018.june.phy.p)
Anova(lm.glo)
AIC(lm.glo)
r.squaredGLMM(lm.glo)

lm.glo<-lmer(otu.richness.Glomero.p~Group.6+(1|Group.7),rich.2018.oct.phy.p)
Anova(lm.glo)
AIC(lm.glo)
r.squaredGLMM(lm.glo)
lm.glo<-lmer(otu.evenness.Glomero.p~Group.6+(1|Group.7),rich.2018.oct.phy.p)
Anova(lm.glo)
AIC(lm.glo)
r.squaredGLMM(lm.glo)

lm.glo<-lmer(otu.richness.Glomero.p~Group.6+(1|Group.7),rich.2019.phy.p)
Anova(lm.glo)
AIC(lm.glo)
r.squaredGLMM(lm.glo)
lm.glo<-lmer(otu.evenness.Glomero.p~Group.6+(1|Group.7),rich.2019.phy.p)
Anova(lm.glo)
AIC(lm.glo)
r.squaredGLMM(lm.glo)

#### Corridor effect on Zygomycota
lm.zyg<-lm(otu.richness.Zygo.p~Group.4*Group.6,richness.pool.phylum.tri)
lm.zyg.xp<-step(lm.zyg, trace=F, direction="both")
zr.p<-anova(lm.zyg)
summary(lm.zyg)
AIC(lm.zyg)

lm.zyg<-lm(otu.evenness.Zygo.p~Group.4*Group.6,richness.pool.phylum.tri)
lm.zyg.xp<-step(lm.zyg, trace=F, direction="both")
ze.p<-anova(lm.zyg)
summary(lm.zyg)
AIC(lm.zyg)

phy.div.stat.p<-cbind(ar.p,br.p,cr.p,gr.p,zr.p)
write.table(phy.div.stat.p,file="phy.div.stat.p.txt",sep="\t")

phy.eve.stat.p<-cbind(ae.p,be.p,ce.p,ge.p,ze.p)
write.table(phy.eve.stat.p,file="phy.eve.stat.p.txt",sep="\t")


lm.zyg<-lmer(otu.richness.Zygo.p~Group.6+(1|Group.7),rich.2017.phy.p)
Anova(lm.zyg)
AIC(lm.zyg)
r.squaredGLMM(lm.zyg)
lm.zyg<-lmer(otu.evenness.Zygo.p~Group.6+(1|Group.7),rich.2017.phy.p)
Anova(lm.zyg)
AIC(lm.zyg)
r.squaredGLMM(lm.zyg)

lm.zyg<-lmer(otu.richness.Zygo.p~Group.6+(1|Group.7),rich.2018.may.phy.p)
Anova(lm.zyg)
AIC(lm.zyg)
r.squaredGLMM(lm.zyg)
lm.zyg<-lmer(otu.evenness.Zygo.p~Group.6+(1|Group.7),rich.2018.may.phy.p)
Anova(lm.zyg)
AIC(lm.zyg)
r.squaredGLMM(lm.zyg)

lm.zyg<-lmer(otu.richness.Zygo.p~Group.6+(1|Group.7),rich.2018.june.phy.p)
Anova(lm.zyg)
AIC(lm.zyg)
r.squaredGLMM(lm.zyg)
lm.zyg<-lmer(otu.evenness.Zygo.p~Group.6+(1|Group.7),rich.2018.june.phy.p)
Anova(lm.zyg)
AIC(lm.zyg)
r.squaredGLMM(lm.zyg)

lm.zyg<-lmer(otu.richness.Zygo.p~Group.6+(1|Group.7),rich.2018.oct.phy.p)
Anova(lm.zyg)
AIC(lm.zyg)
r.squaredGLMM(lm.zyg)
lm.zyg<-lmer(otu.evenness.Zygo.p~Group.6+(1|Group.7),rich.2018.oct.phy.p)
Anova(lm.zyg)
AIC(lm.zyg)
r.squaredGLMM(lm.zyg)

lm.zyg<-lmer(otu.richness.Zygo.p~Group.6+(1|Group.7),rich.2019.phy.p)
Anova(lm.zyg)
AIC(lm.zyg)
r.squaredGLMM(lm.zyg)
lm.zyg<-lmer(otu.evenness.Zygo.p~Group.6+(1|Group.7),rich.2019.phy.p)
Anova(lm.zyg)
AIC(lm.zyg)
r.squaredGLMM(lm.zyg)


#### calculate diversity index inside each phylum at meta-population level ####
ngs.pp.p<-aggregate(ngs.te[,18:283], list(ngs.te$mesocosm,ngs.te$code,ngs.te$corridor), FUN=mean)    ## pool data
env.pp<-ngs.pp.p[1:3]
ngs.pp.t<-as.data.frame(t(ngs.ppool[,4:269]))

ngs.pp.phylum<-cbind(taxo$blast_perc_identity,ngs.pp.t)

ngs.pp.Asco<-subset(ngs.pp.phylum,ngs.pp.phylum$`taxo$blast_perc_identity`=="Ascomycota")
ngs.pp.Asco.r<-ngs.pp.Asco[,-1]
library(vegan)
otu.shannon.asco.pp<- diversity(t(ngs.pp.Asco.r),index="shannon") 
otu.richness.asco.pp<-rowSums(t(ngs.pp.Asco.r)>0)
otu.evenness.asco.pp<- otu.shannon.asco.pp/log(otu.richness.asco.pp)  ## Pielou's evenness

########Basidiomycota
ngs.pp.Basidio<-subset(ngs.pp.phylum,ngs.pp.phylum$`taxo$blast_perc_identity`=="Basidiomycota")
ngs.pp.Basidio.r<-ngs.pp.Basidio[,-1]
otu.shannon.Basidio.pp<- diversity(t(ngs.pp.Basidio.r),index="shannon") 
otu.richness.Basidio.pp<-rowSums(t(ngs.pp.Basidio.r)>0)
otu.evenness.Basidio.pp<- otu.shannon.Basidio.pp/log(otu.richness.Basidio.pp)  ## Pielou's evenness

##### Chytridiomycota
ngs.pp.Chytridio<-subset(ngs.pp.phylum,ngs.pp.phylum$`taxo$blast_perc_identity`=="Chytridiomycota")
ngs.pp.Chytridio.r<-ngs.pp.Chytridio[,-1]
otu.shannon.Chytridio.pp<- diversity(t(ngs.pp.Chytridio.r),index="shannon") 
otu.richness.Chytridio.pp<-rowSums(t(ngs.pp.Chytridio.r)>0)
otu.evenness.Chytridio.pp<- otu.shannon.Chytridio.pp/log(otu.richness.Chytridio.pp)  ## Pielou's evenness

######Glomeromycota
ngs.pp.Glomero<-subset(ngs.pp.phylum,ngs.pp.phylum$`taxo$blast_perc_identity`=="Glomeromycota")
ngs.pp.Glomero.r<-ngs.pp.Glomero[,-1]
otu.shannon.Glomero.pp<- diversity(t(ngs.pp.Glomero.r),index="shannon") 
otu.richness.Glomero.pp<-rowSums(t(ngs.pp.Glomero.r)>0)
otu.evenness.Glomero.pp<- otu.shannon.Glomero.pp/log(otu.richness.Glomero.pp)  ## Pielou's evenness

######Zygomycota
ngs.pp.Zygo<-subset(ngs.pp.phylum,ngs.pp.phylum$`taxo$blast_perc_identity`=="Zygomycota")
ngs.pp.Zygo.r<-ngs.pp.Zygo[,-1]
otu.shannon.Zygo.pp<- diversity(t(ngs.pp.Zygo.r),index="shannon") 
otu.richness.Zygo.pp<-rowSums(t(ngs.pp.Zygo.r)>0)
otu.evenness.Zygo.pp<- otu.shannon.Zygo.pp/log(otu.richness.Zygo.pp)  ## Pielou's evenness

richness.pp.phylum<-cbind(env.pp,otu.richness.asco.pp,otu.shannon.asco.pp,otu.evenness.asco.pp,
                            otu.richness.Basidio.pp,otu.shannon.Basidio.pp,otu.evenness.Basidio.pp,
                            otu.richness.Chytridio.pp,otu.shannon.Chytridio.pp,otu.evenness.Chytridio.pp,
                            otu.richness.Glomero.pp,otu.shannon.Glomero.pp,otu.evenness.Glomero.pp,
                            otu.richness.Zygo.pp,otu.shannon.Zygo.pp,otu.evenness.Zygo.pp)

library(ggplot2)
library(gridExtra)
pdf("Richness in phylum tri species at meta population level.pdf",height=6,width=15)
a.pp<-ggplot(data=richness.pp.phylum,aes(x=Group.2,y=otu.richness.asco.pp))+geom_boxplot()+geom_jitter(width=0.2)+
  ylim(0,62)+
  facet_wrap(~Group.3,ncol=1)

b.pp<-ggplot(data=richness.pp.phylum,aes(x=Group.2,y=otu.richness.Basidio.pp))+geom_boxplot()+geom_jitter(width=0.2)+
  ylim(0,62)+
  facet_wrap(~Group.3,ncol=1)

c.pp<-ggplot(data=richness.pp.phylum,aes(x=Group.2,y=otu.richness.Chytridio.pp))+geom_boxplot()+geom_jitter(width=0.2)+
  ylim(0,62)+
  facet_wrap(~Group.3,ncol=1)

g.pp<-ggplot(data=richness.pp.phylum,aes(x=Group.2,y=otu.richness.Glomero.pp))+geom_boxplot()+geom_jitter(width=0.2)+
  ylim(0,62)+
  facet_wrap(~Group.3,ncol=1)

z.pp<-ggplot(data=richness.pp.phylum,aes(x=Group.2,y=otu.richness.Zygo.pp))+geom_boxplot()+geom_jitter(width=0.2)+
  ylim(0,62)+
  facet_wrap(~Group.3,ncol=1)

grid.arrange(a.pp,b.pp,c.pp,g.pp,z.pp,ncol=5)
dev.off()

pdf("Richness in phylum pool mb species.pdf",height=6,width=15)
amb.p<-ggplot(data=richness.pool.phylum.mb,aes(x=Group.4,y=otu.richness.asco.p))+geom_jitter(width=0.2)+geom_boxplot()+
  facet_wrap(~Group.6,ncol=1)

bmb.p<-ggplot(data=richness.pool.phylum.mb,aes(x=Group.4,y=otu.richness.asco.p))+geom_jitter(width=0.2)+geom_boxplot()+
  facet_wrap(~Group.6,ncol=1)

cmb.p<-ggplot(data=richness.pool.phylum.mb,aes(x=Group.4,y=otu.richness.Chytridio.p))+geom_jitter(width=0.2)+geom_boxplot()+
  facet_wrap(~Group.6,ncol=1)

gmb.p<-ggplot(data=richness.pool.phylum.mb,aes(x=Group.4,y=otu.richness.Glomero.p))+geom_jitter(width=0.2)+geom_boxplot()+
  facet_wrap(~Group.6,ncol=1)

zmb.p<-ggplot(data=richness.pool.phylum.mb,aes(x=Group.4,y=otu.richness.Zygo.p))+geom_jitter(width=0.2)+geom_boxplot()+
  facet_wrap(~Group.6,ncol=1)
grid.arrange(amb.p,bmb.p,cmb.p,gmb.p,zmb.p,ncol=5)
dev.off()


#### statistics for diversity  inside each phylum at meta-population level ####
richness.pp.phylum

rich.2017.phy.mp<-subset(richness.pp.phylum,richness.pp.phylum$Group.2=="17O")
rich.2018.may.phy.mp<-subset(richness.pp.phylum,richness.pp.phylum$Group.2=="18 M")
rich.2018.june.phy.mp<-subset(richness.pp.phylum,richness.pp.phylum$Group.2=="18J")
rich.2018.oct.phy.mp<-subset(richness.pp.phylum,richness.pp.phylum$Group.2=="18O")
rich.2019.phy.mp<-subset(richness.pp.phylum,richness.pp.phylum$Group.2=="19M")

#### #### Corridor effect for Ascomycota
lm.asco<-glm(otu.richness.asco.pp~Group.3,rich.2017.phy.mp,family=poisson(link = "log"))
summary(lm.asco)
Anova(lm.asco)
rsq(lm.asco)
hist(lm.asco)
lm.asco<-lm(otu.evenness.asco.pp~Group.3,rich.2017.phy.mp)
summary(lm.asco)
Anova(lm.asco)
rsq(lm.asco)
hist(lm.asco)

lm.asco<-lm(otu.richness.asco.pp~Group.3,rich.2018.may.phy.mp)
summary(lm.asco)
Anova(lm.asco)
rsq(lm.asco)
hist(lm.asco)
lm.asco<-glm(otu.evenness.asco.pp~Group.3,rich.2018.may.phy.mp,family=gaussian(link = "identity"))
summary(lm.asco)
Anova(lm.asco)
rsq(lm.asco)
AIC(lm.asco)

lm.asco<-lm(otu.richness.asco.pp~Group.3,rich.2018.june.phy.mp)
summary(lm.asco)
Anova(lm.asco)
rsq(lm.asco)
AIC(lm.asco)
lm.asco<-lm(otu.evenness.asco.pp~Group.3,rich.2018.june.phy.mp)
summary(lm.asco)
Anova(lm.asco)
rsq(lm.asco)
AIC(lm.asco)

lm.asco<-glm(otu.richness.asco.pp~Group.3,rich.2018.oct.phy.mp,family=gaussian(link = "identity"))
summary(lm.asco)
Anova(lm.asco)
rsq(lm.asco)
AIC(lm.asco)
lm.asco<-glm(otu.evenness.asco.pp~Group.3,rich.2018.oct.phy.mp,family=gaussian(link = "identity"))
summary(lm.asco)
Anova(lm.asco)
rsq(lm.asco)
AIC(lm.asco)

lm.asco<-glm(otu.richness.asco.pp~Group.3,rich.2019.phy.mp,family=gaussian(link = "identity"))
summary(lm.asco)
Anova(lm.asco)
rsq(lm.asco)
AIC(lm.asco)
lm.asco<-glm(otu.evenness.asco.pp~Group.3,rich.2019.phy.mp,family=gaussian(link = "identity"))
summary(lm.asco)
Anova(lm.asco)
rsq(lm.asco)
AIC(lm.asco)

#### Corridor effect for Basidiomycota
lm.Basidio<-glm(otu.richness.Basidio.pp~Group.3,rich.2017.phy.mp,family=gaussian(link = "identity"))
summary(lm.Basidio)
Anova(lm.Basidio)
rsq(lm.Basidio)
AIC(lm.Basidio)
lm.Basidio<-lm(otu.evenness.Basidio.pp~Group.3,rich.2017.phy.mp)
summary(lm.Basidio)
Anova(lm.Basidio)
rsq(lm.Basidio)
AIC(lm.Basidio)

lm.Basidio<-lm(otu.richness.Basidio.pp~Group.3,rich.2018.may.phy.mp)
summary(lm.Basidio)
Anova(lm.Basidio)
rsq(lm.Basidio)
AIC(lm.Basidio)
lm.Basidio<-lm(otu.evenness.Basidio.pp~Group.3,rich.2018.may.phy.mp)
summary(lm.Basidio)
Anova(lm.Basidio)
rsq(lm.Basidio)
AIC(lm.Basidio)

lm.Basidio<-lm(otu.richness.Basidio.pp~Group.3,rich.2018.june.phy.mp)
summary(lm.Basidio)
Anova(lm.Basidio)
rsq(lm.Basidio)
AIC(lm.Basidio)
lm.Basidio<-lm(otu.evenness.Basidio.pp~Group.3,rich.2018.june.phy.mp)
summary(lm.Basidio)
Anova(lm.Basidio)
rsq(lm.Basidio)
AIC(lm.Basidio)

lm.Basidio<-glm(otu.richness.Basidio.pp~Group.3,rich.2018.oct.phy.mp,family=gaussian(link = "identity"))
summary(lm.Basidio)
Anova(lm.Basidio)
rsq(lm.Basidio)
AIC(lm.Basidio)
#####Significant results
lm.Basidio<-glm(otu.evenness.Basidio.pp~Group.3,rich.2018.oct.phy.mp,family=gaussian(link = "identity"))
summary(lm.Basidio)
Anova(lm.Basidio)
rsq(lm.Basidio)
AIC(lm.Basidio)

lm.Basidio<-lm(otu.richness.Basidio.pp~Group.3,rich.2019.phy.mp)
summary(lm.Basidio)
Anova(lm.Basidio)
rsq(lm.Basidio)
AIC(lm.Basidio)
lm.Basidio<-lm(otu.evenness.Basidio.pp~Group.3,rich.2019.phy.mp)
summary(lm.Basidio)
Anova(lm.Basidio)
rsq(lm.Basidio)
AIC(lm.Basidio)

#### Corridor effect for Chytridiomycota
lm.Chytridio<-lm(otu.richness.Chytridio.pp~Group.3,rich.2017.phy.mp)
summary(lm.Chytridio)
Anova(lm.Chytridio)
rsq(lm.Chytridio)
AIC(lm.Chytridio)
lm.Chytridio<-lm(otu.evenness.Chytridio.pp~Group.3,rich.2017.phy.mp)
summary(lm.Chytridio)
Anova(lm.Chytridio)
rsq(lm.Chytridio)
AIC(lm.Chytridio)

lm.Chytridio<-lm(otu.richness.Chytridio.pp~Group.3,rich.2018.may.phy.mp)
summary(lm.Chytridio)
Anova(lm.Chytridio)
rsq(lm.Chytridio)
AIC(lm.Chytridio)
lm.Chytridio<-lm(otu.evenness.Chytridio.pp~Group.3,rich.2018.may.phy.mp)
summary(lm.Chytridio)
Anova(lm.Chytridio)
rsq(lm.Chytridio)
AIC(lm.Chytridio)

lm.Chytridio<-lm(otu.richness.Chytridio.pp~Group.3,rich.2018.june.phy.mp)
summary(lm.Chytridio)
Anova(lm.Chytridio)
rsq(lm.Chytridio)
AIC(lm.Chytridio)
lm.Chytridio<-lm(otu.evenness.Chytridio.pp~Group.3,rich.2018.june.phy.mp)
summary(lm.Chytridio)
Anova(lm.Chytridio)
rsq(lm.Chytridio)
AIC(lm.Chytridio)

lm.Chytridio<-glm(otu.richness.Chytridio.pp~Group.3,rich.2018.oct.phy.mp,family=gaussian(link = "identity"))
summary(lm.Chytridio)
Anova(lm.Chytridio)
rsq(lm.Chytridio)
AIC(lm.Chytridio)
lm.Chytridio<-lm(otu.evenness.Chytridio.pp~Group.3,rich.2018.oct.phy.mp)
summary(lm.Chytridio)
Anova(lm.Chytridio)
rsq(lm.Chytridio)
AIC(lm.Chytridio)

lm.Chytridio<-lm(otu.richness.Chytridio.pp~Group.3,rich.2019.phy.mp)
summary(lm.Chytridio)
Anova(lm.Chytridio)
rsq(lm.Chytridio)
AIC(lm.Chytridio)
lm.Chytridio<-lm(otu.evenness.Chytridio.pp~Group.3,rich.2019.phy.mp)
summary(lm.Chytridio)
Anova(lm.Chytridio)
rsq(lm.Chytridio)
AIC(lm.Chytridio)

#### Corridor effect for Glomeromycota
lm.Glomero<-lm(otu.richness.Glomero.pp~Group.3,rich.2017.phy.mp)
summary(lm.Glomero)
Anova(lm.Glomero)
rsq(lm.Glomero)
AIC(lm.Glomero)
lm.Glomero<-lm(otu.evenness.Glomero.pp~Group.3,rich.2017.phy.mp)
summary(lm.Glomero)
Anova(lm.Glomero)
rsq(lm.Glomero)
AIC(lm.Glomero)

lm.Glomero<-lm(otu.richness.Glomero.pp~Group.3,rich.2018.may.phy.mp)
summary(lm.Glomero)
Anova(lm.Glomero)
rsq(lm.Glomero)
AIC(lm.Glomero)
lm.Glomero<-lm(otu.evenness.Glomero.pp~Group.3,rich.2018.may.phy.mp)
summary(lm.Glomero)
Anova(lm.Glomero)
rsq(lm.Glomero)
AIC(lm.Glomero)

lm.Glomero<-lm(otu.richness.Glomero.pp~Group.3,rich.2018.june.phy.mp)
summary(lm.Glomero)
Anova(lm.Glomero)
rsq(lm.Glomero)
AIC(lm.Glomero)
lm.Glomero<-lm(otu.evenness.Glomero.pp~Group.3,rich.2018.june.phy.mp)
summary(lm.Glomero)
Anova(lm.Glomero)
rsq(lm.Glomero)
AIC(lm.Glomero)

lm.Glomero<-glm(otu.richness.Glomero.pp~Group.3,rich.2018.oct.phy.mp,family=gaussian(link = "identity"))
summary(lm.Glomero)
Anova(lm.Glomero)
rsq(lm.Glomero)
AIC(lm.Glomero)
lm.Glomero<-glm(otu.evenness.Glomero.pp~Group.3,rich.2018.oct.phy.mp,family=gaussian(link = "identity"))
summary(lm.Glomero)
Anova(lm.Glomero)
rsq(lm.Glomero)
AIC(lm.Glomero)

lm.Glomero<-lm(otu.richness.Glomero.pp~Group.3,rich.2019.phy.mp)
summary(lm.Glomero)
Anova(lm.Glomero)
rsq(lm.Glomero)
AIC(lm.Glomero)
lm.Glomero<-lm(otu.evenness.Glomero.pp~Group.3,rich.2019.phy.mp)
summary(lm.Glomero)
Anova(lm.Glomero)
rsq(lm.Glomero)
AIC(lm.Glomero)

#### Corridor effect on Zygomycota

phy.div.stat.mp<-cbind(ar.mp,br.mp,cr.mp,gr.mp,zr.mp)
write.table(phy.div.stat.mp,file="phy.div.stat.mp.txt",sep="\t")

phy.eve.stat.mp<-cbind(ae.mp,be.mp,ce.mp,ge.mp,ze.mp)
write.table(phy.eve.stat.mp,file="phy.eve.stat.mp.txt",sep="\t")


lm.zygo<-lm(otu.richness.Zygo.pp~Group.3,rich.2017.phy.mp)
summary(lm.zygo)
Anova(lm.zygo)
rsq(lm.zygo)
AIC(lm.zygo)
lm.zygo<-lm(otu.evenness.Zygo.pp~Group.3,rich.2017.phy.mp)
summary(lm.zygo)
Anova(lm.zygo)
rsq(lm.zygo)
AIC(lm.zygo)

lm.zygo<-lm(otu.richness.Zygo.pp~Group.3,rich.2018.may.phy.mp)
summary(lm.zygo)
Anova(lm.zygo)
rsq(lm.zygo)
AIC(lm.zygo)
lm.zygo<-lm(otu.evenness.Zygo.pp~Group.3,rich.2018.may.phy.mp)
summary(lm.zygo)
Anova(lm.zygo)
rsq(lm.zygo)
AIC(lm.zygo)

lm.zygo<-lm(otu.richness.Zygo.pp~Group.3,rich.2018.june.phy.mp)
summary(lm.zygo)
Anova(lm.zygo)
rsq(lm.zygo)
AIC(lm.zygo)
lm.zygo<-lm(otu.evenness.Zygo.pp~Group.3,rich.2018.june.phy.mp)
summary(lm.zygo)
Anova(lm.zygo)
rsq(lm.zygo)
AIC(lm.zygo)

lm.zygo<-glm(otu.richness.Zygo.pp~Group.3,rich.2018.oct.phy.mp,family=gaussian(link = "identity"))
summary(lm.zygo)
Anova(lm.zygo)
rsq(lm.zygo)
AIC(lm.zygo)
lm.zygo<-lm(otu.evenness.Zygo.pp~Group.3,rich.2018.oct.phy.mp)
summary(lm.zygo)
Anova(lm.zygo)
rsq(lm.zygo)
AIC(lm.zygo)

lm.zygo<-lm(otu.richness.Zygo.pp~Group.3,rich.2019.phy.mp)
summary(lm.zygo)
Anova(lm.zygo)
rsq(lm.zygo)
AIC(lm.zygo)
lm.zygo<-lm(otu.evenness.Zygo.pp~Group.3,rich.2019.phy.mp)
summary(lm.zygo)
Anova(lm.zygo)
rsq(lm.zygo)
AIC(lm.zygo)

