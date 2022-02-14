###Statistical test
library(ggplot2)
library(ggprism)
library(gridExtra)
# 1 -----------------------------------------------------------------------
ONE_Tukey_HSD1 <- function(data,group,compare,value){
  
  library(multcomp)#Tukey检验需要用到这个包来标显著性字母标记
  
  a <- data.frame(stringsAsFactors = F)#做一个空的数据框
  type <- unique(data[,group])#统计需要运行多重比较的次数
  for (i in type)#进行type次多重比较
  {
    g1=compare
    sub_dat <- data[data[,group]==i,]#根据指定的i去取相应的数据集出来
    #fit <- aov(sub_dat[,value] ~ sub_dat[,compare] )
    names(sub_dat)[names(sub_dat)==compare] <- 'g1' ## 重命名方便后面使用
    names(sub_dat)[names(sub_dat)==value] <- 'value' ## 重命名方便后面使用
    sub_dat$g1 <- factor(sub_dat$g1)#将列转化成因子以进行多重比较
    
    fit <- aov(value ~ g1,data = sub_dat )#方差分析
    #Tukey_HSD = TukeyHSD(fit, ordered = TRUE, conf.level = 0.95)
    options(warn = -1)
    tuk <- cld(glht(fit, alternative = 'two.sided', linfct = mcp(g1 = 'Tukey')), decreasing = TRUE)#Tukey检验多重比较
    Tukey.labels <- data.frame(Letters=tuk$mcletters$Letters, stringsAsFactors = FALSE)#获取多重比较字母标注
    Tukey.labels$compare = rownames(Tukey.labels)## 提取字母分组行名为group组名
    Tukey.labels$type <- i
    
    mean_sd <- merge(aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=sd),#获取数据标准差
                     aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=mean),by="Group.1"#获取数据均值
    )
    names(mean_sd) <- c('compare','std','mean')#列名重命名
    
    a <- rbind(a,merge(mean_sd,Tukey.labels,by='compare'))#合并数据
  }
  
  names(a) <- c(compare,'std','mean','Letters',group)#列名重命名
  return(a)
}

# 2 -----------------------------------------------------------------------

ONE_Tukey_HSD2 <- function(data,group,compare,value){
  library(multcompView)
  
  a <- data.frame(stringsAsFactors = F)
  type <- unique(data[,group])
  for (i in type)
  {
    g1=compare
    sub_dat <- data[data[,group]==i,]
    #fit <- aov(sub_dat[,value] ~ sub_dat[,compare] )
    ## 重命名方便后面使用
    names(sub_dat)[names(sub_dat)==compare] <- 'g1'
    names(sub_dat)[names(sub_dat)==value] <- 'value'
    sub_dat$g1 <- factor(sub_dat$g1)
    
    fit <- aov(value ~ g1,data = sub_dat )
    Tukey_HSD = TukeyHSD(fit, ordered = TRUE, conf.level = 0.95)
    options(warn = -1)
    tuk <- multcompLetters2(value ~ g1, Tukey_HSD$g1[,"p adj"], sub_dat)
    
    
    #tuk <- cld(glht(fit, alternative = 'two.sided', linfct = mcp(g1 = 'Tukey')), decreasing = TRUE)
    Tukey.labels <- data.frame(tuk['Letters'], stringsAsFactors = FALSE)
    ## 提取字母分组行名为group组名
    Tukey.labels$compare = rownames(Tukey.labels)
    Tukey.labels$type <- i
    
    mean_sd <- merge(aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=sd),
                     aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=mean),by="Group.1"
    )
    names(mean_sd) <- c('compare','std','mean')
    
    a <- rbind(a,merge(mean_sd,Tukey.labels,by='compare'))
  }
  
  names(a) <- c(compare,'std','mean','Letters',group)
  return(a)
}


# 3 -----------------------------------------------------------------------
ONE_LSD <- function(data,group,compare,value){
  library(agricolae)
  
  a <- data.frame(stringsAsFactors = F)
  type <- unique(data[,group])
  for (i in type)
  {
    # sub_dat <- subset(data,group == i)
    sub_dat <- data[data[,group]==i,]
    # fit <- aov(value ~ compare,sub_dat)
    fit <- aov(sub_dat[,value] ~ sub_dat[,compare] )
    out <- LSD.test(fit,'sub_dat[, compare]',p.adj='BH')#进行了p值校正
    #out$groups就可获取多重比较字母列表
    out$groups$type <- i
    out$groups$compare <- rownames(out$groups)
    
    a <- rbind(a,merge(out$means[,1:2], out$groups,by='sub_dat[, value]'))
  }
  names(a) <- c('mean','std','lsd',group,compare)
  return(a)
}

### abundance in clusters at individual level ####
ngs.cluster.change<-cbind(env,ngs.t)
ngs.cluster.change.pp<-subset(ngs.cluster.change,ngs.cluster.change$position != "cor" & ngs.cluster.change$position != "mb")
ngs.cluster.change.pp1<-subset(ngs.cluster.change.pp,ngs.cluster.change.pp$corridor==1)
ngs.cluster.change.pp0<-subset(ngs.cluster.change.pp,ngs.cluster.change.pp$corridor==0)
cluster.pp1<-ngs.cluster.change.pp1[,c(18:283)]
cluster.pp0<-ngs.cluster.change.pp0[,c(18:283)]

ngs.cluster.change.p1<-ngs.cluster.change.pp[,-c(58,78)]

letters<-as.data.frame(matrix(NA,nrow=10,ncol=50))
std<-as.data.frame(matrix(NA,nrow=10,ncol=50))
mean<-as.data.frame(matrix(NA,nrow=10,ncol=50))
for (i in 1:50)
{
letters[i]<-ONE_Tukey_HSD1(ngs.cluster.change.p1,'corridor','code',colnames(ngs.cluster.change.p1[17+i]))[,4]
std[i]<-ONE_Tukey_HSD1(ngs.cluster.change.p1,'corridor','code',colnames(ngs.cluster.change.p1[17+i]))[,2]
mean[i]<-ONE_Tukey_HSD1(ngs.cluster.change.p1,'corridor','code',colnames(ngs.cluster.change.p1[17+i]))[,3]
}

ngs.cluster.change.f<-do.call(rbind, replicate(10, ngs.cluster.change.p1[,c("code","corridor")], simplify=FALSE)) 
abundance<-stack(ngs.cluster.change.p1[,c(18:27)])
test.ab<-cbind(ngs.cluster.change.f,abundance)
tests.statis.ab<-ONE_Tukey_HSD1(test.ab,'corridor','code','values')

ggplot(ngs.cluster.change.f,aes(x=code,y=abundance))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  #ylim(50,180)+
  geom_text(data=stat.ab,aes(x=code,y=mean+1.3*std,label=Letters))+
  facet_wrap(.~corridor~cluster,scales = "free_y")+ labs(x='Different campaign',y='Abundance in cluster at individual level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

test.cluster<-ONE_Tukey_HSD1(ngs.cluster.change.p1,'corridor','code',colnames(ngs.cluster.change.pp[17+2]))[,4]

#### richness and evenness of all fungal cluster at individual,population and meta-population level####
#df1 <- ONE_Tukey_HSD1(data=dat,group='alpha',compare='Group',value='v')

richness.tri<-subset(richness,richness$plant=="Trifolium"&richness$position!="cor")

dfir.time <- ONE_Tukey_HSD1(richness.tri,'corridor','code','otu.richness')
dfir.corridor <- ONE_Tukey_HSD1(richness.tri,'code','corridor','otu.richness')
dfie.time <- ONE_Tukey_HSD1(richness.tri,'corridor','code','otu.evenness')
dfie.corridor <- ONE_Tukey_HSD1(richness.tri,'code','corridor','otu.evenness')

tri.r<-ggplot(data=richness.tri,aes(x=code,y=otu.richness))+
  geom_boxplot(aes(fill=as.factor(corridor)))+geom_text(data=dfir.time,aes(x=code,y=mean+1.3*std,label=Letters))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  xlab("Different Campaign") + 
  ylab("Individual Sample OTU richness")

df1.time <- ONE_Tukey_HSD1(richness.pp,'corridor','code','otu.richness')
df1.corridor <- ONE_Tukey_HSD1(richness.pp,'code','corridor','otu.richness')
df2.time <- ONE_Tukey_HSD1(richness.pp,'corridor','code','otu.evenness')
df2.corridor <- ONE_Tukey_HSD1(richness.pp,'code','corridor','otu.evenness')

pool.e<-ggplot(data=richness.pp,aes(x=Group.4,y=otu.richness))+
  geom_boxplot(aes(fill=as.factor(Group.6)))+geom_text(data=df1.time,aes(x=Group.4,y=mean+1.3*std,label=Letters))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.6,color=as.factor(Group.6)))+ 
  xlab("Different Campaign") + 
  ylab("Pool Sample OTU richness")

dfm1.time <- ONE_Tukey_HSD1(richness.ppool,'corridor','campaign','otu.richness.pp')
dfm1.corridor <- ONE_Tukey_HSD1(richness.ppool,'campaign','corridor','otu.richness.pp')
dfm2.time <- ONE_Tukey_HSD1(richness.ppool,'corridor','campaign','otu.evenness.pp')
dfm2.corridor <- ONE_Tukey_HSD1(richness.ppool,'campaign','corridor','otu.evenness.pp')

pdf("Figure S3.Community Richness and Evenness at three levels with statis.pdf",height=11,width=18.5)
p.rpp = ggplot(richness.ppool,aes(x=campaign,y=otu.richness.pp))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(50,180)+
  geom_text(data=dfm1.time,aes(x=campaign,y=mean+1.3*std,label=Letters))+
  #facet_wrap(.~corridor,scales = "free_y")+ 
  labs(x='Different campaign',y='Mycobiota OTU richness at meta-population level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

p.epp = ggplot(richness.ppool,aes(x=campaign,y=otu.evenness.pp))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0.4,0.7)+
  geom_text(data=dfm2.time,aes(x=campaign,y=mean+1.3*std,label=Letters))+
  #facet_wrap(.~corridor,scales = "free_y")+ 
  labs(x='Different campaign',y='Mycobiota OTU evenness at meta-population level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

p.rp = ggplot(richness.pp,aes(x=code,y=otu.richness))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(50,180)+
  geom_text(data=df1.time,aes(x=code,y=mean+1.3*std,label=Letters))+
  #facet_wrap(.~corridor,scales = "free_y")+ 
  labs(x='Different campaign',y='Mycobiota OTU richness at population level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

p.ep = ggplot(richness.pp,aes(x=code,y=otu.evenness))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0.4,0.7)+
  geom_text(data=df2.time,aes(x=code,y=mean+1.3*std,label=Letters))+
  #facet_wrap(.~corridor,scales = "free_y")+ 
  labs(x='Different campaign',y='Mycobiota OTU evenness at population level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

p.r = ggplot(richness.tri,aes(x=code,y=otu.richness))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(50,180)+
  geom_text(data=dfir.time,aes(x=code,y=mean+1.3*std,label=Letters))+
  #facet_wrap(.~corridor,scales = "free_y")+ 
  labs(x='Different campaign',y='Mycobiota OTU richness at individual level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

p.e = ggplot(richness.tri,aes(x=code,y=otu.evenness))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0.4,0.7)+
  geom_text(data=dfie.time,aes(x=code,y=mean+1.3*std,label=Letters))+
  #facet_wrap(.~corridor,scales = "free_y")+ 
  labs(x='Different campaign',y='Mycobiota OTU evenness at individual level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

grid.arrange(p.r,p.rp,p.rpp,p.e,p.ep,p.epp,ncol=3)

dev.off()  

corridor.sta<-cbind(dfir.corridor,df1.corridor,dfm1.corridor,dfie.corridor,df2.corridor,dfm1.corridor)
write.table(corridor.sta,file="corridor.sta.re.txt",sep="\t")

##### dissimilarity between patch 1 and patch 2 statistics ####

dist.corridor <- ONE_Tukey_HSD1(bray.i,'code','corridor','distance')
dist.time <- ONE_Tukey_HSD1(bray.i,'corridor','code','distance')

pdf("Figure.S7 Community dissimilarity at pairwise individual level.pdf",height=5,width=7)
dis.i = ggplot(bray.i,aes(x=code,y=distance))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0.1,1)+
  geom_text(data=dist.time,aes(x=code,y=mean+1.3*std,label=Letters))+
  #facet_wrap(.~corridor,scales = "free_y")+ 
  labs(x='Different campaign',y='Bray-Curtis dissimilarity between T. repens patch 1 and 2 at individual level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))
dis.i
dev.off()

statc.corridor <- ONE_Tukey_HSD1(dist.allc,'code','corridor','distance')
statc.time <- ONE_Tukey_HSD1(dist.allc,'corridor','code','distance')
state.corridor <- ONE_Tukey_HSD1(dist.alle,'code','corridor','distance')
state.time <- ONE_Tukey_HSD1(dist.alle,'corridor','code','distance')

install.packages("ggpubr")
library(ggpubr)
compare_means(distance~corridor, data = dist.allc, 
              group.by = "code")

dis.e = ggplot(dist.alle,aes(x=code,y=distance))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0.2,0.9)+
  geom_text(data=state.time,aes(x=code,y=mean+1.3*std,label=Letters))+
  facet_wrap(.~corridor,scales = "free_y")+ labs(x='Different campaign',y='Bray-Curtis dissimilarity between T. repens patch 1 and 2')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

pdf("Community dissimilarity at pairwise population level.pdf",height=6,width=8)
dis.c = ggplot(dist.allc,aes(x=code,y=distance))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),size=2.5,aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0.2,0.9)+
  geom_text(data=statc.time,aes(x=code,y=mean+1.3*std,label=Letters))+
  #facet_wrap(.~corridor,scales = "free_y")+
  labs(x='Different campaign',y='Bray-Curtis dissimilarity between patch 1 and 2')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle=45))
dis.c
dev.off() 

##### richness inside each phylum at meta-population level statistics####
statc.mac <- ONE_Tukey_HSD1(richness.pp.phylum,'Group.2','Group.3','otu.richness.asco.pp')
statc.mbc <- ONE_Tukey_HSD1(richness.pp.phylum,'Group.2','Group.3','otu.richness.Basidio.pp')
statc.mcc <- ONE_Tukey_HSD1(richness.pp.phylum,'Group.2','Group.3','otu.richness.Chytridio.pp')
statc.mgc <- ONE_Tukey_HSD1(richness.pp.phylum,'Group.2','Group.3','otu.richness.Glomero.pp')
statc.mzc <- ONE_Tukey_HSD1(richness.pp.phylum,'Group.2','Group.3','otu.richness.Zygo.pp')

statc.matime <- ONE_Tukey_HSD2(richness.pp.phylum,'Group.3','Group.2','otu.richness.asco.pp')
statc.mbtime <- ONE_Tukey_HSD2(richness.pp.phylum,'Group.3','Group.2','otu.richness.Basidio.pp')
statc.mctime <- ONE_Tukey_HSD2(richness.pp.phylum,'Group.3','Group.2','otu.richness.Chytridio.pp')
statc.mgtime <- ONE_Tukey_HSD2(richness.pp.phylum,'Group.3','Group.2','otu.richness.Glomero.pp')
statc.mztime <- ONE_Tukey_HSD2(richness.pp.phylum,'Group.3','Group.2','otu.richness.Zygo.pp')

pdf("Figure S4.Richness in phylum tri species at meta population level with statis.pdf",height=5,width=24)
a.pps<-ggplot(data=richness.pp.phylum,aes(x=Group.2,y=otu.richness.asco.pp))+geom_boxplot(aes(fill=as.factor(Group.3)))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.3,color=as.factor(Group.3)))+ 
  ylim(0,62)+ 
  geom_text(data=statc.matime,aes(x=Group.2,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Richness in Ascomycota at meta-population level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))
  
b.pps<-ggplot(data=richness.pp.phylum,aes(x=Group.2,y=otu.richness.Basidio.pp))+geom_boxplot(aes(fill=as.factor(Group.3)))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.3,color=as.factor(Group.3)))+ 
  ylim(0,62)+ 
  geom_text(data=statc.mbtime,aes(x=Group.2,y=mean+3*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Richness in Basidiomycota at meta-population level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

c.pps<-ggplot(data=richness.pp.phylum,aes(x=Group.2,y=otu.richness.Chytridio.pp))+geom_boxplot(aes(fill=as.factor(Group.3)))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.3,color=as.factor(Group.3)))+ 
  ylim(0,62)+ 
  geom_text(data=statc.mctime,aes(x=Group.2,y=mean+3*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Richness in Chytridiomycota at meta-population level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

g.pps<-ggplot(data=richness.pp.phylum,aes(x=Group.2,y=otu.richness.Glomero.pp))+geom_boxplot(aes(fill=as.factor(Group.3)))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.3,color=as.factor(Group.3)))+ 
  ylim(0,62)+ 
  geom_text(data=statc.mgtime,aes(x=Group.2,y=mean+4*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Richness in Glomeromycota at meta-population level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

z.pps<-ggplot(data=richness.pp.phylum,aes(x=Group.2,y=otu.richness.Zygo.pp))+geom_boxplot(aes(fill=as.factor(Group.3)))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.3,color=as.factor(Group.3)))+ 
  ylim(0,62)+ 
  geom_text(data=statc.mztime,aes(x=Group.2,y=mean+4*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Richness in Zygomycota at meta-population level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

grid.arrange(a.pps,b.pps,c.pps,g.pps,z.pps,ncol=5)
dev.off()

##### richness inside each phylum at population level statistics####
statc.pac <- ONE_Tukey_HSD1(richness.pool.phylum.tri,'Group.4','Group.6','otu.richness.asco.p')
statc.pbc <- ONE_Tukey_HSD1(richness.pool.phylum.tri,'Group.4','Group.6','otu.richness.Basidio.p')
statc.pcc <- ONE_Tukey_HSD1(richness.pool.phylum.tri,'Group.4','Group.6','otu.richness.Chytridio.p')
statc.pgc <- ONE_Tukey_HSD1(richness.pool.phylum.tri,'Group.4','Group.6','otu.richness.Glomero.p')
statc.pzc <- ONE_Tukey_HSD1(richness.pool.phylum.tri,'Group.4','Group.6','otu.richness.Zygo.p')

statc.patime <- ONE_Tukey_HSD1(richness.pool.phylum.tri,'Group.6','Group.4','otu.richness.asco.p')
statc.pbtime <- ONE_Tukey_HSD1(richness.pool.phylum.tri,'Group.6','Group.4','otu.richness.Basidio.p')
statc.pctime <- ONE_Tukey_HSD1(richness.pool.phylum.tri,'Group.6','Group.4','otu.richness.Chytridio.p')
statc.pgtime <- ONE_Tukey_HSD1(richness.pool.phylum.tri,'Group.6','Group.4','otu.richness.Glomero.p')
statc.pztime <- ONE_Tukey_HSD1(richness.pool.phylum.tri,'Group.6','Group.4','otu.richness.Zygo.p')

pdf("Figure S4.Richness in phylum tri species at population level with statis.pdf",height=5,width=24)

a.ps<-ggplot(data=richness.pool.phylum.tri,aes(x=Group.4,y=otu.richness.asco.p))+geom_boxplot(aes(fill=as.factor(Group.6)))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.6,color=as.factor(Group.6)))+ 
  ylim(0,55)+ 
  geom_text(data=statc.patime,aes(x=Group.4,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.6,ncol=1)+
  labs(x='Different campaign',y='Richness in Ascomycota at population level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

b.ps<-ggplot(data=richness.pool.phylum.tri,aes(x=Group.4,y=otu.richness.Basidio.p))+geom_boxplot(aes(fill=as.factor(Group.6)))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.6,color=as.factor(Group.6)))+ 
  ylim(0,55)+ 
  geom_text(data=statc.pbtime,aes(x=Group.4,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.6,ncol=1)+
  labs(x='Different campaign',y='Richness in Basidiomycota at population level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

c.ps<-ggplot(data=richness.pool.phylum.tri,aes(x=Group.4,y=otu.richness.Chytridio.p))+geom_boxplot(aes(fill=as.factor(Group.6)))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.6,color=as.factor(Group.6)))+ 
  ylim(0,55)+ 
  geom_text(data=statc.pctime,aes(x=Group.4,y=mean+3*std,label=Letters))+
  #facet_wrap(~Group.6,ncol=1)+
  labs(x='Different campaign',y='Richness in Chytridiomycota at population level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

g.ps<-ggplot(data=richness.pool.phylum.tri,aes(x=Group.4,y=otu.richness.Glomero.p))+geom_boxplot(aes(fill=as.factor(Group.6)))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.6,color=as.factor(Group.6)))+ 
  ylim(0,55)+ 
  geom_text(data=statc.pgtime,aes(x=Group.4,y=mean+2*std,label=Letters))+
  #facet_wrap(~Group.6,ncol=1)+
  labs(x='Different campaign',y='Richness in Glomeromycota at population level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

z.ps<-ggplot(data=richness.pool.phylum.tri,aes(x=Group.4,y=otu.richness.Zygo.p))+geom_boxplot(aes(fill=as.factor(Group.6)))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.6,color=as.factor(Group.6)))+ 
  ylim(0,55)+ 
  geom_text(data=statc.pztime,aes(x=Group.4,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.6,ncol=1)+
  labs(x='Different campaign',y='Richness in Zygomycota at population level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

grid.arrange(a.ps,b.ps,c.ps,g.ps,z.ps,ncol=5)
dev.off()

##### richness inside each phylum at individual level statistics####
statc.ac <- ONE_Tukey_HSD1(richness.with.phylum.nomb,'code','corridor','otu.richness.asco')
statc.bc <- ONE_Tukey_HSD1(richness.with.phylum.nomb,'code','corridor','otu.richness.Basidio')
statc.cc <- ONE_Tukey_HSD1(richness.with.phylum.nomb,'code','corridor','otu.richness.Chytridio')
statc.gc <- ONE_Tukey_HSD1(richness.with.phylum.nomb,'code','corridor','otu.richness.Glomero')
statc.zc <- ONE_Tukey_HSD1(richness.with.phylum.nomb,'code','corridor','otu.richness.Zygo')

statc.atime <- ONE_Tukey_HSD1(richness.with.phylum.nomb,'corridor','code','otu.richness.asco')
statc.btime <- ONE_Tukey_HSD1(richness.with.phylum.nomb,'corridor','code','otu.richness.Basidio')
statc.ctime <- ONE_Tukey_HSD1(richness.with.phylum.nomb,'corridor','code','otu.richness.Chytridio')
statc.gtime <- ONE_Tukey_HSD1(richness.with.phylum.nomb,'corridor','code','otu.richness.Glomero')
statc.ztime <- ONE_Tukey_HSD1(richness.with.phylum.nomb,'corridor','code','otu.richness.Zygo')

pdf("Figure S4.Richness in phylum at individual level with statis.pdf",height=5,width=24)

a.r<-ggplot(data=richness.with.phylum.nomb,aes(x=code,y=otu.richness.asco))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,55)+ 
  geom_text(data=statc.atime,aes(x=code,y=mean+2*std,label=Letters))+
  #facet_wrap(~corridor,ncol=1)+
  labs(x='Different campaign',y='Richness in Ascomycota at individual level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

b.r<-ggplot(data=richness.with.phylum.nomb,aes(x=code,y=otu.richness.Basidio))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,55)+ 
  geom_text(data=statc.btime,aes(x=code,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~corridor,ncol=1)+
  labs(x='Different campaign',y='Richness in Basidiomycota at individual level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

c.r<-ggplot(data=richness.with.phylum.nomb,aes(x=code,y=otu.richness.Chytridio))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,55)+ 
  geom_text(data=statc.ctime,aes(x=code,y=mean+3*std,label=Letters))+
  #facet_wrap(~corridor,ncol=1)+
  labs(x='Different campaign',y='Richness in Chytridiomycota at individual level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

g.r<-ggplot(data=richness.with.phylum.nomb,aes(x=code,y=otu.richness.Glomero))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,55)+ 
  geom_text(data=statc.gtime,aes(x=code,y=mean+4*std,label=Letters))+
  #facet_wrap(~corridor,ncol=1)+
  labs(x='Different campaign',y='Richness in Glomeromycota at individual level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

z.r<-ggplot(data=richness.with.phylum.nomb,aes(x=code,y=otu.richness.Zygo))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,55)+ 
  geom_text(data=statc.atime,aes(x=code,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~corridor,ncol=1)+
  labs(x='Different campaign',y='Richness in Zygomycota at individual level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

grid.arrange(a.r,b.r,c.r,g.r,z.r,ncol=5)
dev.off()

pdf("Figure S4.Richness in phylum with statis.pdf",height=10,width=22)

grid.arrange(a.r,b.r,c.r,g.r,z.r,a.ps,b.ps,c.ps,g.ps,z.ps,a.pps,b.pps,c.pps,g.pps,z.pps,ncol=5)

dev.off()


cor.sta.phy<-cbind(statc.ac,statc.bc,statc.cc,statc.gc,statc.zc,
                   statc.pac,statc.pbc,statc.pcc,statc.pgc,statc.pzc,
                   statc.mac,statc.mbc,statc.mcc,statc.mgc,statc.mzc)

write.table(cor.sta.phy,file="cor.sta.phy.txt",sep="\t")



##### evenness inside each phylum at patch level statistics####
statc.pac.e <- ONE_Tukey_HSD1(richness.pool.phylum.tri,'Group.4','Group.6','otu.evenness.asco.p')
statc.pbc.e <- ONE_Tukey_HSD1(richness.pool.phylum.tri,'Group.4','Group.6','otu.evenness.Basidio.p')
statc.pcc.e <- ONE_Tukey_HSD1(richness.pool.phylum.tri,'Group.4','Group.6','otu.evenness.Chytridio.p')
statc.pgc.e <- ONE_Tukey_HSD1(richness.pool.phylum.tri,'Group.4','Group.6','otu.evenness.Glomero.p')
statc.pzc.e <- ONE_Tukey_HSD1(richness.pool.phylum.tri,'Group.4','Group.6','otu.evenness.Zygo.p')

statc.patime.e <- ONE_Tukey_HSD1(richness.pool.phylum.tri,'Group.6','Group.4','otu.evenness.asco.p')
statc.pbtime.e <- ONE_Tukey_HSD1(richness.pool.phylum.tri,'Group.6','Group.4','otu.evenness.Basidio.p')
statc.pctime.e <- ONE_Tukey_HSD1(richness.pool.phylum.tri,'Group.6','Group.4','otu.evenness.Chytridio.p')
statc.pgtime.e <- ONE_Tukey_HSD1(richness.pool.phylum.tri,'Group.6','Group.4','otu.evenness.Glomero.p')
statc.pztime.e <- ONE_Tukey_HSD1(richness.pool.phylum.tri,'Group.6','Group.4','otu.evenness.Zygo.p')

pdf("Figure S4.Evenness in phylum tri species at population level with statis.pdf",height=5,width=24)

a.ps.e<-ggplot(data=richness.pool.phylum.tri,aes(x=Group.4,y=otu.evenness.asco.p))+geom_boxplot(aes(fill=as.factor(Group.6)))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.6,color=as.factor(Group.6)))+ 
  ylim(0,1)+ 
  geom_text(data=statc.pac.e,aes(x=Group.4,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.6,ncol=1)+
  labs(x='Different campaign',y='Evenness in Ascomycota at patch level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

b.ps.e<-ggplot(data=richness.pool.phylum.tri,aes(x=Group.4,y=otu.evenness.Basidio.p))+geom_boxplot(aes(fill=as.factor(Group.6)))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.6,color=as.factor(Group.6)))+ 
  ylim(0,1)+ 
  geom_text(data=statc.pbc.e,aes(x=Group.4,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.6,ncol=1)+
  labs(x='Different campaign',y='Evenness in Basidiomycota at patch level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

c.ps.e<-ggplot(data=richness.pool.phylum.tri,aes(x=Group.4,y=otu.evenness.Chytridio.p))+geom_boxplot(aes(fill=as.factor(Group.6)))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.6,color=as.factor(Group.6)))+ 
  ylim(0,1)+ 
  geom_text(data=statc.pcc.e,aes(x=Group.4,y=mean+3*std,label=Letters))+
  #facet_wrap(~Group.6,ncol=1)+
  labs(x='Different campaign',y='Evenness in Chytridiomycota at patch level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

g.ps.e<-ggplot(data=richness.pool.phylum.tri,aes(x=Group.4,y=otu.evenness.Glomero.p))+geom_boxplot(aes(fill=as.factor(Group.6)))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.6,color=as.factor(Group.6)))+ 
  ylim(0,1)+ 
  geom_text(data=statc.pgc.e,aes(x=Group.4,y=mean+2*std,label=Letters))+
  #facet_wrap(~Group.6,ncol=1)+
  labs(x='Different campaign',y='Evenness in Glomeromycota at patch level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

z.ps.e<-ggplot(data=richness.pool.phylum.tri,aes(x=Group.4,y=otu.evenness.Zygo.p))+geom_boxplot(aes(fill=as.factor(Group.6)))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.6,color=as.factor(Group.6)))+ 
  ylim(0,1)+ 
  geom_text(data=statc.pzc.e,aes(x=Group.4,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.6,ncol=1)+
  labs(x='Different campaign',y='Evenness in Zygomycota at patch level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

grid.arrange(a.ps.e,b.ps.e,c.ps.e,g.ps.e,z.ps.e,ncol=5)
dev.off()

##### richness inside each phylum at individual level statistics####
statc.ac.e <- ONE_Tukey_HSD1(richness.with.phylum.nomb,'code','corridor','otu.evenness.asco')
statc.bc.e <- ONE_Tukey_HSD1(richness.with.phylum.nomb,'code','corridor','otu.evenness.Basidio')
statc.cc.e <- ONE_Tukey_HSD1(richness.with.phylum.nomb,'code','corridor','otu.evenness.Chytridio')
statc.gc.e <- ONE_Tukey_HSD1(richness.with.phylum.nomb,'code','corridor','otu.evenness.Glomero')
statc.zc.e <- ONE_Tukey_HSD1(richness.with.phylum.nomb,'code','corridor','otu.evenness.Zygo')

statc.atime.e <- ONE_Tukey_HSD1(richness.with.phylum.nomb,'corridor','code','otu.evenness.asco')
statc.btime.e <- ONE_Tukey_HSD1(richness.with.phylum.nomb,'corridor','code','otu.evenness.Basidio')
statc.ctime.e <- ONE_Tukey_HSD1(richness.with.phylum.nomb,'corridor','code','otu.evenness.Chytridio')
statc.gtime.e <- ONE_Tukey_HSD1(richness.with.phylum.nomb,'corridor','code','otu.evenness.Glomero')
statc.ztime.e <- ONE_Tukey_HSD1(richness.with.phylum.nomb,'corridor','code','otu.evenness.Zygo')

pdf("Figure S4.Evenness in phylum at individual level with statis.pdf",height=5,width=24)

a.r.e<-ggplot(data=richness.with.phylum.nomb,aes(x=code,y=otu.evenness.asco))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,1)+ 
  geom_text(data=statc.ac.e,aes(x=code,y=mean+2*std,label=Letters))+
  #facet_wrap(~corridor,ncol=1)+
  labs(x='Different campaign',y='Evenness in Ascomycota at individual level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

b.r.e<-ggplot(data=richness.with.phylum.nomb,aes(x=code,y=otu.evenness.Basidio))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,1)+ 
  geom_text(data=statc.bc.e,aes(x=code,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~corridor,ncol=1)+
  labs(x='Different campaign',y='Evenness in Basidiomycota at individual level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

c.r.e<-ggplot(data=richness.with.phylum.nomb,aes(x=code,y=otu.evenness.Chytridio))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,1)+ 
  geom_text(data=statc.cc.e,aes(x=code,y=mean+3*std,label=Letters))+
  #facet_wrap(~corridor,ncol=1)+
  labs(x='Different campaign',y='Evenness in Chytridiomycota at individual level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

g.r.e<-ggplot(data=richness.with.phylum.nomb,aes(x=code,y=otu.evenness.Glomero))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,1)+ 
  geom_text(data=statc.gc.e,aes(x=code,y=mean+4*std,label=Letters))+
  #facet_wrap(~corridor,ncol=1)+
  labs(x='Different campaign',y='Evenness in Glomeromycota at individual level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

z.r.e<-ggplot(data=richness.with.phylum.nomb,aes(x=code,y=otu.evenness.Zygo))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,1)+ 
  geom_text(data=statc.ac.e,aes(x=code,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~corridor,ncol=1)+
  labs(x='Different campaign',y='Evenness in Zygomycota at individual level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

grid.arrange(a.r.e,b.r.e,c.r.e,g.r.e,z.r.e,ncol=5)
dev.off()

pdf("Figure S6+.Evenness in phylum with statis.pdf",height=6,width=22)

grid.arrange(a.r.e,b.r.e,c.r.e,g.r.e,z.r.e,a.ps.e,b.ps.e,c.ps.e,g.ps.e,z.ps.e,ncol=5)

dev.off()


cor.sta.phy<-cbind(statc.ac,statc.bc,statc.cc,statc.gc,statc.zc,
                   statc.pac,statc.pbc,statc.pcc,statc.pgc,statc.pzc,
                   statc.mac,statc.mbc,statc.mcc,statc.mgc,statc.mzc)

write.table(cor.sta.phy,file="cor.sta.phy.txt",sep="\t")


### abundance in phylum at individual, population and meta-population level

##### richness and abundance inside each FUNGuild statistics####
#####Richness at individual level
statc.pathc <- ONE_Tukey_HSD1(richness.with.guild,'corridor','code','otu.richness.path')
statc.amfc <- ONE_Tukey_HSD1(richness.with.guild,'corridor','code','otu.richness.amf')
statc.sapc <- ONE_Tukey_HSD1(richness.with.guild,'corridor','code','otu.richness.sap')

statc.pathtime <- ONE_Tukey_HSD1(richness.with.guild,'code','corridor','otu.richness.path')
statc.amftime <- ONE_Tukey_HSD1(richness.with.guild,'code','corridor','otu.richness.amf')
statc.saptime <- ONE_Tukey_HSD1(richness.with.guild,'code','corridor','otu.richness.sap')

pdf("Figure S8.Richness in FUNGUilds at individual level with statis.pdf",height=5,width=14)
path<-ggplot(data=richness.with.guild,aes(x=code,y=otu.richness.path))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,6)+ 
  geom_text(data=statc.pathtime,aes(x=code,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Richness in Plant Pathogen at individual level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

amf<-ggplot(data=richness.with.guild,aes(x=code,y=otu.richness.amf))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,6)+ 
  geom_text(data=statc.amftime,aes(x=code,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Richness in AMF at individual level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

sap<-ggplot(data=richness.with.guild,aes(x=code,y=otu.richness.sap))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,6)+ 
  geom_text(data=statc.saptime,aes(x=code,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Richness in Saprotroph at individual level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

grid.arrange(path,amf,sap,ncol=3)
dev.off()

#####Richness at patch level
statc.pathc.p <- ONE_Tukey_HSD1(richness.with.guild.p,'Group.6','Group.4','otu.richness.path.p')
statc.amfc.p <- ONE_Tukey_HSD1(richness.with.guild.p,'Group.6','Group.4','otu.richness.amf.p')
statc.sapc.p <- ONE_Tukey_HSD1(richness.with.guild.p,'Group.6','Group.4','otu.richness.sap.p')

statc.pathtime.p <- ONE_Tukey_HSD1(richness.with.guild.p,'Group.4','Group.6','otu.richness.path.p')
statc.amftime.p <- ONE_Tukey_HSD1(richness.with.guild.p,'Group.4','Group.6','otu.richness.amf.p')
statc.saptime.p <- ONE_Tukey_HSD1(richness.with.guild.p,'Group.4','Group.6','otu.richness.sap.p')

pdf("Figure S8.Richness in FUNGUilds at patch level with statis.pdf",height=5,width=14)
path.p<-ggplot(data=richness.with.guild.p,aes(x=Group.4,y=otu.richness.path.p))+geom_boxplot(aes(fill=as.factor(Group.6)))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.6,color=as.factor(Group.6)))+ 
  ylim(0,6)+ 
  geom_text(data=statc.pathtime.p,aes(x=Group.4,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Richness in Plant Pathogen at patch level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

amf.p<-ggplot(data=richness.with.guild.p,aes(x=Group.4,y=otu.richness.amf.p))+geom_boxplot(aes(fill=as.factor(Group.6)))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.6,color=as.factor(Group.6)))+ 
  ylim(0,6)+ 
  geom_text(data=statc.amftime.p,aes(x=Group.4,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Richness in AMF at patch level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

sap.p<-ggplot(data=richness.with.guild.p,aes(x=Group.4,y=otu.richness.sap.p))+geom_boxplot(aes(fill=as.factor(Group.6)))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.6,color=as.factor(Group.6)))+ 
  ylim(0,6)+ 
  geom_text(data=statc.saptime.p,aes(x=Group.4,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Richness in Saprotroph at patch level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

grid.arrange(path.p,amf.p,sap.p,ncol=3)
dev.off()

#####Richness at landscape level
statc.pathc.l <- ONE_Tukey_HSD1(richness.with.guild.l,'campaign','corridor','otu.richness.path.l')
statc.amfc.l <- ONE_Tukey_HSD1(richness.with.guild.l,'campaign','corridor','otu.richness.amf.l')
statc.sapc.l <- ONE_Tukey_HSD1(richness.with.guild.l,'campaign','corridor','otu.richness.sap.l')

statc.pathtime.l <- ONE_Tukey_HSD1(richness.with.guild.l,'corridor','campaign','otu.richness.path.l')
statc.amftime.l <- ONE_Tukey_HSD1(richness.with.guild.l,'corridor','campaign','otu.richness.amf.l')
statc.saptime.l <- ONE_Tukey_HSD1(richness.with.guild.l,'corridor','campaign','otu.richness.sap.l')

pdf("Figure S8.Richness in FUNGUilds at landscape level with statis.pdf",height=5,width=14)
path.l<-ggplot(data=richness.with.guild.l,aes(x=campaign,y=otu.richness.path.l))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,6)+ 
  geom_text(data=statc.pathtime.l,aes(x=campaign,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Richness in Plant Pathogen at landscape level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

amf.l<-ggplot(data=richness.with.guild.l,aes(x=campaign,y=otu.richness.amf.l))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,6)+ 
  geom_text(data=statc.amftime.l,aes(x=campaign,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Richness in AMF at patch level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

sap.l<-ggplot(data=richness.with.guild.l,aes(x=campaign,y=otu.richness.sap.l))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,6)+ 
  geom_text(data=statc.saptime.l,aes(x=campaign,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Richness in Saprotroph at landscape level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

grid.arrange(path.l,amf.l,sap.l,ncol=3)
dev.off()


##### Abundance at individual level ####
statc.pathac <- ONE_Tukey_HSD1(ngs.gui.t1,'corridor','code','Plant Pathogen')
statc.amfac <- ONE_Tukey_HSD1(ngs.gui.t1,'corridor','code','Arbuscular Mycorrhizal')
statc.sapac <- ONE_Tukey_HSD1(ngs.gui.t1,'corridor','code','Saprotroph')

statc.pathatime <- ONE_Tukey_HSD1(ngs.gui.t1,'code','corridor','Plant Pathogen')
statc.amfatime <- ONE_Tukey_HSD1(ngs.gui.t1,'code','corridor','Arbuscular Mycorrhizal')
statc.sapatime <- ONE_Tukey_HSD1(ngs.gui.t1,'code','corridor','Saprotroph')

pdf("Figure S9.Richness and abundance in FUNGUilds at individual level with statis.pdf",height=6,width=12)
path.a<-ggplot(data=ngs.gui.t1,aes(x=code,y=`Plant Pathogen`))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,100)+ 
  geom_text(data=statc.pathatime,aes(x=code,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Abundance in Plant Pathogen at individual level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

amf.a<-ggplot(data=ngs.gui.t1,aes(x=code,y=`Arbuscular Mycorrhizal`))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,400)+ 
  geom_text(data=statc.amfatime,aes(x=code,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Abundance in AMF at individual level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

sap.a<-ggplot(data=ngs.gui.t1,aes(x=code,y=`Saprotroph`))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,600)+ 
  geom_text(data=statc.sapatime,aes(x=code,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Abundance in Saprotroph at individual level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

grid.arrange(path,amf,sap,path.a,amf.a,sap.a,ncol=3)
dev.off()


##### Abundance at patch level ####
statc.pathac.p <- ONE_Tukey_HSD1(ngs.gui.p,'Group.3','Group.2','Plant Pathogen')
statc.amfac.p <- ONE_Tukey_HSD1(ngs.gui.p,'Group.3','Group.2','Arbuscular Mycorrhizal')
statc.sapac.p <- ONE_Tukey_HSD1(ngs.gui.p,'Group.3','Group.2','Saprotroph')

statc.pathatime.p <- ONE_Tukey_HSD1(ngs.gui.p,'Group.2','Group.3','Plant Pathogen')
statc.amfatime.p <- ONE_Tukey_HSD1(ngs.gui.p,'Group.2','Group.3','Arbuscular Mycorrhizal')
statc.sapatime.p <- ONE_Tukey_HSD1(ngs.gui.p,'Group.2','Group.3','Saprotroph')

pdf("Figure S9.Richness and abundance in FUNGUilds at patch level with statis.pdf",height=6,width=12)
path.ap<-ggplot(data=ngs.gui.p,aes(x=Group.2,y=`Plant Pathogen`))+geom_boxplot(aes(fill=as.factor(Group.3)))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.3,color=as.factor(Group.3)))+ 
  ylim(0,75)+ 
  geom_text(data=statc.pathatime.p,aes(x=Group.2,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Abundance in Plant Pathogen at patch level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

amf.ap<-ggplot(data=ngs.gui.l,aes(x=Group.2,y=`Arbuscular Mycorrhizal`))+geom_boxplot(aes(fill=as.factor(Group.3)))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.3,color=as.factor(Group.3)))+ 
  ylim(0,400)+ 
  geom_text(data=statc.amfatime.p,aes(x=Group.2,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Abundance in AMF at patch level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

sap.ap<-ggplot(data=ngs.gui.l,aes(x=Group.2,y=`Saprotroph`))+geom_boxplot(aes(fill=as.factor(Group.3)))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.3,color=as.factor(Group.3)))+ 
  ylim(0,600)+ 
  geom_text(data=statc.sapatime.p,aes(x=Group.2,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Abundance in Saprotroph at patch level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

grid.arrange(path,amf,sap,path.ap,amf.ap,sap.ap,ncol=3)
dev.off()


##### Abundance at landscape level ####
statc.pathac.l <- ONE_Tukey_HSD1(ngs.gui.l,'Group.3','Group.2','Plant Pathogen')
statc.amfac.l <- ONE_Tukey_HSD1(ngs.gui.l,'Group.3','Group.2','Arbuscular Mycorrhizal')
statc.sapac.l <- ONE_Tukey_HSD1(ngs.gui.l,'Group.3','Group.2','Saprotroph')

statc.pathatime.l <- ONE_Tukey_HSD1(ngs.gui.l,'Group.2','Group.3','Plant Pathogen')
statc.amfatime.l <- ONE_Tukey_HSD1(ngs.gui.l,'Group.2','Group.3','Arbuscular Mycorrhizal')
statc.sapatime.l <- ONE_Tukey_HSD1(ngs.gui.l,'Group.2','Group.3','Saprotroph')

pdf("Figure S9.Richness and abundance in FUNGUilds at landscape level with statis.pdf",height=6,width=12)
path.al<-ggplot(data=ngs.gui.l,aes(x=Group.2,y=`Plant Pathogen`))+geom_boxplot(aes(fill=as.factor(Group.3)))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.3,color=as.factor(Group.3)))+ 
  ylim(0,75)+ 
  geom_text(data=statc.pathatime.l,aes(x=Group.2,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Abundance in Plant Pathogen at landscape level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

amf.al<-ggplot(data=ngs.gui.l,aes(x=Group.2,y=`Arbuscular Mycorrhizal`))+geom_boxplot(aes(fill=as.factor(Group.3)))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.3,color=as.factor(Group.3)))+ 
  ylim(0,400)+ 
  geom_text(data=statc.amfatime.l,aes(x=Group.2,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Abundance in AMF at landscape level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

sap.al<-ggplot(data=ngs.gui.l,aes(x=Group.2,y=`Saprotroph`))+geom_boxplot(aes(fill=as.factor(Group.3)))+
  geom_point(position=position_dodge(width=0.75),aes(group=Group.3,color=as.factor(Group.3)))+ 
  ylim(0,600)+ 
  geom_text(data=statc.sapatime.l,aes(x=Group.2,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Abundance in Saprotroph at landscape level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

grid.arrange(path.a,amf.a,sap.a,path.ap,amf.ap,sap.ap,path.al,amf.al,sap.al,ncol=3)
dev.off()

pdf("Figure S9.Richness and abundance in Plant Pathogen at all levels with statis.pdf",height=8,width=8)
grid.arrange(path,path.a,path.p,path.ap,path.l,path.al,ncol=2)
dev.off()

#####Dissimilarity for each phylum at individual and patch levels####
dis.phylum.all1<-read.table("Bray curtis distance at individual level All phylum1.txt", header=T, sep="\t")  ## otu abundance rawdata

#### At individual level ####
statc.dispc.ia <- ONE_Tukey_HSD1(dis.phylum.all1,'time','corridor','Ascomycota')
statc.dispt.ia <- ONE_Tukey_HSD1(dis.phylum.all1,'corridor','time','Ascomycota')

dis.ai<-ggplot(data=dis.phylum.all1,aes(x=time,y=Ascomycota))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,1)+ 
  geom_text(data=statc.dispc.ia,aes(x=time,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Bray-Curtis didsimilarity for Ascomycota at individual level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

statc.dispc.ib <- ONE_Tukey_HSD1(dis.phylum.all1,'time','corridor','Basidiomycota')
statc.dispt.ib <- ONE_Tukey_HSD1(dis.phylum.all1,'corridor','time','Basidiomycota')

dis.bi<-ggplot(data=dis.phylum.all1,aes(x=time,y=Basidiomycota))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,1)+ 
  geom_text(data=statc.dispc.ib,aes(x=time,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Bray-Curtis didsimilarity for Basidiomycota at individual level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

statc.dispc.ic <- ONE_Tukey_HSD1(dis.phylum.all1,'time','corridor','Chytridiomycota')
statc.dispt.ic <- ONE_Tukey_HSD1(dis.phylum.all1,'corridor','time','Chytridiomycota')

dis.ci<-ggplot(data=dis.phylum.all1,aes(x=time,y=Chytridiomycota))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,1)+ 
  geom_text(data=statc.dispc.ic,aes(x=time,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Bray-Curtis didsimilarity for Chytridiomycota at individual level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))


statc.dispc.ig <- ONE_Tukey_HSD1(dis.phylum.all1,'time','corridor','Glomeromycota')
statc.dispt.ig <- ONE_Tukey_HSD1(dis.phylum.all1,'corridor','time','Glomeromycota')

dis.gi<-ggplot(data=dis.phylum.all1,aes(x=time,y=Glomeromycota))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,1)+ 
  geom_text(data=statc.dispc.ig,aes(x=time,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Bray-Curtis didsimilarity for Glomeromycota at individual level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))


statc.dispc.iz <- ONE_Tukey_HSD1(dis.phylum.all1,'time','corridor','Zygomycota')
statc.dispt.iz <- ONE_Tukey_HSD1(dis.phylum.all1,'corridor','time','Zygomycota')

dis.zi<-ggplot(data=dis.phylum.all1,aes(x=time,y=Zygomycota))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,1)+ 
  geom_text(data=statc.dispc.iz,aes(x=time,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Bray-Curtis didsimilarity for Zygomycota at individual level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

pdf("Figure.S7+.Bray-Curtis dissimilarity for all phylum at individual level with statis.pdf",height=5,width=25)
grid.arrange(dis.ai,dis.bi,dis.ci,dis.gi,dis.zi,ncol=5)
dev.off()

#### At patch level ####
dis.phylum.allp1<-read.table("Bray curtis distance at patch level All phylum1.txt", header=T, sep="\t")  ## otu abundance rawdata

statc.dispc.pa <- ONE_Tukey_HSD1(dis.phylum.allp1,'time','corridor','Ascomycota')
statc.dispt.pa <- ONE_Tukey_HSD1(dis.phylum.allp1,'corridor','time','Ascomycota')

dis.ap<-ggplot(data=dis.phylum.allp1,aes(x=time,y=Ascomycota))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,1)+ 
  geom_text(data=statc.dispc.pa,aes(x=time,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Bray-Curtis didsimilarity for Ascomycota at patch level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

statc.dispc.pb <- ONE_Tukey_HSD1(dis.phylum.allp1,'time','corridor','Basidiomycota')
statc.dispt.pb <- ONE_Tukey_HSD1(dis.phylum.allp1,'corridor','time','Basidiomycota')

dis.bp<-ggplot(data=dis.phylum.allp1,aes(x=time,y=Basidiomycota))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,1)+ 
  geom_text(data=statc.dispc.pb,aes(x=time,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Bray-Curtis didsimilarity for Basidiomycota at patch level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

statc.dispc.pc <- ONE_Tukey_HSD1(dis.phylum.allp1,'time','corridor','Chytridiomycota')
statc.dispt.pc <- ONE_Tukey_HSD1(dis.phylum.allp1,'corridor','time','Chytridiomycota')

dis.cp<-ggplot(data=dis.phylum.allp1,aes(x=time,y=Chytridiomycota))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,1)+ 
  geom_text(data=statc.dispc.pc,aes(x=time,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Bray-Curtis didsimilarity for Chytridiomycota at patch level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))


statc.dispc.pg <- ONE_Tukey_HSD1(dis.phylum.allp1,'time','corridor','Glomeromycota')
statc.dispt.pg <- ONE_Tukey_HSD1(dis.phylum.allp1,'corridor','time','Glomeromycota')

dis.gp<-ggplot(data=dis.phylum.allp1,aes(x=time,y=Glomeromycota))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,1)+ 
  geom_text(data=statc.dispc.pg,aes(x=time,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Bray-Curtis didsimilarity for Glomeromycota at patch level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))


statc.dispc.pz <- ONE_Tukey_HSD1(dis.phylum.allp1,'time','corridor','Zygomycota')
statc.dispt.pz <- ONE_Tukey_HSD1(dis.phylum.allp1,'corridor','time','Zygomycota')

dis.zp<-ggplot(data=dis.phylum.allp1,aes(x=time,y=Zygomycota))+geom_boxplot(aes(fill=as.factor(corridor)))+
  geom_point(position=position_dodge(width=0.75),aes(group=corridor,color=as.factor(corridor)))+ 
  ylim(0,1)+ 
  geom_text(data=statc.dispc.pz,aes(x=time,y=mean+1.5*std,label=Letters))+
  #facet_wrap(~Group.3,ncol=1)+
  labs(x='Different campaign',y='Bray-Curtis didsimilarity for Zygomycota at patch level')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))

pdf("Figure S7+.Bray-Curtis dissimilarity for all phylum at individual and patch level with statis.pdf",height=10,width=25)
grid.arrange(dis.ai,dis.bi,dis.ci,dis.gi,dis.zi,dis.ap,dis.bp,dis.cp,dis.gp,dis.zp,ncol=5)
dev.off()

