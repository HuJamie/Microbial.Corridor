#############
#############
library(vegan)
##### Bray-Curtis dissimilarity of all mycobiota At individual level####
ngs.rarefied

## calculate the Bray-curtis distance for each sampling campaign
ngs.t<-as.data.frame(ngs.rarefied.t)
big<-cbind(ngs.t,env)

#### 2017 October
ngs.2017.p1<-subset(big, big$corridor=="1" & big$code=="17O" & big$position.code=="patch")
dist.17o1<-as.matrix(vegdist(ngs.2017.p1[,1:266],method="bray"))

ngs.2017.p1$pool.code
dist.17o1[4:6,1:3]
dist.17o1[10:12,7:9]
dist.17o1[16:18,13:15]
dist.17o1[22:24,19:21]
dist.17o1[28:30,25:27]
dist.17o1[34:36,31:33]
dist.17o1[40:41,37:39]
dist.17o1[45:47,42:44]
dist.17o1[52:54,49:51]
dist.17o1[58,55:57]

doct.17.1<-as.vector(rbind(dist.17o1[4:6,1:3],dist.17o1[10:12,7:9],dist.17o1[16:18,13:15],dist.17o1[22:24,19:21],
                           dist.17o1[28:30,25:27],dist.17o1[34:36,31:33],dist.17o1[40:41,37:39],
                           dist.17o1[45:47,42:44],dist.17o1[52:54,49:51],dist.17o1[58,55:57]))

ngs.2017.p0<-subset(big, big$corridor=="0" & big$code=="17O"& big$position.code=="patch" )
dist.17o0<-as.matrix(vegdist(ngs.2017.p0[,1:266],method="bray"))

dist.17o0[4:6,1:3]
dist.17o0[10:12,7:9]
dist.17o0[16:18,13:15]
dist.17o0[22,19:21]
dist.17o0[26:28,23:25]
dist.17o0[32:34,29:31]
dist.17o0[38:40,35:37]
dist.17o0[44:46,41:43]
dist.17o0[50:52,47:49]

doct.17.0<-as.vector(rbind(dist.17o0[4:6,1:3],dist.17o0[10:12,7:9],dist.17o0[16:18,13:15],dist.17o0[22,19:21],
                           dist.17o0[26:28,23:25],dist.17o0[32:34,29:31],dist.17o0[38:40,35:37],
                           dist.17o0[44:46,41:43],dist.17o0[50:52,47:49]))

##### 2018 May
ngs.2018m.p1<-subset(big, big$corridor=="1" & big$code=="18 M" & big$position.code=="patch")
dist.18m1<-as.matrix(vegdist(ngs.2018m.p1[,1:266],method="bray"))

dist.18m1[1:3,4:6]
dist.18m1[7:9,10:11]
dist.18m1[12:14,15:17]
dist.18m1[18:20,21:23]
dist.18m1[24:26,27:29]
dist.18m1[30:31,32:34]
dist.18m1[35:37,38:40]
dist.18m1[41:43,44:46]
dist.18m1[47:49,50:52]
dist.18m1[53:55,56:58]

dmay.18.1<-as.vector(rbind(dist.18m1[1:3,4:6],dist.18m1[12:14,15:17],t(dist.18m1[7:9,10:11]),dist.18m1[18:20,21:23],
                 dist.18m1[24:26,27:29],dist.18m1[30:31,32:34],dist.18m1[35:37,38:40],dist.18m1[41:43,44:46],
                 dist.18m1[47:49,50:52],dist.18m1[53:55,56:58]))

ngs.2018m.p0<-subset(big, big$corridor=="0" & big$code=="18 M"& big$position.code=="patch" )
dist.18m0<-as.matrix(vegdist(ngs.2018m.p0[,1:266],method="bray"))

dist.18m0[1:3,4:6]
dist.18m0[7:9,10:12]
dist.18m0[13:15,16:18]
dist.18m0[19:21,22:24]
dist.18m0[25:27,28:30]
dist.18m0[31:32,33:35]
dist.18m0[36:38,39:41]
dist.18m0[42:44,45:47]
dist.18m0[48:50,51:53]

dmay.18.0<-as.vector(rbind(dist.18m0[1:3,4:6],dist.18m0[7:9,10:12],dist.18m0[13:15,16:18],dist.18m0[19:21,22:24],
                 dist.18m0[25:27,28:30],dist.18m0[31:32,33:35],dist.18m0[36:38,39:41],
                 dist.18m0[42:44,45:47],dist.18m0[48:50,51:53]))

##### 2018 June
ngs.2018j.p1<-subset(big, big$corridor=="1" & big$code=="18J" & big$position.code=="patch")
dist.18j1<-as.matrix(vegdist(ngs.2018j.p1[,1:266],method="bray"))

dist.18j1[1:3,4:6]
dist.18j1[7:9,10:12]
dist.18j1[13:15,16:18]
dist.18j1[19:21,22:23]
dist.18j1[24:26,27:29]
dist.18j1[30:31,32:34]
dist.18j1[35:37,38:40]
dist.18j1[41:43,44:46]
dist.18j1[47:49,50:52]
dist.18j1[53:55,56:58]

djune.18.1<-as.vector(cbind(dist.18j1[1:3,4:6],dist.18j1[7:9,10:12],dist.18j1[13:15,16:18],dist.18j1[19:21,22:23],
                  dist.18j1[24:26,27:29],dist.18j1[35:37,38:40],t(dist.18j1[30:31,32:34]),
                  dist.18j1[41:43,44:46],dist.18j1[47:49,50:52],dist.18j1[53:55,56:58]))

ngs.2018j.p0<-subset(big, big$corridor=="0" & big$code=="18J"& big$position.code=="patch" )
dist.18j0<-as.matrix(vegdist(ngs.2018j.p0[,1:266],method="bray"))

dist.18j0[1:3,4:6]
dist.18j0[7:9,10:12]
dist.18j0[13:15,16:18]
dist.18j0[19:21,22:24]
dist.18j0[25:27,28:30]
dist.18j0[31:32,33:35]
dist.18j0[36:38,39:41]
dist.18j0[42:44,45:47]
dist.18j0[48:50,51:53]

djune.18.0<-as.vector(rbind(dist.18j0[1:3,4:6],dist.18j0[7:9,10:12],dist.18j0[13:15,16:18],dist.18j0[19:21,22:24],
                  dist.18j0[25:27,28:30],dist.18j0[31:32,33:35],dist.18j0[36:38,39:41],
                  dist.18j0[42:44,45:47],dist.18j0[48:50,51:53]))

##### 2018 October
ngs.2018o.p1<-subset(big, big$corridor=="1" & big$code=="18O" & big$position.code=="patch")
dist.18o1<-as.matrix(vegdist(ngs.2018o.p1[,1:266],method="bray"))

doct.18.1<-as.vector(cbind(dist.18o1[1:3,4:6],dist.18o1[7:9,10:11],dist.18o1[12:14,15:17],dist.18o1[18:20,21:23],
                 dist.18o1[24:26,27:29],dist.18o1[30:32,33:35],dist.18o1[36:38,39:41],
                 dist.18o1[42:44,45:46],dist.18o1[47:49,50:52],dist.18o1[53:55,56:58]))

ngs.2018o.p0<-subset(big, big$corridor=="0" & big$code=="18O"& big$position.code=="patch" )
dist.18o0<-as.matrix(vegdist(ngs.2018o.p0[,1:266],method="bray"))

dist.18o0[1:3,4:6]
dist.18o0[7:8,9:11]
dist.18o0[12:14,15:17]
dist.18o0[18:20,21:23]
dist.18o0[24:26,27:29]
dist.18o0[30:32,33:35]
dist.18o0[36:38,39:41]
dist.18o0[42:44,45:47]
dist.18o0[48:50,51:53]

doct.18.0<-as.vector(rbind(dist.18o0[1:3,4:6],dist.18o0[7:8,9:11],dist.18o0[12:14,15:17],dist.18o0[18:20,21:23],
                dist.18o0[24:26,27:29],dist.18o0[30:32,33:35],dist.18o0[36:38,39:41],
                dist.18o0[42:44,45:47],dist.18o0[48:50,51:53]))

##### 2019 May
ngs.2019m.p1<-subset(big, big$corridor=="1" & big$code=="19M" & big$position.code=="patch")
dist.19m1<-as.matrix(vegdist(ngs.2019m.p1[,1:266],method="bray"))

dist.19m1[1:3,4:6]
dist.19m1[7:9,10:12]
dist.19m1[13:15,16:18]
##dist.19m1[19:21,]
dist.19m1[22:24,25:27]
dist.19m1[28:30,31:33]
dist.19m1[34:36,37:39]
dist.19m1[40:42,43:45]
dist.19m1[46:48,49:51]
dist.19m1[52:54,55:56]

dmay.19.1<-as.vector(cbind(dist.19m1[1:3,4:6],dist.19m1[7:9,10:12],dist.19m1[13:15,16:18],dist.19m1[22:24,25:27],
                dist.19m1[28:30,31:33],dist.19m1[34:36,37:39],dist.19m1[40:42,43:45],
                dist.19m1[46:48,49:51],dist.19m1[52:54,55:56]))
  
ngs.2019m.p0<-subset(big, big$corridor=="0" & big$code=="19M"& big$position.code=="patch" )
dist.19m0<-as.matrix(vegdist(ngs.2019m.p0[,1:266],method="bray"))

dist.19m0[1:3,4:6]
dist.19m0[7:9,10:12]
dist.19m0[13:15,16:18]
dist.19m0[19:21,22:24]
dist.19m0[25:27,28:30]
dist.19m0[31:33,34]
dist.19m0[35:37,38:40]
dist.19m0[41:43,44:46]

dmay.19.0<-as.vector(rbind(dist.19m0[1:3,4:6],dist.19m0[7:9,10:12],dist.19m0[13:15,16:18],dist.19m0[19:21,22:24],
      dist.19m0[25:27,28:30],dist.19m0[31:33,34],dist.19m0[35:37,38:40],dist.19m0[41:43,44:46]))


test<-as.data.frame(cbind(doct.17.0,doct.17.1,dmay.18.0,dmay.18.1,djune.18.0,djune.18.1,doct.18.0,doct.18.1,dmay.19.0,dmay.19.1))
test1<-stack(test)
write.table(test1,file="Bray curtis distance at individual level.txt",sep="\t")
bray.i<-read.table("Bray curtis distance at individual level.txt", header=T, row.names=1, sep="\t")  ## otu abundance rawdata


#### Bray-Curtis dissimilarity at each phylum level at individual level ####

####Ascomycota####
ngs.Asco

## calculate the Bray-Curtis distance for each sampling campaign
ngs.at<-as.data.frame(t(ngs.Asco[,-1]))
big.a<-cbind(ngs.at,env)

#### 2017 October
Asco.2017.p1<-subset(big.a, big.a$corridor=="1" & big.a$code=="17O" & big.a$position.code=="patch")
dist.a17o1<-as.matrix(vegdist(Asco.2017.p1[,1:79],method="bray"))

aoct.17.1<-as.vector(rbind(dist.a17o1[4:6,1:3],dist.a17o1[10:12,7:9],dist.a17o1[16:18,13:15],dist.a17o1[22:24,19:21],
                           dist.a17o1[28:30,25:27],dist.a17o1[34:36,31:33],dist.a17o1[40:41,37:39],
                           dist.a17o1[45:47,42:44],dist.a17o1[52:54,49:51],dist.a17o1[58,55:57]))

Asco.2017.p0<-subset(big.a, big.a$corridor=="0" & big.a$code=="17O"& big.a$position.code=="patch" )
dist.a17o0<-as.matrix(vegdist(Asco.2017.p0[,1:79],method="bray"))

aoct.17.0<-as.vector(rbind(dist.a17o0[4:6,1:3],dist.a17o0[10:12,7:9],dist.a17o0[16:18,13:15],dist.a17o0[22,19:21],
                           dist.a17o0[26:28,23:25],dist.a17o0[32:34,29:31],dist.a17o0[38:40,35:37],
                           dist.a17o0[44:46,41:43],dist.a17o0[50:52,47:49]))

##### 2018 May
Asco.2018m.p1<-subset(big.a, big.a$corridor=="1" & big.a$code=="18 M" & big.a$position.code=="patch")
dist.a18m1<-as.matrix(vegdist(Asco.2018m.p1[,1:79],method="bray"))

amay.18.1<-as.vector(rbind(dist.a18m1[1:3,4:6],dist.a18m1[12:14,15:17],t(dist.a18m1[7:9,10:11]),dist.a18m1[18:20,21:23],
                           dist.a18m1[24:26,27:29],dist.a18m1[30:31,32:34],dist.a18m1[35:37,38:40],dist.a18m1[41:43,44:46],
                           dist.a18m1[47:49,50:52],dist.a18m1[53:55,56:58]))

Asco.2018m.p0<-subset(big.a, big.a$corridor=="0" & big.a$code=="18 M"& big.a$position.code=="patch" )
dist.a18m0<-as.matrix(vegdist(Asco.2018m.p0[,1:79],method="bray"))

amay.18.0<-as.vector(rbind(dist.a18m0[1:3,4:6],dist.a18m0[7:9,10:12],dist.a18m0[13:15,16:18],dist.a18m0[19:21,22:24],
                           dist.a18m0[25:27,28:30],dist.a18m0[31:32,33:35],dist.a18m0[36:38,39:41],
                           dist.a18m0[42:44,45:47],dist.a18m0[48:50,51:53]))

##### 2018 June
Asco.2018j.p1<-subset(big.a, big.a$corridor=="1" & big.a$code=="18J" & big.a$position.code=="patch")
dist.a18j1<-as.matrix(vegdist(Asco.2018j.p1[,1:79],method="bray"))

ajune.18.1<-as.vector(cbind(dist.a18j1[1:3,4:6],dist.a18j1[7:9,10:12],dist.a18j1[13:15,16:18],dist.a18j1[19:21,22:23],
                            dist.a18j1[24:26,27:29],dist.a18j1[35:37,38:40],t(dist.a18j1[30:31,32:34]),
                            dist.a18j1[41:43,44:46],dist.a18j1[47:49,50:52],dist.a18j1[53:55,56:58]))

Asco.2018j.p0<-subset(big.a, big.a$corridor=="0" & big.a$code=="18J"& big.a$position.code=="patch" )
dist.a18j0<-as.matrix(vegdist(Asco.2018j.p0[,1:79],method="bray"))

ajune.18.0<-as.vector(rbind(dist.a18j0[1:3,4:6],dist.a18j0[7:9,10:12],dist.a18j0[13:15,16:18],dist.a18j0[19:21,22:24],
                            dist.a18j0[25:27,28:30],dist.a18j0[31:32,33:35],dist.a18j0[36:38,39:41],
                            dist.a18j0[42:44,45:47],dist.a18j0[48:50,51:53]))

##### 2018 October
Asco.2018o.p1<-subset(big.a, big.a$corridor=="1" & big.a$code=="18O" & big.a$position.code=="patch")
dist.a18o1<-as.matrix(vegdist(Asco.2018o.p1[,1:79],method="bray"))

aoct.18.1<-as.vector(cbind(dist.a18o1[1:3,4:6],dist.a18o1[7:9,10:11],dist.a18o1[12:14,15:17],dist.a18o1[18:20,21:23],
                           dist.a18o1[24:26,27:29],dist.a18o1[30:32,33:35],dist.a18o1[36:38,39:41],
                           dist.a18o1[42:44,45:46],dist.a18o1[47:49,50:52],dist.a18o1[53:55,56:58]))

Asco.2018o.p0<-subset(big.a, big.a$corridor=="0" & big.a$code=="18O"& big.a$position.code=="patch" )
dist.a18o0<-as.matrix(vegdist(Asco.2018o.p0[,1:79],method="bray"))

aoct.18.0<-as.vector(rbind(dist.a18o0[1:3,4:6],dist.a18o0[7:8,9:11],dist.a18o0[12:14,15:17],dist.a18o0[18:20,21:23],
                           dist.a18o0[24:26,27:29],dist.a18o0[30:32,33:35],dist.a18o0[36:38,39:41],
                           dist.a18o0[42:44,45:47],dist.a18o0[48:50,51:53]))

##### 2019 May
Asco.2019m.p1<-subset(big.a, big.a$corridor=="1" & big.a$code=="19M" & big.a$position.code=="patch")
dist.a19m1<-as.matrix(vegdist(Asco.2019m.p1[,1:79],method="bray"))

amay.19.1<-as.vector(cbind(dist.a19m1[1:3,4:6],dist.a19m1[7:9,10:12],dist.a19m1[13:15,16:18],dist.a19m1[22:24,25:27],
                           dist.a19m1[28:30,31:33],dist.a19m1[34:36,37:39],dist.a19m1[40:42,43:45],
                           dist.a19m1[46:48,49:51],dist.a19m1[52:54,55:56]))

Asco.2019m.p0<-subset(big.a, big.a$corridor=="0" & big.a$code=="19M"& big.a$position.code=="patch" )
dist.a19m0<-as.matrix(vegdist(Asco.2019m.p0[,1:79],method="bray"))

amay.19.0<-as.vector(rbind(dist.a19m0[1:3,4:6],dist.a19m0[7:9,10:12],dist.a19m0[13:15,16:18],dist.a19m0[19:21,22:24],
                           dist.a19m0[25:27,28:30],dist.a19m0[31:33,34],dist.a19m0[35:37,38:40],dist.a19m0[41:43,44:46]))


dis.Asco<-as.data.frame(cbind(aoct.17.0,aoct.17.1,amay.18.0,amay.18.1,ajune.18.0,ajune.18.1,aoct.18.0,aoct.18.1,amay.19.0,amay.19.1))
### number of rows of result is not a multiple of vector length (arg 1)
dis.Asco1<-stack(dis.Asco)
write.table(dis.Asco1,file="Bray curtis distance at individual level Asco.txt",sep="\t")
dis.Asco.i<-read.table("Bray curtis distance at individual level Asco.txt", header=T, row.names=1, sep="\t")  ## otu abundance rawdata

####Basidiomycota####
ngs.Basidio

## calculate the Bray-Curtis distance for each sampling campaign
ngs.bt<-as.data.frame(t(ngs.Basidio[,-1]))
big.b<-cbind(ngs.bt,env)

#### 2017 October
Basidio.2017.p1<-subset(big.b, big.b$corridor=="1" & big.b$code=="17O" & big.b$position.code=="patch")
dist.b17o1<-as.matrix(vegdist(Basidio.2017.p1[,1:64],method="bray"))

boct.17.1<-as.vector(rbind(dist.b17o1[4:6,1:3],dist.b17o1[10:12,7:9],dist.b17o1[16:18,13:15],dist.b17o1[22:24,19:21],
                           dist.b17o1[28:30,25:27],dist.b17o1[34:36,31:33],dist.b17o1[40:41,37:39],
                           dist.b17o1[45:47,42:44],dist.b17o1[52:54,49:51],dist.b17o1[58,55:57]))

Basidio.2017.p0<-subset(big.b, big.b$corridor=="0" & big.b$code=="17O"& big.b$position.code=="patch" )
dist.b17o0<-as.matrix(vegdist(Basidio.2017.p0[,1:64],method="bray"))

boct.17.0<-as.vector(rbind(dist.b17o0[4:6,1:3],dist.b17o0[10:12,7:9],dist.b17o0[16:18,13:15],dist.b17o0[22,19:21],
                           dist.b17o0[26:28,23:25],dist.b17o0[32:34,29:31],dist.b17o0[38:40,35:37],
                           dist.b17o0[44:46,41:43],dist.b17o0[50:52,47:49]))

##### 2018 May
Basidio.2018m.p1<-subset(big.b, big.b$corridor=="1" & big.b$code=="18 M" & big.b$position.code=="patch")
dist.b18m1<-as.matrix(vegdist(Basidio.2018m.p1[,1:64],method="bray"))

bmay.18.1<-as.vector(rbind(dist.b18m1[1:3,4:6],dist.b18m1[12:14,15:17],t(dist.b18m1[7:9,10:11]),dist.b18m1[18:20,21:23],
                           dist.b18m1[24:26,27:29],dist.b18m1[30:31,32:34],dist.b18m1[35:37,38:40],dist.b18m1[41:43,44:46],
                           dist.b18m1[47:49,50:52],dist.b18m1[53:55,56:58]))

Basidio.2018m.p0<-subset(big.b, big.b$corridor=="0" & big.b$code=="18 M"& big.b$position.code=="patch" )
dist.b18m0<-as.matrix(vegdist(Basidio.2018m.p0[,1:64],method="bray"))

bmay.18.0<-as.vector(rbind(dist.b18m0[1:3,4:6],dist.b18m0[7:9,10:12],dist.b18m0[13:15,16:18],dist.b18m0[19:21,22:24],
                           dist.b18m0[25:27,28:30],dist.b18m0[31:32,33:35],dist.b18m0[36:38,39:41],
                           dist.b18m0[42:44,45:47],dist.b18m0[48:50,51:53]))

##### 2018 June
Basidio.2018j.p1<-subset(big.b, big.b$corridor=="1" & big.b$code=="18J" & big.b$position.code=="patch")
dist.b18j1<-as.matrix(vegdist(Basidio.2018j.p1[,1:64],method="bray"))

bjune.18.1<-as.vector(cbind(dist.b18j1[1:3,4:6],dist.b18j1[7:9,10:12],dist.b18j1[13:15,16:18],dist.b18j1[19:21,22:23],
                            dist.b18j1[24:26,27:29],dist.b18j1[35:37,38:40],t(dist.b18j1[30:31,32:34]),
                            dist.b18j1[41:43,44:46],dist.b18j1[47:49,50:52],dist.b18j1[53:55,56:58]))

Basidio.2018j.p0<-subset(big.b, big.b$corridor=="0" & big.b$code=="18J"& big.b$position.code=="patch" )
dist.b18j0<-as.matrix(vegdist(Basidio.2018j.p0[,1:64],method="bray"))

bjune.18.0<-as.vector(rbind(dist.b18j0[1:3,4:6],dist.b18j0[7:9,10:12],dist.b18j0[13:15,16:18],dist.b18j0[19:21,22:24],
                            dist.b18j0[25:27,28:30],dist.b18j0[31:32,33:35],dist.b18j0[36:38,39:41],
                            dist.b18j0[42:44,45:47],dist.b18j0[48:50,51:53]))

##### 2018 October
Basidio.2018o.p1<-subset(big.b, big.b$corridor=="1" & big.b$code=="18O" & big.b$position.code=="patch")
dist.b18o1<-as.matrix(vegdist(Basidio.2018o.p1[,1:64],method="bray"))

boct.18.1<-as.vector(cbind(dist.b18o1[1:3,4:6],dist.b18o1[7:9,10:11],dist.b18o1[12:14,15:17],dist.b18o1[18:20,21:23],
                           dist.b18o1[24:26,27:29],dist.b18o1[30:32,33:35],dist.b18o1[36:38,39:41],
                           dist.b18o1[42:44,45:46],dist.b18o1[47:49,50:52],dist.b18o1[53:55,56:58]))

Basidio.2018o.p0<-subset(big.b, big.b$corridor=="0" & big.b$code=="18O"& big.b$position.code=="patch" )
dist.b18o0<-as.matrix(vegdist(Basidio.2018o.p0[,1:64],method="bray"))

boct.18.0<-as.vector(rbind(dist.b18o0[1:3,4:6],dist.b18o0[7:8,9:11],dist.b18o0[12:14,15:17],dist.b18o0[18:20,21:23],
                           dist.b18o0[24:26,27:29],dist.b18o0[30:32,33:35],dist.b18o0[36:38,39:41],
                           dist.b18o0[42:44,45:47],dist.b18o0[48:50,51:53]))

##### 2019 May
Basidio.2019m.p1<-subset(big.b, big.b$corridor=="1" & big.b$code=="19M" & big.b$position.code=="patch")
dist.b19m1<-as.matrix(vegdist(Basidio.2019m.p1[,1:64],method="bray"))

bmay.19.1<-as.vector(cbind(dist.b19m1[1:3,4:6],dist.b19m1[7:9,10:12],dist.b19m1[13:15,16:18],dist.b19m1[22:24,25:27],
                           dist.b19m1[28:30,31:33],dist.b19m1[34:36,37:39],dist.b19m1[40:42,43:45],
                           dist.b19m1[46:48,49:51],dist.b19m1[52:54,55:56]))

Basidio.2019m.p0<-subset(big.b, big.b$corridor=="0" & big.b$code=="19M"& big.b$position.code=="patch" )
dist.b19m0<-as.matrix(vegdist(Basidio.2019m.p0[,1:64],method="bray"))

bmay.19.0<-as.vector(rbind(dist.b19m0[1:3,4:6],dist.b19m0[7:9,10:12],dist.b19m0[13:15,16:18],dist.b19m0[19:21,22:24],
                           dist.b19m0[25:27,28:30],dist.b19m0[31:33,34],dist.b19m0[35:37,38:40],dist.b19m0[41:43,44:46]))


dis.Basidio<-as.data.frame(cbind(boct.17.0,boct.17.1,bmay.18.0,bmay.18.1,bjune.18.0,bjune.18.1,boct.18.0,boct.18.1,bmay.19.0,bmay.19.1))
dis.Basidio1<-stack(dis.Basidio)
write.table(dis.Basidio1,file="Bray curtis distance at individual level Basidio.txt",sep="\t")
dis.Basidio.i<-read.table("Bray curtis distance at individual level Basidio.txt", header=T, row.names=1, sep="\t")  ## otu abundance rawdata


####Chytridiomycota####
ngs.Chytridio

## calculate the Bray-Curtis distance for each sampling campaign
ngs.ct<-as.data.frame(t(ngs.Chytridio[,-1]))
big.c<-cbind(ngs.ct,env)

#### 2017 October
Chytridio.2017.p1<-subset(big.c, big.c$corridor=="1" & big.c$code=="17O" & big.c$position.code=="patch")
dist.c17o1<-as.matrix(vegdist(Chytridio.2017.p1[,1:68],method="bray"))

coct.17.1<-as.vector(rbind(dist.c17o1[4:6,1:3],dist.c17o1[10:12,7:9],dist.c17o1[16:18,13:15],dist.c17o1[22:24,19:21],
                           dist.c17o1[28:30,25:27],dist.c17o1[34:36,31:33],dist.c17o1[40:41,37:39],
                           dist.c17o1[45:47,42:44],dist.c17o1[52:54,49:51],dist.c17o1[58,55:57]))

Chytridio.2017.p0<-subset(big.c, big.c$corridor=="0" & big.c$code=="17O"& big.c$position.code=="patch" )
dist.c17o0<-as.matrix(vegdist(Chytridio.2017.p0[,1:68],method="bray"))

coct.17.0<-as.vector(rbind(dist.c17o0[4:6,1:3],dist.c17o0[10:12,7:9],dist.c17o0[16:18,13:15],dist.c17o0[22,19:21],
                           dist.c17o0[26:28,23:25],dist.c17o0[32:34,29:31],dist.c17o0[38:40,35:37],
                           dist.c17o0[44:46,41:43],dist.c17o0[50:52,47:49]))

##### 2018 May
Chytridio.2018m.p1<-subset(big.c, big.c$corridor=="1" & big.c$code=="18 M" & big.c$position.code=="patch")
dist.c18m1<-as.matrix(vegdist(Chytridio.2018m.p1[,1:68],method="bray"))

cmay.18.1<-as.vector(rbind(dist.c18m1[1:3,4:6],dist.c18m1[12:14,15:17],t(dist.c18m1[7:9,10:11]),dist.c18m1[18:20,21:23],
                           dist.c18m1[24:26,27:29],dist.c18m1[30:31,32:34],dist.c18m1[35:37,38:40],dist.c18m1[41:43,44:46],
                           dist.c18m1[47:49,50:52],dist.c18m1[53:55,56:58]))

Chytridio.2018m.p0<-subset(big.c, big.c$corridor=="0" & big.c$code=="18 M"& big.c$position.code=="patch" )
dist.c18m0<-as.matrix(vegdist(Chytridio.2018m.p0[,1:68],method="bray"))

cmay.18.0<-as.vector(rbind(dist.c18m0[1:3,4:6],dist.c18m0[7:9,10:12],dist.c18m0[13:15,16:18],dist.c18m0[19:21,22:24],
                           dist.c18m0[25:27,28:30],dist.c18m0[31:32,33:35],dist.c18m0[36:38,39:41],
                           dist.c18m0[42:44,45:47],dist.c18m0[48:50,51:53]))

##### 2018 June
Chytridio.2018j.p1<-subset(big.c, big.c$corridor=="1" & big.c$code=="18J" & big.c$position.code=="patch")
dist.c18j1<-as.matrix(vegdist(Chytridio.2018j.p1[,1:68],method="bray"))

cjune.18.1<-as.vector(cbind(dist.c18j1[1:3,4:6],dist.c18j1[7:9,10:12],dist.c18j1[13:15,16:18],dist.c18j1[19:21,22:23],
                            dist.c18j1[24:26,27:29],dist.c18j1[35:37,38:40],t(dist.c18j1[30:31,32:34]),
                            dist.c18j1[41:43,44:46],dist.c18j1[47:49,50:52],dist.c18j1[53:55,56:58]))

Chytridio.2018j.p0<-subset(big.c, big.c$corridor=="0" & big.c$code=="18J"& big.c$position.code=="patch" )
dist.c18j0<-as.matrix(vegdist(Chytridio.2018j.p0[,1:68],method="bray"))

cjune.18.0<-as.vector(rbind(dist.c18j0[1:3,4:6],dist.c18j0[7:9,10:12],dist.c18j0[13:15,16:18],dist.c18j0[19:21,22:24],
                            dist.c18j0[25:27,28:30],dist.c18j0[31:32,33:35],dist.c18j0[36:38,39:41],
                            dist.c18j0[42:44,45:47],dist.c18j0[48:50,51:53]))

##### 2018 October
Chytridio.2018o.p1<-subset(big.c, big.c$corridor=="1" & big.c$code=="18O" & big.c$position.code=="patch")
dist.c18o1<-as.matrix(vegdist(Chytridio.2018o.p1[,1:68],method="bray"))

coct.18.1<-as.vector(cbind(dist.c18o1[1:3,4:6],dist.c18o1[7:9,10:11],dist.c18o1[12:14,15:17],dist.c18o1[18:20,21:23],
                           dist.c18o1[24:26,27:29],dist.c18o1[30:32,33:35],dist.c18o1[36:38,39:41],
                           dist.c18o1[42:44,45:46],dist.c18o1[47:49,50:52],dist.c18o1[53:55,56:58]))

Chytridio.2018o.p0<-subset(big.c, big.c$corridor=="0" & big.c$code=="18O"& big.c$position.code=="patch" )
dist.c18o0<-as.matrix(vegdist(Chytridio.2018o.p0[,1:68],method="bray"))

coct.18.0<-as.vector(rbind(dist.c18o0[1:3,4:6],dist.c18o0[7:8,9:11],dist.c18o0[12:14,15:17],dist.c18o0[18:20,21:23],
                           dist.c18o0[24:26,27:29],dist.c18o0[30:32,33:35],dist.c18o0[36:38,39:41],
                           dist.c18o0[42:44,45:47],dist.c18o0[48:50,51:53]))

##### 2019 May
Chytridio.2019m.p1<-subset(big.c, big.c$corridor=="1" & big.c$code=="19M" & big.c$position.code=="patch")
dist.c19m1<-as.matrix(vegdist(Chytridio.2019m.p1[,1:68],method="bray"))

cmay.19.1<-as.vector(cbind(dist.c19m1[1:3,4:6],dist.c19m1[7:9,10:12],dist.c19m1[13:15,16:18],dist.c19m1[22:24,25:27],
                           dist.c19m1[28:30,31:33],dist.c19m1[34:36,37:39],dist.c19m1[40:42,43:45],
                           dist.c19m1[46:48,49:51],dist.c19m1[52:54,55:56]))

Chytridio.2019m.p0<-subset(big.c, big.c$corridor=="0" & big.c$code=="19M"& big.c$position.code=="patch" )
dist.c19m0<-as.matrix(vegdist(Chytridio.2019m.p0[,1:68],method="bray"))

cmay.19.0<-as.vector(rbind(dist.c19m0[1:3,4:6],dist.c19m0[7:9,10:12],dist.c19m0[13:15,16:18],dist.c19m0[19:21,22:24],
                           dist.c19m0[25:27,28:30],dist.c19m0[31:33,34],dist.c19m0[35:37,38:40],dist.c19m0[41:43,44:46]))


dis.Chytridio<-as.data.frame(cbind(coct.17.0,coct.17.1,cmay.18.0,cmay.18.1,cjune.18.0,cjune.18.1,coct.18.0,coct.18.1,cmay.19.0,cmay.19.1))
dis.Chytridio1<-stack(dis.Chytridio)
write.table(dis.Chytridio,file="Bray curtis distance at individual level Chytridio.txt",sep="\t")
dis.Chytridio.i<-read.table("Bray curtis distance at individual level Chytridio.txt", header=T, row.names=1, sep="\t")  ## otu abundance rawdata

####Glomeromycota####
ngs.Glomero

## calculate the Bray-Curtis distance for each sampling campaign
ngs.gt<-as.data.frame(t(ngs.Glomero[,-1]))
big.g<-cbind(ngs.gt,env)

#### 2017 October
Glomero.2017.p1<-subset(big.g, big.g$corridor=="1" & big.g$code=="17O" & big.g$position.code=="patch")
dist.g17o1<-as.matrix(vegdist(Glomero.2017.p1[,1:18],method="bray"))

goct.17.1<-as.vector(rbind(dist.g17o1[4:6,1:3],dist.g17o1[10:12,7:9],dist.g17o1[16:18,13:15],dist.g17o1[22:24,19:21],
                           dist.g17o1[28:30,25:27],dist.g17o1[34:36,31:33],dist.g17o1[40:41,37:39],
                           dist.g17o1[45:47,42:44],dist.g17o1[52:54,49:51],dist.g17o1[58,55:57]))

Glomero.2017.p0<-subset(big.g, big.g$corridor=="0" & big.g$code=="17O"& big.g$position.code=="patch" )
dist.g17o0<-as.matrix(vegdist(Glomero.2017.p0[,1:18],method="bray"))

goct.17.0<-as.vector(rbind(dist.g17o0[4:6,1:3],dist.g17o0[10:12,7:9],dist.g17o0[16:18,13:15],dist.g17o0[22,19:21],
                           dist.g17o0[26:28,23:25],dist.g17o0[32:34,29:31],dist.g17o0[38:40,35:37],
                           dist.g17o0[44:46,41:43],dist.g17o0[50:52,47:49]))

##### 2018 May
Glomero.2018m.p1<-subset(big.g, big.g$corridor=="1" & big.g$code=="18 M" & big.g$position.code=="patch")
dist.g18m1<-as.matrix(vegdist(Glomero.2018m.p1[,1:18],method="bray"))

gmay.18.1<-as.vector(rbind(dist.g18m1[1:3,4:6],dist.g18m1[12:14,15:17],t(dist.g18m1[7:9,10:11]),dist.g18m1[18:20,21:23],
                           dist.g18m1[24:26,27:29],dist.g18m1[30:31,32:34],dist.g18m1[35:37,38:40],dist.g18m1[41:43,44:46],
                           dist.g18m1[47:49,50:52],dist.g18m1[53:55,56:58]))

Glomero.2018m.p0<-subset(big.g, big.g$corridor=="0" & big.g$code=="18 M"& big.g$position.code=="patch" )
dist.g18m0<-as.matrix(vegdist(Glomero.2018m.p0[,1:18],method="bray"))

gmay.18.0<-as.vector(rbind(dist.g18m0[1:3,4:6],dist.g18m0[7:9,10:12],dist.g18m0[13:15,16:18],dist.g18m0[19:21,22:24],
                           dist.g18m0[25:27,28:30],dist.g18m0[31:32,33:35],dist.g18m0[36:38,39:41],
                           dist.g18m0[42:44,45:47],dist.g18m0[48:50,51:53]))

##### 2018 June
Glomero.2018j.p1<-subset(big.c, big.c$corridor=="1" & big.c$code=="18J" & big.c$position.code=="patch")
dist.g18j1<-as.matrix(vegdist(Glomero.2018j.p1[,1:18],method="bray"))

gjune.18.1<-as.vector(cbind(dist.g18j1[1:3,4:6],dist.g18j1[7:9,10:12],dist.g18j1[13:15,16:18],dist.g18j1[19:21,22:23],
                            dist.g18j1[24:26,27:29],dist.g18j1[35:37,38:40],t(dist.g18j1[30:31,32:34]),
                            dist.g18j1[41:43,44:46],dist.g18j1[47:49,50:52],dist.g18j1[53:55,56:58]))

Glomero.2018j.p0<-subset(big.g, big.g$corridor=="0" & big.g$code=="18J"& big.g$position.code=="patch" )
dist.g18j0<-as.matrix(vegdist(Glomero.2018j.p0[,1:18],method="bray"))

gjune.18.0<-as.vector(rbind(dist.g18j0[1:3,4:6],dist.g18j0[7:9,10:12],dist.g18j0[13:15,16:18],dist.g18j0[19:21,22:24],
                            dist.g18j0[25:27,28:30],dist.g18j0[31:32,33:35],dist.g18j0[36:38,39:41],
                            dist.g18j0[42:44,45:47],dist.g18j0[48:50,51:53]))

##### 2018 October
Glomero.2018o.p1<-subset(big.g, big.g$corridor=="1" & big.g$code=="18O" & big.g$position.code=="patch")
dist.g18o1<-as.matrix(vegdist(Glomero.2018o.p1[,1:18],method="bray"))

goct.18.1<-as.vector(cbind(dist.g18o1[1:3,4:6],dist.g18o1[7:9,10:11],dist.g18o1[12:14,15:17],dist.g18o1[18:20,21:23],
                           dist.g18o1[24:26,27:29],dist.g18o1[30:32,33:35],dist.g18o1[36:38,39:41],
                           dist.g18o1[42:44,45:46],dist.g18o1[47:49,50:52],dist.g18o1[53:55,56:58]))

Glomero.2018o.p0<-subset(big.g, big.g$corridor=="0" & big.g$code=="18O"& big.g$position.code=="patch" )
dist.g18o0<-as.matrix(vegdist(Glomero.2018o.p0[,1:18],method="bray"))

goct.18.0<-as.vector(rbind(dist.g18o0[1:3,4:6],dist.g18o0[7:8,9:11],dist.g18o0[12:14,15:17],dist.g18o0[18:20,21:23],
                           dist.g18o0[24:26,27:29],dist.g18o0[30:32,33:35],dist.g18o0[36:38,39:41],
                           dist.g18o0[42:44,45:47],dist.g18o0[48:50,51:53]))

##### 2019 May
Glomero.2019m.p1<-subset(big.g, big.g$corridor=="1" & big.g$code=="19M" & big.g$position.code=="patch")
dist.g19m1<-as.matrix(vegdist(Glomero.2019m.p1[,1:18],method="bray"))

gmay.19.1<-as.vector(cbind(dist.g19m1[1:3,4:6],dist.g19m1[7:9,10:12],dist.g19m1[13:15,16:18],dist.g19m1[22:24,25:27],
                           dist.g19m1[28:30,31:33],dist.g19m1[34:36,37:39],dist.g19m1[40:42,43:45],
                           dist.g19m1[46:48,49:51],dist.g19m1[52:54,55:56]))

Glomero.2019m.p0<-subset(big.g, big.g$corridor=="0" & big.g$code=="19M"& big.g$position.code=="patch" )
dist.g19m0<-as.matrix(vegdist(Glomero.2019m.p0[,1:18],method="bray"))

gmay.19.0<-as.vector(rbind(dist.g19m0[1:3,4:6],dist.g19m0[7:9,10:12],dist.g19m0[13:15,16:18],dist.g19m0[19:21,22:24],
                           dist.g19m0[25:27,28:30],dist.g19m0[31:33,34],dist.g19m0[35:37,38:40],dist.g19m0[41:43,44:46]))


dis.Glomero<-as.data.frame(cbind(goct.17.0,goct.17.1,gmay.18.0,gmay.18.1,gjune.18.0,gjune.18.1,goct.18.0,goct.18.1,gmay.19.0,gmay.19.1))
dis.Glomero1<-stack(dis.Glomero)
write.table(dis.Glomero,file="Bray curtis distance at individual level Glomero.txt",sep="\t")
dis.Glomero.i<-read.table("Bray curtis distance at individual level Glomero.txt", header=T, row.names=1, sep="\t")  ## otu abundance rawdata

####Zygomycota####
ngs.Zygo

## calculate the Bray-Curtis distance for each sampling campaign
ngs.zt<-as.data.frame(t(ngs.Zygo[,-1]))
big.z<-cbind(ngs.zt,env)

#### 2017 October
Zygo.2017.p1<-subset(big.c, big.c$corridor=="1" & big.c$code=="17O" & big.c$position.code=="patch")
dist.z17o1<-as.matrix(vegdist(Zygo.2017.p1[,1:30],method="bray"))

zoct.17.1<-as.vector(rbind(dist.z17o1[4:6,1:3],dist.z17o1[10:12,7:9],dist.z17o1[16:18,13:15],dist.z17o1[22:24,19:21],
                           dist.z17o1[28:30,25:27],dist.z17o1[34:36,31:33],dist.z17o1[40:41,37:39],
                           dist.z17o1[45:47,42:44],dist.z17o1[52:54,49:51],dist.z17o1[58,55:57]))

Zygo.2017.p0<-subset(big.z, big.z$corridor=="0" & big.z$code=="17O"& big.z$position.code=="patch" )
dist.z17o0<-as.matrix(vegdist(Zygo.2017.p0[,1:30],method="bray"))

zoct.17.0<-as.vector(rbind(dist.z17o0[4:6,1:3],dist.z17o0[10:12,7:9],dist.z17o0[16:18,13:15],dist.z17o0[22,19:21],
                           dist.z17o0[26:28,23:25],dist.z17o0[32:34,29:31],dist.z17o0[38:40,35:37],
                           dist.z17o0[44:46,41:43],dist.z17o0[50:52,47:49]))

##### 2018 May
Zygo.2018m.p1<-subset(big.z, big.z$corridor=="1" & big.z$code=="18 M" & big.z$position.code=="patch")
dist.z18m1<-as.matrix(vegdist(Zygo.2018m.p1[,1:30],method="bray"))

zmay.18.1<-as.vector(rbind(dist.z18m1[1:3,4:6],dist.z18m1[12:14,15:17],t(dist.z18m1[7:9,10:11]),dist.z18m1[18:20,21:23],
                           dist.z18m1[24:26,27:29],dist.z18m1[30:31,32:34],dist.z18m1[35:37,38:40],dist.z18m1[41:43,44:46],
                           dist.z18m1[47:49,50:52],dist.z18m1[53:55,56:58]))

Zygo.2018m.p0<-subset(big.z, big.z$corridor=="0" & big.z$code=="18 M"& big.z$position.code=="patch" )
dist.z18m0<-as.matrix(vegdist(Zygo.2018m.p0[,1:30],method="bray"))

zmay.18.0<-as.vector(rbind(dist.z18m0[1:3,4:6],dist.z18m0[7:9,10:12],dist.z18m0[13:15,16:18],dist.z18m0[19:21,22:24],
                           dist.z18m0[25:27,28:30],dist.z18m0[31:32,33:35],dist.z18m0[36:38,39:41],
                           dist.z18m0[42:44,45:47],dist.z18m0[48:50,51:53]))

##### 2018 June
Zygo.2018j.p1<-subset(big.z, big.z$corridor=="1" & big.z$code=="18J" & big.z$position.code=="patch")
dist.z18j1<-as.matrix(vegdist(Zygo.2018j.p1[,1:30],method="bray"))

zjune.18.1<-as.vector(cbind(dist.z18j1[1:3,4:6],dist.z18j1[7:9,10:12],dist.z18j1[13:15,16:18],dist.z18j1[19:21,22:23],
                            dist.z18j1[24:26,27:29],dist.z18j1[35:37,38:40],t(dist.z18j1[30:31,32:34]),
                            dist.z18j1[41:43,44:46],dist.z18j1[47:49,50:52],dist.z18j1[53:55,56:58]))

Zygo.2018j.p0<-subset(big.z, big.z$corridor=="0" & big.z$code=="18J"& big.z$position.code=="patch" )
dist.z18j0<-as.matrix(vegdist(Zygo.2018j.p0[,1:30],method="bray"))

zjune.18.0<-as.vector(rbind(dist.z18j0[1:3,4:6],dist.z18j0[7:9,10:12],dist.z18j0[13:15,16:18],dist.z18j0[19:21,22:24],
                            dist.z18j0[25:27,28:30],dist.z18j0[31:32,33:35],dist.z18j0[36:38,39:41],
                            dist.z18j0[42:44,45:47],dist.z18j0[48:50,51:53]))

##### 2018 October
Zygo.2018o.p1<-subset(big.z, big.z$corridor=="1" & big.z$code=="18O" & big.z$position.code=="patch")
dist.z18o1<-as.matrix(vegdist(Zygo.2018o.p1[,1:30],method="bray"))

zoct.18.1<-as.vector(cbind(dist.z18o1[1:3,4:6],dist.z18o1[7:9,10:11],dist.z18o1[12:14,15:17],dist.z18o1[18:20,21:23],
                           dist.z18o1[24:26,27:29],dist.z18o1[30:32,33:35],dist.z18o1[36:38,39:41],
                           dist.z18o1[42:44,45:46],dist.z18o1[47:49,50:52],dist.z18o1[53:55,56:58]))

Zygo.2018o.p0<-subset(big.z, big.z$corridor=="0" & big.z$code=="18O"& big.z$position.code=="patch" )
dist.z18o0<-as.matrix(vegdist(Zygo.2018o.p0[,1:30],method="bray"))

zoct.18.0<-as.vector(rbind(dist.z18o0[1:3,4:6],dist.z18o0[7:8,9:11],dist.z18o0[12:14,15:17],dist.z18o0[18:20,21:23],
                           dist.z18o0[24:26,27:29],dist.z18o0[30:32,33:35],dist.z18o0[36:38,39:41],
                           dist.z18o0[42:44,45:47],dist.z18o0[48:50,51:53]))

##### 2019 May
Zygo.2019m.p1<-subset(big.z, big.z$corridor=="1" & big.z$code=="19M" & big.z$position.code=="patch")
dist.z19m1<-as.matrix(vegdist(Zygo.2019m.p1[,1:30],method="bray"))

zmay.19.1<-as.vector(cbind(dist.z19m1[1:3,4:6],dist.z19m1[7:9,10:12],dist.z19m1[13:15,16:18],dist.z19m1[22:24,25:27],
                           dist.z19m1[28:30,31:33],dist.z19m1[34:36,37:39],dist.z19m1[40:42,43:45],
                           dist.z19m1[46:48,49:51],dist.z19m1[52:54,55:56]))

Zygo.2019m.p0<-subset(big.z, big.z$corridor=="0" & big.z$code=="19M"& big.z$position.code=="patch" )
dist.z19m0<-as.matrix(vegdist(Zygo.2019m.p0[,1:30],method="bray"))

zmay.19.0<-as.vector(rbind(dist.z19m0[1:3,4:6],dist.z19m0[7:9,10:12],dist.z19m0[13:15,16:18],dist.z19m0[19:21,22:24],
                           dist.z19m0[25:27,28:30],dist.z19m0[31:33,34],dist.z19m0[35:37,38:40],dist.z19m0[41:43,44:46]))


dis.Zygo<-as.data.frame(cbind(zoct.17.0,zoct.17.1,zmay.18.0,zmay.18.1,zjune.18.0,zjune.18.1,zoct.18.0,zoct.18.1,zmay.19.0,zmay.19.1))
dis.Zygo1<-stack(dis.Zygo)
write.table(dis.Zygo,file="Bray curtis distance at individual level Zygo.txt",sep="\t")
dis.Zygo.i<-read.table("Bray curtis distance at individual level Zygo.txt", header=T, row.names=1, sep="\t")  ## otu abundance rawdata


dis.Asco2<-dis.Asco1[-c(76:84,166:168,247:252,415:420,583:588,739:756,835:840),]
dis.Basidio2<-dis.Basidio1[-c(76:84,166:168,247:252,415:420,583:588,739:756,835:840),]
dis.Chytridio2<-dis.Chytridio1[-c(76:84,166:168,247:252,415:420,583:588,739:756,835:840),]
dis.Glomero2<-dis.Glomero1[-c(76:84,166:168,247:252,415:420,583:588,739:756,835:840),]
dis.Zygo2<-dis.Zygo1[-c(76:84,166:168,247:252,415:420,583:588,739:756,835:840),]

dis.phylum.all<-cbind(dis.Asco2,dis.Basidio2,dis.Chytridio2,dis.Glomero2,dis.Zygo2)
write.table(dis.phylum.all,file="Bray curtis distance at individual level All phylum.txt",sep="\t")
dis.phylum.all1<-read.table("Bray curtis distance at individual level All phylum1.txt", header=T, row.names=1, sep="\t")  ## otu abundance rawdata

####Bray-Curtis dissimilarity at patch level for each phylum root mycobiota ####
#### Distance among population 1 and 2 ####
ngs.poolf<-cbind(richness.poolf[,1:6],ngs.pool1)
ngs.poolpp<-subset(ngs.poolf,ngs.poolf$position!="cor"&ngs.poolf$position!="mb")
env.p<-subset(env.pool,env.pool$Group.5!="cor"&env.pool$Group.5!="mb")
  
#### Ascomycota ####
ngs.Asco.p<-subset(ngs.poolp.phylum,ngs.poolp.phylum$`taxo$blast_perc_identity`=="Ascomycota")
ngs.Asco.p1<-cbind(env.p,t(ngs.Asco.p[,-1]))

#### 2017 October ####
Asco.p.2017<-subset(ngs.Asco.p1,ngs.Asco.p1$Group.4=="17O")
Asco.p.2017.1<-subset(Asco.p.2017,Asco.p.2017$Group.6==1)
Asco.p.2017.0<-subset(Asco.p.2017,Asco.p.2017$Group.6==0)

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
apdist.17o<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(Asco.p.2017.1)[1]-1)
{
  p1<-as.vector(t(Asco.p.2017.1[i,-c(1:7)]))
  p2<-as.vector(t(Asco.p.2017.1[i+1,-c(1:7)]))
  apdist.17o[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

for (i in 1:dim(Asco.p.2017.0)[1]-1)
{
  p1<-as.vector(t(Asco.p.2017.0[i,-c(1:7)]))
  p2<-as.vector(t(Asco.p.2017.0[i+1,-c(1:7)]))
  apdist.17o[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

apdist.17of<-apdist.17o[ c(TRUE,FALSE), ]  # rows
colnames(apdist.17of)<-c("Asco 2017 October With corridor","Asco 2017 October Without corridor")

#### 2018 May ####
Asco.p.18m<-subset(ngs.Asco.p1,ngs.Asco.p1$Group.4=="18 M")
Asco.p.18m.1<-subset(Asco.p.18m,Asco.p.18m$Group.6==1)
Asco.p.18m.0<-subset(Asco.p.18m,Asco.p.18m$Group.6==0)

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
apdist.18m<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(Asco.p.18m.1)[1]-1)
{
  p1<-as.vector(t(Asco.p.18m.1[i,-c(1:7)]))
  p2<-as.vector(t(Asco.p.18m.1[i+1,-c(1:7)]))
  apdist.18m[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

for (i in 1:dim(Asco.p.18m.0)[1]-1)
{
  p1<-as.vector(t(Asco.p.18m.0[i,-c(1:7)]))
  p2<-as.vector(t(Asco.p.18m.0[i+1,-c(1:7)]))
  apdist.18m[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

apdist.18mf<-apdist.18m[ c(TRUE,FALSE), ]  # rows
colnames(apdist.18mf)<-c("Asco 2018 May With corridor","Asco 2018 May Without corridor")

#### 2018 June ####
Asco.p.18j<-subset(ngs.Asco.p1,ngs.Asco.p1$Group.4=="18J")
Asco.p.18j.1<-subset(Asco.p.18j,Asco.p.18j$Group.6==1)
Asco.p.18j.0<-subset(Asco.p.18j,Asco.p.18j$Group.6==0)

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
apdist.18j<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(Asco.p.18j.1)[1]-1)
{
  p1<-as.vector(t(Asco.p.18j.1[i,-c(1:7)]))
  p2<-as.vector(t(Asco.p.18j.1[i+1,-c(1:7)]))
  apdist.18j[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

for (i in 1:dim(Asco.p.18j.0)[1]-1)
{
  p1<-as.vector(t(Asco.p.18j.0[i,-c(1:7)]))
  p2<-as.vector(t(Asco.p.18j.0[i+1,-c(1:7)]))
  apdist.18j[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

apdist.18jf<-apdist.18j[ c(TRUE,FALSE), ]  # rows
colnames(apdist.18jf)<-c("Asco 2018 June With corridor","Asco 2018 June Without corridor")

#### 2018 October ####
Asco.p.18o<-subset(ngs.Asco.p1,ngs.Asco.p1$Group.4=="18O")
Asco.p.18o.1<-subset(Asco.p.18o,Asco.p.18o$Group.6==1)
Asco.p.18o.0<-subset(Asco.p.18o,Asco.p.18o$Group.6==0)

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
apdist.18o<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(Asco.p.18o.1)[1]-1)
{
  p1<-as.vector(t(Asco.p.18o.1[i,-c(1:7)]))
  p2<-as.vector(t(Asco.p.18o.1[i+1,-c(1:7)]))
  apdist.18o[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

for (i in 1:dim(Asco.p.18o.0)[1]-1)
{
  p1<-as.vector(t(Asco.p.18o.0[i,-c(1:7)]))
  p2<-as.vector(t(Asco.p.18o.0[i+1,-c(1:7)]))
  apdist.18o[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

apdist.18of<-apdist.18o[ c(TRUE,FALSE), ]  # rows
colnames(apdist.18of)<-c("Asco 2018 October With corridor","Asco 2018 October Without corridor")

#### 2019 May ####
Asco.p.19m<-subset(ngs.Asco.p1,ngs.Asco.p1$Group.4=="19M")
Asco.p.19m.1<-subset(Asco.p.19m,Asco.p.19m$Group.6==1)
Asco.p.19m.0<-subset(Asco.p.19m,Asco.p.19m$Group.6==0)

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
apdist.19m<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(Asco.p.19m.1)[1]-1)
{
  p1<-as.vector(t(Asco.p.19m.1[i,-c(1:7)]))
  p2<-as.vector(t(Asco.p.19m.1[i+1,-c(1:7)]))
  apdist.19m[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

for (i in 1:dim(Asco.p.19m.0)[1]-1)
{
  p1<-as.vector(t(Asco.p.19m.0[i,-c(1:7)]))
  p2<-as.vector(t(Asco.p.19m.0[i+1,-c(1:7)]))
  apdist.19m[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

apdist.19mf<-apdist.19m[ c(TRUE,FALSE), ]  # rows
colnames(apdist.19mf)<-c("Asco 2019 May With corridor","Asco 2019 May Without corridor")

asco.p.da<-cbind(apdist.17of,apdist.18mf,apdist.18jf,apdist.18of,apdist.19mf)

asco.p.daf<-stack(asco.p.da)
boxplot(values~ind,asco.p.daf)

write.table(asco.p.daf,file="asco.p.daf.txt",sep="\t")

#dist.alle<-read.table("dist.all20210421.txt", sep="\t", header=T)        ## sample properties, we need this

#### Basidiomycota ####
ngs.Basidio.p<-subset(ngs.poolp.phylum,ngs.poolp.phylum$`taxo$blast_perc_identity`=="Basidiomycota")
ngs.Basidio.p1<-cbind(env.p,t(ngs.Basidio.p[,-1]))

#### 2017 October ####
Basidio.p.2017<-subset(ngs.Basidio.p1,ngs.Basidio.p1$Group.4=="17O")
Basidio.p.2017.1<-subset(Basidio.p.2017,Basidio.p.2017$Group.6==1)
Basidio.p.2017.0<-subset(Basidio.p.2017,Basidio.p.2017$Group.6==0)

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
bpdist.17o<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(Basidio.p.2017.1)[1]-1)
{
  p1<-as.vector(t(Basidio.p.2017.1[i,-c(1:7)]))
  p2<-as.vector(t(Basidio.p.2017.1[i+1,-c(1:7)]))
  bpdist.17o[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

for (i in 1:dim(Basidio.p.2017.0)[1]-1)
{
  p1<-as.vector(t(Basidio.p.2017.0[i,-c(1:7)]))
  p2<-as.vector(t(Basidio.p.2017.0[i+1,-c(1:7)]))
  bpdist.17o[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

bpdist.17of<-bpdist.17o[ c(TRUE,FALSE), ]  # rows
colnames(bpdist.17of)<-c("Basidio 2017 October With corridor","Basidio 2017 October Without corridor")

#### 2018 May ####
Basidio.p.18m<-subset(ngs.Basidio.p1,ngs.Basidio.p1$Group.4=="18 M")
Basidio.p.18m.1<-subset(Basidio.p.18m,Basidio.p.18m$Group.6==1)
Basidio.p.18m.0<-subset(Basidio.p.18m,Basidio.p.18m$Group.6==0)

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
bpdist.18m<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(Basidio.p.18m.1)[1]-1)
{
  p1<-as.vector(t(Basidio.p.18m.1[i,-c(1:7)]))
  p2<-as.vector(t(Basidio.p.18m.1[i+1,-c(1:7)]))
  bpdist.18m[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

for (i in 1:dim(Basidio.p.18m.0)[1]-1)
{
  p1<-as.vector(t(Basidio.p.18m.0[i,-c(1:7)]))
  p2<-as.vector(t(Basidio.p.18m.0[i+1,-c(1:7)]))
  bpdist.18m[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

bpdist.18mf<-bpdist.18m[ c(TRUE,FALSE), ]  # rows
colnames(bpdist.18mf)<-c("Basidio 2018 May With corridor","Basidio 2018 May Without corridor")

#### 2018 June ####
Basidio.p.18j<-subset(ngs.Basidio.p1,ngs.Basidio.p1$Group.4=="18J")
Basidio.p.18j.1<-subset(Basidio.p.18j,Basidio.p.18j$Group.6==1)
Basidio.p.18j.0<-subset(Basidio.p.18j,Basidio.p.18j$Group.6==0)

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
bpdist.18j<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(Basidio.p.18j.1)[1]-1)
{
  p1<-as.vector(t(Basidio.p.18j.1[i,-c(1:7)]))
  p2<-as.vector(t(Basidio.p.18j.1[i+1,-c(1:7)]))
  bpdist.18j[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

for (i in 1:dim(Basidio.p.18j.0)[1]-1)
{
  p1<-as.vector(t(Basidio.p.18j.0[i,-c(1:7)]))
  p2<-as.vector(t(Basidio.p.18j.0[i+1,-c(1:7)]))
  bpdist.18j[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

bpdist.18jf<-bpdist.18j[ c(TRUE,FALSE), ]  # rows
colnames(bpdist.18jf)<-c("Basidio 2018 June With corridor","Basidio 2018 June Without corridor")

#### 2018 October ####
Basidio.p.18o<-subset(ngs.Basidio.p1,ngs.Basidio.p1$Group.4=="18O")
Basidio.p.18o.1<-subset(Basidio.p.18o,Basidio.p.18o$Group.6==1)
Basidio.p.18o.0<-subset(Basidio.p.18o,Basidio.p.18o$Group.6==0)

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
bpdist.18o<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(Basidio.p.18o.1)[1]-1)
{
  p1<-as.vector(t(Basidio.p.18o.1[i,-c(1:7)]))
  p2<-as.vector(t(Basidio.p.18o.1[i+1,-c(1:7)]))
  bpdist.18o[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

for (i in 1:dim(Basidio.p.18o.0)[1]-1)
{
  p1<-as.vector(t(Basidio.p.18o.0[i,-c(1:7)]))
  p2<-as.vector(t(Basidio.p.18o.0[i+1,-c(1:7)]))
  bpdist.18o[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

bpdist.18of<-bpdist.18o[ c(TRUE,FALSE), ]  # rows
colnames(bpdist.18of)<-c("Basidio 2018 October With corridor","Basidio 2018 October Without corridor")

#### 2019 May ####
Basidio.p.19m<-subset(ngs.Basidio.p1,ngs.Basidio.p1$Group.4=="19M")
Basidio.p.19m.1<-subset(Basidio.p.19m,Basidio.p.19m$Group.6==1)
Basidio.p.19m.0<-subset(Basidio.p.19m,Basidio.p.19m$Group.6==0)

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
bpdist.19m<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(Basidio.p.19m.1)[1]-1)
{
  p1<-as.vector(t(Basidio.p.19m.1[i,-c(1:7)]))
  p2<-as.vector(t(Basidio.p.19m.1[i+1,-c(1:7)]))
  bpdist.19m[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

for (i in 1:dim(Basidio.p.19m.0)[1]-1)
{
  p1<-as.vector(t(Basidio.p.19m.0[i,-c(1:7)]))
  p2<-as.vector(t(Basidio.p.19m.0[i+1,-c(1:7)]))
  bpdist.19m[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

bpdist.19mf<-bpdist.19m[ c(TRUE,FALSE), ]  # rows
colnames(bpdist.19mf)<-c("Basidio 2019 May With corridor","Basidio 2019 May Without corridor")

Basidio.p.da<-cbind(bpdist.17of,bpdist.18mf,bpdist.18jf,bpdist.18of,bpdist.19mf)

Basidio.p.daf<-stack(Basidio.p.da)
boxplot(values~ind,Basidio.p.daf)

write.table(Basidio.p.daf,file="Basidio.p.daf.txt",sep="\t")

#dist.alle<-read.table("dist.all20210421.txt", sep="\t", header=T)        ## sample properties, we need this

#### Chytridiomycota ####
ngs.Chytridio.p<-subset(ngs.poolp.phylum,ngs.poolp.phylum$`taxo$blast_perc_identity`=="Chytridiomycota")
ngs.Chytridio.p1<-cbind(env.p,t(ngs.Chytridio.p[,-1]))

#### 2017 October ####
Chytridio.p.2017<-subset(ngs.Chytridio.p1,ngs.Chytridio.p1$Group.4=="17O")
Chytridio.p.2017.1<-subset(Chytridio.p.2017,Chytridio.p.2017$Group.6==1)
Chytridio.p.2017.0<-subset(Chytridio.p.2017,Chytridio.p.2017$Group.6==0)

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
cpdist.17o<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(Chytridio.p.2017.1)[1]-1)
{
  p1<-as.vector(t(Chytridio.p.2017.1[i,-c(1:7)]))
  p2<-as.vector(t(Chytridio.p.2017.1[i+1,-c(1:7)]))
  cpdist.17o[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

for (i in 1:dim(Chytridio.p.2017.0)[1]-1)
{
  p1<-as.vector(t(Chytridio.p.2017.0[i,-c(1:7)]))
  p2<-as.vector(t(Chytridio.p.2017.0[i+1,-c(1:7)]))
  cpdist.17o[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

cpdist.17of<-cpdist.17o[ c(TRUE,FALSE), ]  # rows
colnames(cpdist.17of)<-c("Chytridio 2017 October With corridor","Chytridio 2017 October Without corridor")

#### 2018 May ####
Chytridio.p.18m<-subset(ngs.Chytridio.p1,ngs.Chytridio.p1$Group.4=="18 M")
Chytridio.p.18m.1<-subset(Chytridio.p.18m,Chytridio.p.18m$Group.6==1)
Chytridio.p.18m.0<-subset(Chytridio.p.18m,Chytridio.p.18m$Group.6==0)

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
cpdist.18m<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(Chytridio.p.18m.1)[1]-1)
{
  p1<-as.vector(t(Chytridio.p.18m.1[i,-c(1:7)]))
  p2<-as.vector(t(Chytridio.p.18m.1[i+1,-c(1:7)]))
  cpdist.18m[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

for (i in 1:dim(Chytridio.p.18m.0)[1]-1)
{
  p1<-as.vector(t(Chytridio.p.18m.0[i,-c(1:7)]))
  p2<-as.vector(t(Chytridio.p.18m.0[i+1,-c(1:7)]))
  cpdist.18m[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

cpdist.18mf<-cpdist.18m[ c(TRUE,FALSE), ]  # rows
colnames(cpdist.18mf)<-c("Chytridio 2018 May With corridor","Chytridio 2018 May Without corridor")

#### 2018 June ####
Chytridio.p.18j<-subset(ngs.Chytridio.p1,ngs.Chytridio.p1$Group.4=="18J")
Chytridio.p.18j.1<-subset(Chytridio.p.18j,Chytridio.p.18j$Group.6==1)
Chytridio.p.18j.0<-subset(Chytridio.p.18j,Chytridio.p.18j$Group.6==0)

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
cpdist.18j<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(Chytridio.p.18j.1)[1]-1)
{
  p1<-as.vector(t(Chytridio.p.18j.1[i,-c(1:7)]))
  p2<-as.vector(t(Chytridio.p.18j.1[i+1,-c(1:7)]))
  cpdist.18j[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

for (i in 1:dim(Chytridio.p.18j.0)[1]-1)
{
  p1<-as.vector(t(Chytridio.p.18j.0[i,-c(1:7)]))
  p2<-as.vector(t(Chytridio.p.18j.0[i+1,-c(1:7)]))
  cpdist.18j[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

cpdist.18jf<-cpdist.18j[ c(TRUE,FALSE), ]  # rows
colnames(cpdist.18jf)<-c("Chytridio 2018 June With corridor","Chytridio 2018 June Without corridor")

#### 2018 October ####
Chytridio.p.18o<-subset(ngs.Chytridio.p1,ngs.Chytridio.p1$Group.4=="18O")
Chytridio.p.18o.1<-subset(Chytridio.p.18o,Chytridio.p.18o$Group.6==1)
Chytridio.p.18o.0<-subset(Chytridio.p.18o,Chytridio.p.18o$Group.6==0)

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
cpdist.18o<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(Chytridio.p.18o.1)[1]-1)
{
  p1<-as.vector(t(Chytridio.p.18o.1[i,-c(1:7)]))
  p2<-as.vector(t(Chytridio.p.18o.1[i+1,-c(1:7)]))
  cpdist.18o[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

for (i in 1:dim(Chytridio.p.18o.0)[1]-1)
{
  p1<-as.vector(t(Chytridio.p.18o.0[i,-c(1:7)]))
  p2<-as.vector(t(Chytridio.p.18o.0[i+1,-c(1:7)]))
  cpdist.18o[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

cpdist.18of<-cpdist.18o[ c(TRUE,FALSE), ]  # rows
colnames(cpdist.18of)<-c("Chytridio 2018 October With corridor","Chytridio 2018 October Without corridor")

#### 2019 May ####
Chytridio.p.19m<-subset(ngs.Chytridio.p1,ngs.Chytridio.p1$Group.4=="19M")
Chytridio.p.19m.1<-subset(Chytridio.p.19m,Chytridio.p.19m$Group.6==1)
Chytridio.p.19m.0<-subset(Chytridio.p.19m,Chytridio.p.19m$Group.6==0)

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
cpdist.19m<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(Chytridio.p.19m.1)[1]-1)
{
  p1<-as.vector(t(Chytridio.p.19m.1[i,-c(1:7)]))
  p2<-as.vector(t(Chytridio.p.19m.1[i+1,-c(1:7)]))
  cpdist.19m[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

for (i in 1:dim(Chytridio.p.19m.0)[1]-1)
{
  p1<-as.vector(t(Chytridio.p.19m.0[i,-c(1:7)]))
  p2<-as.vector(t(Chytridio.p.19m.0[i+1,-c(1:7)]))
  cpdist.19m[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

cpdist.19mf<-cpdist.19m[ c(TRUE,FALSE), ]  # rows
colnames(cpdist.19mf)<-c("Chytridio 2019 May With corridor","Chytridio 2019 May Without corridor")

Chytridio.p.da<-cbind(cpdist.17of,cpdist.18mf,cpdist.18jf,cpdist.18of,cpdist.19mf)

Chytridio.p.daf<-stack(Chytridio.p.da)
boxplot(values~ind,Chytridio.p.daf)

write.table(Chytridio.p.daf,file="Chytridio.p.daf.txt",sep="\t")

#dist.alle<-read.table("dist.all20210421.txt", sep="\t", header=T)        ## sample properties, we need this


#### Glomeromycota ####
ngs.Glomero.p<-subset(ngs.poolp.phylum,ngs.poolp.phylum$`taxo$blast_perc_identity`=="Glomeromycota")
ngs.Glomero.p1<-cbind(env.p,t(ngs.Glomero.p[,-1]))

#### 2017 October ####
Glomero.p.2017<-subset(ngs.Glomero.p1,ngs.Glomero.p1$Group.4=="17O")
Glomero.p.2017.1<-subset(Glomero.p.2017,Glomero.p.2017$Group.6==1)
Glomero.p.2017.0<-subset(Glomero.p.2017,Glomero.p.2017$Group.6==0)

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
gpdist.17o<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(Glomero.p.2017.1)[1]-1)
{
  p1<-as.vector(t(Glomero.p.2017.1[i,-c(1:7)]))
  p2<-as.vector(t(Glomero.p.2017.1[i+1,-c(1:7)]))
  gpdist.17o[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

for (i in 1:dim(Glomero.p.2017.0)[1]-1)
{
  p1<-as.vector(t(Glomero.p.2017.0[i,-c(1:7)]))
  p2<-as.vector(t(Glomero.p.2017.0[i+1,-c(1:7)]))
  gpdist.17o[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

gpdist.17of<-gpdist.17o[ c(TRUE,FALSE), ]  # rows
colnames(gpdist.17of)<-c("Glomero 2017 October With corridor","Glomero 2017 October Without corridor")

#### 2018 May ####
Glomero.p.18m<-subset(ngs.Glomero.p1,ngs.Glomero.p1$Group.4=="18 M")
Glomero.p.18m.1<-subset(Glomero.p.18m,Glomero.p.18m$Group.6==1)
Glomero.p.18m.0<-subset(Glomero.p.18m,Glomero.p.18m$Group.6==0)

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
gpdist.18m<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(Glomero.p.18m.1)[1]-1)
{
  p1<-as.vector(t(Glomero.p.18m.1[i,-c(1:7)]))
  p2<-as.vector(t(Glomero.p.18m.1[i+1,-c(1:7)]))
  gpdist.18m[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

for (i in 1:dim(Glomero.p.18m.0)[1]-1)
{
  p1<-as.vector(t(Glomero.p.18m.0[i,-c(1:7)]))
  p2<-as.vector(t(Glomero.p.18m.0[i+1,-c(1:7)]))
  gpdist.18m[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

gpdist.18mf<-gpdist.18m[ c(TRUE,FALSE), ]  # rows
colnames(gpdist.18mf)<-c("Glomero 2018 May With corridor","Glomero 2018 May Without corridor")

#### 2018 June ####
Glomero.p.18j<-subset(ngs.Glomero.p1,ngs.Glomero.p1$Group.4=="18J")
Glomero.p.18j.1<-subset(Glomero.p.18j,Glomero.p.18j$Group.6==1)
Glomero.p.18j.0<-subset(Glomero.p.18j,Glomero.p.18j$Group.6==0)

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
gpdist.18j<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(Glomero.p.18j.1)[1]-1)
{
  p1<-as.vector(t(Glomero.p.18j.1[i,-c(1:7)]))
  p2<-as.vector(t(Glomero.p.18j.1[i+1,-c(1:7)]))
  gpdist.18j[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

for (i in 1:dim(Glomero.p.18j.0)[1]-1)
{
  p1<-as.vector(t(Glomero.p.18j.0[i,-c(1:7)]))
  p2<-as.vector(t(Glomero.p.18j.0[i+1,-c(1:7)]))
  gpdist.18j[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

gpdist.18jf<-gpdist.18j[ c(TRUE,FALSE), ]  # rows
colnames(gpdist.18jf)<-c("Glomero 2018 June With corridor","Glomero 2018 June Without corridor")

#### 2018 October ####
Glomero.p.18o<-subset(ngs.Glomero.p1,ngs.Glomero.p1$Group.4=="18O")
Glomero.p.18o.1<-subset(Glomero.p.18o,Glomero.p.18o$Group.6==1)
Glomero.p.18o.0<-subset(Glomero.p.18o,Glomero.p.18o$Group.6==0)

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
gpdist.18o<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(Glomero.p.18o.1)[1]-1)
{
  p1<-as.vector(t(Glomero.p.18o.1[i,-c(1:7)]))
  p2<-as.vector(t(Glomero.p.18o.1[i+1,-c(1:7)]))
  gpdist.18o[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

for (i in 1:dim(Glomero.p.18o.0)[1]-1)
{
  p1<-as.vector(t(Glomero.p.18o.0[i,-c(1:7)]))
  p2<-as.vector(t(Glomero.p.18o.0[i+1,-c(1:7)]))
  gpdist.18o[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

gpdist.18of<-gpdist.18o[ c(TRUE,FALSE), ]  # rows
colnames(gpdist.18of)<-c("Glomero 2018 October With corridor","Glomero 2018 October Without corridor")

#### 2019 May ####
Glomero.p.19m<-subset(ngs.Glomero.p1,ngs.Glomero.p1$Group.4=="19M")
Glomero.p.19m.1<-subset(Glomero.p.19m,Glomero.p.19m$Group.6==1)
Glomero.p.19m.0<-subset(Glomero.p.19m,Glomero.p.19m$Group.6==0)

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
gpdist.19m<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(Glomero.p.19m.1)[1]-1)
{
  p1<-as.vector(t(Glomero.p.19m.1[i,-c(1:7)]))
  p2<-as.vector(t(Glomero.p.19m.1[i+1,-c(1:7)]))
  gpdist.19m[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

for (i in 1:dim(Glomero.p.19m.0)[1]-1)
{
  p1<-as.vector(t(Glomero.p.19m.0[i,-c(1:7)]))
  p2<-as.vector(t(Glomero.p.19m.0[i+1,-c(1:7)]))
  gpdist.19m[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

gpdist.19mf<-apdist.19m[ c(TRUE,FALSE), ]  # rows
colnames(gpdist.19mf)<-c("Glomero 2019 May With corridor","Glomero 2019 May Without corridor")

Glomero.p.da<-cbind(gpdist.17of,gpdist.18mf,gpdist.18jf,gpdist.18of,gpdist.19mf)

Glomero.p.daf<-stack(Glomero.p.da)
boxplot(values~ind,Glomero.p.daf)

write.table(Glomero.p.daf,file="Glomero.p.daf.txt",sep="\t")

#dist.alle<-read.table("dist.all20210421.txt", sep="\t", header=T)        ## sample properties, we need this

#### Zygomycota ####
ngs.Zygo.p<-subset(ngs.poolp.phylum,ngs.poolp.phylum$`taxo$blast_perc_identity`=="Zygomycota")
ngs.Zygo.p1<-cbind(env.p,t(ngs.Zygo.p[,-1]))

#### 2017 October ####
Zygo.p.2017<-subset(ngs.Zygo.p1,ngs.Zygo.p1$Group.4=="17O")
Zygo.p.2017.1<-subset(Zygo.p.2017,Zygo.p.2017$Group.6==1)
Zygo.p.2017.0<-subset(Zygo.p.2017,Zygo.p.2017$Group.6==0)

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
zpdist.17o<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(Zygo.p.2017.1)[1]-1)
{
  p1<-as.vector(t(Zygo.p.2017.1[i,-c(1:7)]))
  p2<-as.vector(t(Zygo.p.2017.1[i+1,-c(1:7)]))
  zpdist.17o[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

for (i in 1:dim(Zygo.p.2017.0)[1]-1)
{
  p1<-as.vector(t(Zygo.p.2017.0[i,-c(1:7)]))
  p2<-as.vector(t(Zygo.p.2017.0[i+1,-c(1:7)]))
  zpdist.17o[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

zpdist.17of<-zpdist.17o[ c(TRUE,FALSE), ]  # rows
colnames(zpdist.17of)<-c("Zygo 2017 October With corridor","Zygo 2017 October Without corridor")

#### 2018 May ####
Zygo.p.18m<-subset(ngs.Zygo.p1,ngs.Zygo.p1$Group.4=="18 M")
Zygo.p.18m.1<-subset(Zygo.p.18m,Zygo.p.18m$Group.6==1)
Zygo.p.18m.0<-subset(Zygo.p.18m,Zygo.p.18m$Group.6==0)

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
zpdist.18m<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(Zygo.p.18m.1)[1]-1)
{
  p1<-as.vector(t(Zygo.p.18m.1[i,-c(1:7)]))
  p2<-as.vector(t(Zygo.p.18m.1[i+1,-c(1:7)]))
  zpdist.18m[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

for (i in 1:dim(Zygo.p.18m.0)[1]-1)
{
  p1<-as.vector(t(Zygo.p.18m.0[i,-c(1:7)]))
  p2<-as.vector(t(Zygo.p.18m.0[i+1,-c(1:7)]))
  zpdist.18m[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

zpdist.18mf<-zpdist.18m[ c(TRUE,FALSE), ]  # rows
colnames(zpdist.18mf)<-c("Zygo 2018 May With corridor","Zygo 2018 May Without corridor")

#### 2018 June ####
Zygo.p.18j<-subset(ngs.Zygo.p1,ngs.Zygo.p1$Group.4=="18J")
Zygo.p.18j.1<-subset(Zygo.p.18j,Zygo.p.18j$Group.6==1)
Zygo.p.18j.0<-subset(Zygo.p.18j,Zygo.p.18j$Group.6==0)

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
zpdist.18j<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(Zygo.p.18j.1)[1]-1)
{
  p1<-as.vector(t(Zygo.p.18j.1[i,-c(1:7)]))
  p2<-as.vector(t(Zygo.p.18j.1[i+1,-c(1:7)]))
  zpdist.18j[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

for (i in 1:dim(Zygo.p.18j.0)[1]-1)
{
  p1<-as.vector(t(Zygo.p.18j.0[i,-c(1:7)]))
  p2<-as.vector(t(Zygo.p.18j.0[i+1,-c(1:7)]))
  zpdist.18j[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

zpdist.18jf<-zpdist.18j[ c(TRUE,FALSE), ]  # rows
colnames(zpdist.18jf)<-c("Zygo 2018 June With corridor","Zygo 2018 June Without corridor")

#### 2018 October ####
Zygo.p.18o<-subset(ngs.Zygo.p1,ngs.Zygo.p1$Group.4=="18O")
Zygo.p.18o.1<-subset(Zygo.p.18o,Zygo.p.18o$Group.6==1)
Zygo.p.18o.0<-subset(Zygo.p.18o,Zygo.p.18o$Group.6==0)

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
zpdist.18o<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(Zygo.p.18o.1)[1]-1)
{
  p1<-as.vector(t(Zygo.p.18o.1[i,-c(1:7)]))
  p2<-as.vector(t(Zygo.p.18o.1[i+1,-c(1:7)]))
  zpdist.18o[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

for (i in 1:dim(Zygo.p.18o.0)[1]-1)
{
  p1<-as.vector(t(Zygo.p.18o.0[i,-c(1:7)]))
  p2<-as.vector(t(Zygo.p.18o.0[i+1,-c(1:7)]))
  zpdist.18o[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

zpdist.18of<-zpdist.18o[ c(TRUE,FALSE), ]  # rows
colnames(zpdist.18of)<-c("Zygo 2018 October With corridor","Zygo 2018 October Without corridor")

#### 2019 May ####
Zygo.p.19m<-subset(ngs.Zygo.p1,ngs.Zygo.p1$Group.4=="19M")
Zygo.p.19m.1<-subset(Zygo.p.19m,Zygo.p.19m$Group.6==1)
Zygo.p.19m.0<-subset(Zygo.p.19m,Zygo.p.19m$Group.6==0)

## here is the loop to calculate all the point, CAUTION 1,2;3,4;5,6; NOT 1,2;2,3;3,4;
zpdist.19m<-as.data.frame(matrix(NA,nrow=20,ncol=2))
for (i in 1:dim(Zygo.p.19m.1)[1]-1)
{
  p1<-as.vector(t(Zygo.p.19m.1[i,-c(1:7)]))
  p2<-as.vector(t(Zygo.p.19m.1[i+1,-c(1:7)]))
  zpdist.19m[i,1]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

for (i in 1:dim(Zygo.p.19m.0)[1]-1)
{
  p1<-as.vector(t(Zygo.p.19m.0[i,-c(1:7)]))
  p2<-as.vector(t(Zygo.p.19m.0[i+1,-c(1:7)]))
  zpdist.19m[i,2]<-vegdist(rbind(p1,p2),method="bray",na.rm=TRUE)
}

zpdist.19mf<-zpdist.19m[ c(TRUE,FALSE), ]  # rows
colnames(zpdist.19mf)<-c("Zygo 2019 May With corridor","Zygo 2019 May Without corridor")

Zygo.p.da<-cbind(zpdist.17of,zpdist.18mf,zpdist.18jf,zpdist.18of,zpdist.19mf)

Zygo.p.daf<-stack(Zygo.p.da)
boxplot(values~ind,Zygo.p.daf)

write.table(Zygo.p.daf,file="Zygo.p.daf.txt",sep="\t")

dis.phylum.allp<-cbind(asco.p.daf,Basidio.p.daf,Chytridio.p.daf,Glomero.p.daf,Zygo.p.daf)
write.table(dis.phylum.allp,file="Bray curtis distance at patch level All phylum.txt",sep="\t")
dis.phylum.allp1<-read.table("Bray curtis distance at patch level All phylum1.txt", header=T, row.names=1, sep="\t")  ## otu abundance rawdata






