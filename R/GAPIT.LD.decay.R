`GAPIT.LD.decay` <-function(GI=NULL,X=NULL,chr=NULL, cut.dis=1,n.select=10000,
                                 WS0=NULL,
                                 ws=100,
                                 Aver.Dis=1000,
                                 max.num=10,  ## set the number of max selected range
                                 fre.by=100,  ## set 
                                 MAXfregment=NULL,
                                 max.number=NULL){
# Object: Analysis for Genotype data:Distribution of SNP density,Accumulation,Moving Average of density,result:a pdf of the scree plot
# myG:Genotype data
# chr: chromosome value
# WS0 is the cutoff threshold for marker to display
# ws is used to calculate within windowsize
# Aver.Dis is average display windowsize
# mav1:Moving Average set value length
# Authors:  Zhiwu Zhang and Jiabo Wang
# Last update: AUG 24, 2022 
##############################################################################################

# GI=myGM
# X=myGD[,-1]
# nchar=max.num
map.store=NULL
for(i in 1:length(unique(GI[,2])))
{
  map0=GI[GI[,2]==i,]
  map.store=append(map.store,nrow(map0))
}
max.number=floor(min(map.store)*0.7)
print(paste("There is ",max.number," markers in the shortest chromosome.",sep=""))
# print(max.number)
max.dist=max(GI[,3])
if(is.null(MAXfregment))MAXfregment=10^nchar(max(GI[,3]))

posi=as.numeric(GI[,2])*10^(nchar(max(as.numeric(GI[,3]))))+as.numeric(GI[,3])
set.seed(99163)

rs.index=sample(nrow(GI),n.select)
rs.index=sort(rs.index)
# print(table(rs.index))

Xr=X[,rs.index]
myGMr=myGM[rs.index,]
posi.rs=abs(as.numeric(myGMr[,2])*10^(nchar(max(as.numeric(myGMr[,3]))))+as.numeric(myGMr[,3]))
# if(is.null(WS0)) WS0=(max(posi.rs)%/%1000)*1000
# if(WS0==0)WS0=1
# posi.rs[posi.rs>max.dist]=NA
n=length(posi.rs)
if(n.select>n)n.select=n
if(is.null(max.number)) max.number=nrow(myGMr)/max.num
if(max.number>n.select) max.number=n.select
max.number=100
freg=floor(seq(1,max.number,length.out=7))
freg=c(1,3,7,15,31,63,127)
freg.legend=c(2,4,8,16,32,64,128)
# print(freg)
data.all=list()
dis.all=NULL
R.all=NULL
for(i in 1:7)
{
  dist=abs(posi.rs[-c(1:freg[i])]-posi.rs[-c((n-freg[i]+1):n)])
  # dist.out=GAPIT.Remove.outliers(dist,pro=0.1,size=1.1)
  # print(i)
  index.na=dist>max.dist
  x1=Xr[,-c((n-freg[i]+1):n)]
  x2=Xr[,-c(1:freg[i])]
  R.x=as.numeric(mapply(GAPIT.Cor.matrix,as.data.frame(x1),as.data.frame(x2)))
  R.x[is.na(R.x)]=0
  data.all[[i]]=cbind(dist[!index.na],R.x[!index.na])
  dis.all=append(dis.all,as.numeric(data.all[[i]][,1]))
  R.all=append(R.all,as.numeric(data.all[[i]][,2]))
}


if(is.null(WS0)) WS0=((max(dis.all))%/%1000)*1000
if(WS0>max.dist) WS0=max.dist
# WS0=100000000000
grDevices::pdf("GAPIT.Genotype.LD_decay.pdf", width =13, height = 6)
par(mfcol=c(1,3),mar = c(5,5,2,2))

for(i in 1:7)
{
 plot(data.all[[i]][,1]/Aver.Dis,data.all[[i]][,2], las=1,xlab="", ylim=c(-1,1),
  ylab="", main="",cex=.5,axes=FALSE,col=i,xlim=c(0,WS0/Aver.Dis))
 if(i<7)par(new=T)
}
abline(h=0,col="darkred")
axis(2,col="black",col.ticks="black",col.axis="black",ylim=c(-1,1),las=1)
axis(1,col="black",col.ticks="black",col.axis="black",xlim=c(0,WS0/Aver.Dis),las=1)
# mtext("Distance (Kb)",side=1,line=3,cex=1.2)
# mtext("R",side=2,line=3,cex=1.2)
title(main="a",xlab="Distance (Kb)",ylab="R",line=3,cex=1.2)
# legend("topright",legend=paste("+",freg.legend,sep=""),
# col=1:7,pch=1,lty=0,lwd=1,cex=0.6,
#  bty = "n", bg = par("bg"))
# par(mar = c(5,5,2,2),xpd=TRUE)

for(i in 1:7)
{
 plot(data.all[[i]][,1]/Aver.Dis,data.all[[i]][,2]^2, las=1,xlab="", ylim=c(0,1),
     ylab="", main="",axes=FALSE,cex=.5,col=i,xlim=c(0,WS0/Aver.Dis))
     # ylab="R sqaure", main="",axes=TRUE,cex=.5,col="gray60",xlim=c(0,WS0/Aver.Dis))
 if(i<7)par(new=T)
 # plot(as.numeric(fig.d[,1]),as.numeric(fig.d[,2]), las=1,xlab="Distance (Kb)", ylab="R sqaure", main="d",cex=.5,col="gray60",xlim=c(0,WS0/Aver.Dis))
}

axis(2,col="black",col.ticks="black",col.axis="black",ylim=c(0,1),las=1)
axis(1,col="black",col.ticks="black",col.axis="black",xlim=c(0,WS0/Aver.Dis),las=1)
# mtext("Distance (Kb)",side=1,line=3,cex=1.2)
# mtext("R sqaure",side=2,line=3,cex=1.2)
title(main="b",xlab="Distance (Kb)",ylab="R sqaure",line=3,cex=1.2)
dist2=dis.all
dist2[dist2>WS0]=NA
indOrder=order(dist2)
ma=cbind(as.data.frame(dist2[indOrder]),as.data.frame(R.all[indOrder]))
indRM=ma[,1]==0
maPure=ma[!indRM,]
maPure=maPure[!is.na(maPure[,1]),]
ns=nrow(maPure)
# ws=20
slide=ws
loc=matrix(NA,floor(ns/slide),2)

for (i in 1:floor(ns/slide)){
  pieceD=maPure[ ((i-1)*slide+1):((i-1)*slide+ws), 1]
  pieceR=maPure[ ((i-1)*slide+1):((i-1)*slide+ws), 2]^2
  loc[i,1]=mean(pieceD,na.rm=T)
  loc[i,2]=mean(pieceR,na.rm=T)
}


r0.hist=hist(R.all,  plot=FALSE)
r0=r0.hist$counts
r0.demo=ifelse(nchar(max(r0))<=4,1,ifelse(nchar(max(r0))<=8,1000,ifelse(nchar(max(r0))<=12,10000000,100000000000)))
r0.hist$counts=r0/r0.demo

d.V.hist=hist(dis.all, plot=FALSE)
d.V0=d.V.hist$counts
d.V0.demo=ifelse(nchar(max(d.V0))<=4,1,ifelse(nchar(max(d.V0))<=8,1000,ifelse(nchar(max(d.V0))<=12,10000000,100000000000)))

ylab0=ifelse(nchar(max(d.V0))<=4,1,ifelse(nchar(max(d.V0))<=8,2,ifelse(nchar(max(d.V0))<=12,3,4)))
ylab.store=c("Frequency","Frequency (Thousands)","Frequency (Million)","Frequency (Billion)")
d.V.hist$counts=d.V0/d.V0.demo
par(mar = c(5,2,2,5))
plot(r0.hist, xlab="R", las=1,ylab="",axes=FALSE, main="c",col="gray")

axis(4,col="black",col.ticks="black",col.axis="black")
axis(1,col="black",col.ticks="black",col.axis="black")
# mtext("",side=1,line=3,cex=1.2)
mtext(ylab.store[ylab0],side=4,line=3,cex=0.8)
# title(main="c",line=3,cex=1.2)
# lines(loc[,1]/Aver.Dis,loc[,2],lwd=4,col="gold",xlim=c(0,WS0/Aver.Dis))
legend("topleft",legend=paste("+",freg.legend,sep=""),
col=1:7,pch=1,lty=0,lwd=1,cex=1,
 bty = "n", bg = par("bg"))
grDevices::dev.off()




# print(paste("GAPIT.LD.decay ", "pdf generate.","successfully!" ,sep = ""))

#GAPIT.Genotype.View
}
#=============================================================================================
