`GAPIT.Genotype.View` <-function(myGI=NULL,myGD=NULL,chr=NULL, w1_start=NULL,w1_end=NULL,mav1=NULL,
                                 WS0=1e6,ws=200,Aver.Dis=1000,...){
# Object: Analysis for Genotype data:Distribution of SNP density,Accumulation,Moving Average of density,result:a pdf of the scree plot
# myG:Genotype data
# chr: chromosome value
# WS0 is the cutoff threshold for marker to display
# ws is used to calculate within windowsize
# Aver.Dis is average display windowsize
# mav1:Moving Average set value length
# Authors:  Zhiwu Zhang and Jiabo Wang
# Last update: May 15, 2022 
##############################################################################################

#if(nrow(myGI)<1000) return() #Markers are not enough for this analysis
  
if(is.null(myGI)){stop("Validation Invalid. Please select read valid Map flies  !")}
if(is.null(myGD)){stop("Validation Invalid. Please select read valid Genotype flies  !")}
if(is.null(w1_start)){w1_start=1}
##if(is.null(w1_end)){w1_end=100}
if(is.null(mav1)){mav1=10}

# if(is.null(chr)){chr=1}
chr=as.character(unique(myGI[,2]))
allchr=as.character(myGI[,2])
# chr=chr[order(chr)]
for(i in 1:length(chr))
{
  allchr[allchr==chr[i]]=i
}
myGI[,2]=as.data.frame(allchr)
colnames(myGI)[2]="Chr"
# print(table(myGI[,2]))
# map=myGI
# WS0=1e6
## make an index for marker selection with binsize
print("Filting marker for GAPIT.Genotype.View function ...")
pos.fix=as.numeric(myGI[,2])*10^(nchar(max(as.numeric(myGI[,3]))))+as.numeric(myGI[,3])
set.seed(99163)
bins=ceiling(pos.fix/WS0)
n.bins=length(unique(bins))
uni.bins=unique(bins)

n.markers=nrow(myGI)
n.select=10000
if(n.markers<n.select)n.select=n.markers
n.targ=floor(n.select/n.bins)
if(n.targ<1)
{
  n.targ=1
  uni.bins=sample(uni.bins,n.select)
}
rs.markers=NULL
for(i in uni.bins)
{
  map0=myGI[bins==i,]
  n.targ0=n.targ
  if(nrow(map0)<n.targ)n.targ0=nrow(map0)
  rs.markers=append(rs.markers,as.character(map0[sample(1:(nrow(map0)),n.targ0),1]))
}

rs.markers=unique(rs.markers)
rs.index=as.character(myGI[,1])%in%rs.markers
print(table(rs.index))
X=myGD
x1=X[,-ncol(X)]
x2=X[,-1]
# print("@@@@@@")

# dist=myGI[-1,3]-myGI[-nrow(myGI),3]
dist=as.numeric(myGI[-1,3])-as.numeric(myGI[-nrow(myGI),3])
index=dist<10|dist>WS0
dist[index]=NA

# myF=function(a,b) cor(a,b)

r=mapply(GAPIT.Cor.matrix,as.data.frame(x1),as.data.frame(x2))

grDevices::pdf("GAPIT.Marker.Density.R.sqaure.pdf", width =10, height = 6)
# plot(dist/Aver.Dis,r^2, xlab="Distance (Kb)", ylab="R sqaure", main="f",cex=.5,col="gray60")
d.V=dist/Aver.Dis
par(mfrow=c(2,3),mar = c(5,5,2,2))
plot(r[rs.index], xlab="Marker",las=1, ylab="R", main="a",cex=.5,col="gray60")
plot(d.V[rs.index],las=1, xlab="Marker", ylab="Distance (Kb)", main="b",cex=.5,col="gray60")
plot(d.V[rs.index],r[rs.index], las=1,xlab="Distance (Kb)", ylab="R", main="c",cex=.5,col="gray60")
abline(h=0,col="darkred")
r0.hist=hist(r,  plot=FALSE)
r0=r0.hist$counts
r0.demo=ifelse(nchar(max(r0))<=4,1,ifelse(nchar(max(r0))<=8,1000,ifelse(nchar(max(r0))<=12,10000000,100000000000)))
r0.hist$counts=r0/r0.demo
ylab0=ifelse(nchar(max(r0))<=4,1,ifelse(nchar(max(r0))<=8,2,ifelse(nchar(max(r0))<=12,3,4)))
ylab.store=c("Frequency","Frequency (Thousands)","Frequency (Million)","Frequency (Billion)")
# print(nchar(max(r0)))
# print(ylab.store)
plot(r0.hist, xlab="R", las=1,ylab=ylab.store[ylab0], main="d",col="gray")
# hist(r0, xlab="R", las=1,ylab=ylab.store[ylab0], main="d")
d.V.hist=hist(d.V, plot=FALSE)
d.V0=d.V.hist$counts
d.V0.demo=ifelse(nchar(max(d.V0))<=4,1,ifelse(nchar(max(d.V0))<=8,1000,ifelse(nchar(max(d.V0))<=12,10000000,100000000000)))

ylab0=ifelse(nchar(max(d.V0))<=4,1,ifelse(nchar(max(d.V0))<=8,2,ifelse(nchar(max(d.V0))<=12,3,4)))
ylab.store=c("Frequency","Frequency (Thousands)","Frequency (Million)","Frequency (Billion)")
d.V.hist$counts=d.V0/d.V0.demo
plot(d.V.hist, las=1,xlab="Distance (Kb)",col="gray", ylab=ylab.store[ylab0], main="e",cex=.5)

# hist(d.V0, las=1,xlab="Distance (Kb)", ylab="Frequency", main="e",cex=.5)
plot(d.V[rs.index],(r^2)[rs.index], las=1,xlab="Distance (Kb)", ylab="R sqaure", main="f",cex=.5,col="gray60")

#Moving average
indOrder=order(dist)
ma=cbind(as.data.frame(dist[indOrder]),as.data.frame(r[indOrder]))
indRM=ma[,1]==0
maPure=ma[!indRM,]
ns=nrow(maPure)
# ws=ws
slide=10
loc=matrix(NA,floor(ns/slide),2)

for (i in 1:floor(ns/slide)){
  pieceD=maPure[ ((i-1)*slide+1):((i-1)*slide+ws), 1]
  pieceR=maPure[ ((i-1)*slide+1):((i-1)*slide+ws), 2]^2
  loc[i,1]=mean(pieceD,na.rm=T)
  loc[i,2]=mean(pieceR,na.rm=T)
}
lines(loc[,1]/Aver.Dis,loc[,2],col="darkred")

grDevices::dev.off()

H=1-abs(X-1)
het.ind=apply(H,1,mean)
het.snp=apply(H,2,mean)

ss=apply(X,2,sum)
maf=apply(cbind(.5*ss/(nrow(X)),1-.5*ss/(nrow(X))),1,min)
#Set collor

m=ncol(X)
# theCol=array(0,m-1)
# for (i in 2:(m-1)){
#   if(myGI[i,2]==myGI[i-1,2])theCol[i]=theCol[i-1]
#   else theCol[i]=abs(theCol[i-1]-1)
# }

theCol=as.numeric(myGI[,2])%%2 # here should work, based on the Chr is numeric values
colDisp=array("gray50",m-1)
colIndex=theCol==1
colDisp[colIndex]="goldenrod"
colDisp=colDisp[rs.index]
grDevices::pdf("GAPIT.Marker.MAF.Heterozosity.pdf", width =10, height = 6)

#Display
layout.matrix <- matrix(c(1,2,3), nrow = 3, ncol = 1)
layout(mat = layout.matrix,
       heights = c(100,80,120), # Heights of the two rows
       widths = c(2, 3)) # Widths of the two columns
par(mar = c(1, 5, 1, 1))
plot(het.snp[rs.index],  las=1,ylab="Heterozygosity", cex=.5,col=colDisp,xaxt='n')
par(mar = c(1, 5, 0, 1))
plot(maf[rs.index], las=1,xlab="Marker", ylab="MAF",cex=.5,col=colDisp,xaxt='n')
par(mar = c(5, 5, 0, 1))
plot((r^2)[rs.index],  las=1,ylab="R Sqaure", xlab="Marker", cex=.5,col=colDisp)
grDevices::dev.off()


#Display Het and MAF distribution
grDevices::pdf("GAPIT.Marker.MAF.Heterozosity.distribution.pdf", width =10, height = 3.5)
layout.matrix <- matrix(c(1,2,3), nrow = 1, ncol = 3)
layout(mat = layout.matrix,
       heights = c(100,80,120), # Heights of the two rows
       widths = c(2, 2,2)) # Widths of the two columns
par(mar = c(5, 5, 2, 0))
hist(as.numeric(het.ind[rs.index]), las=1,xlab="Individual heterozygosity",freq=FALSE,ylab="Frequency", cex=.5,main="a")
par(mar = c(5, 4, 2, 1))
hist(het.snp[rs.index], las=1,xlab="Marker heterozygosity", freq=FALSE,ylab="Frequency",cex=.5,main="b")
par(mar = c(5, 4, 2, 1))
hist(maf[rs.index],  las=1,ylab="Frequency", xlab="MAF",freq=FALSE, cex=.5,main="c")

grDevices::dev.off()




print(paste("GAPIT.Genotype.View ", ".Three pdfs generate.","successfully!" ,sep = ""))

#GAPIT.Genotype.View
}
#=============================================================================================
