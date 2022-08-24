`GAPIT.Genotype.View` <-function(GI=NULL,X=NULL,chr=NULL, w1_start=NULL,w1_end=NULL,mav1=NULL,
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
  
if(is.null(GI)){stop("Validation Invalid. Please select read valid Map flies  !")}
if(is.null(X)){stop("Validation Invalid. Please select read valid Genotype flies  !")}

# modified by Jiabo in 20190927. sorted number of chrom by numeric and charicter

chor_taxa=as.character(unique(GI[,2]))
chor_taxa=chor_taxa[order(as.numeric(as.character(chor_taxa)))]
chr_letter=grep("[A-Z]|[a-z]",chor_taxa)
if(!setequal(integer(0),chr_letter))
  {     
  # myGI=as.matrix(myGI)
      for(i in 1:(length(chor_taxa)))
        {
         Chr=as.character(GI[,2])
         index=Chr==chor_taxa[i]
         Chr[index]=i 
        }
      GI[,2]=as.data.frame(Chr)
  }
GI2=GI[order(as.numeric(GI[,3])),]
GI2=GI2[order(as.numeric(GI2[,2])),]
rs2=as.character(GI2[,1])
rs1=as.character(GI[,1])
index=match(rs2,rs1)
X=X[,index]
GI=GI2
# chrom=as.character(unique(GI[,2]))

if(is.null(w1_start)){w1_start=1}
##if(is.null(w1_end)){w1_end=100}
if(is.null(mav1)){mav1=10}

# if(is.null(chr)){chr=1}
chr=as.character(unique(GI[,2]))
allchr=as.character(GI[,2])
# chr=chr[order(chr)]
for(i in 1:length(chr))
{
  allchr[allchr==chr[i]]=i
}
GI[,2]=as.data.frame(allchr)
colnames(GI)[2]="Chr"
# print(table(myGI[,2]))
# map=myGI
# WS0=1e6
## make an index for marker selection with binsize
print("Filting marker for GAPIT.Genotype.View function ...")
pos.fix=as.numeric(GI[,2])*10^(nchar(max(as.numeric(GI[,3]))))+as.numeric(GI[,3])
set.seed(99163)
bins=ceiling(pos.fix/WS0)
n.bins=length(unique(bins))
uni.bins=unique(bins)

n.markers=nrow(GI)
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
  map0=GI[bins==i,]
  n.targ0=n.targ
  if(nrow(map0)<n.targ)n.targ0=nrow(map0)
  rs.markers=append(rs.markers,as.character(map0[sample(1:(nrow(map0)),n.targ0),1]))
}

rs.markers=unique(rs.markers)
rs.index=as.character(GI[,1])%in%rs.markers
print(table(rs.index))
# X=myGD
x1=X[,-ncol(X)]
x2=X[,-1]
# print("@@@@@@")
# set different colors for odd or even chromosome
#Set collor

m=ncol(X)
theCol=as.numeric(GI[,2])%%2 # here should work, based on the Chr is numeric values
colDisp=array("gray50",m-1)
colIndex=theCol==1
colDisp[colIndex]="goldenrod"
colDisp=colDisp[rs.index]
# dist=myGI[-1,3]-myGI[-nrow(myGI),3]
dist=as.numeric(GI[-1,3])-as.numeric(GI[-nrow(GI),3])
dist2=dist
index=dist<10|dist>WS0
dist[index]=NA

# myF=function(a,b) cor(a,b)
GI2=GI[rs.index,]
chr.pos=rep(NA,length(chr))
chr.pos2=rep(1,length(chr)+1)
rownames(GI2)=1:nrow(GI2)
mm=nrow(GI2)
# pos0=0
for(i in 1:length(chr))
{
  chr.pos[i]=floor(median(as.numeric(rownames(GI2[GI2[,2]==chr[i],]))))
  chr.pos2[i+1]=max(as.numeric(rownames(GI2[GI2[,2]==chr[i],])))
}
odd=seq(1,length(chr),2)
r=mapply(GAPIT.Cor.matrix,as.data.frame(x1),as.data.frame(x2))

grDevices::pdf("GAPIT.Genotype.Density_R_sqaure.pdf", width =10, height = 6)
d.V=dist/Aver.Dis
par(mfrow=c(2,3),mar = c(5,5,2,2))
plot(r[rs.index], xlab="Marker",las=1,xlim=c(1,mm), 
    ylab="R",axes=FALSE, main="a",cex=.5,col=colDisp)
axis(1,at=chr.pos2,labels=rep("",length(chr)+1))
axis(1,at=chr.pos[odd],labels=chr[odd],tick=FALSE)
axis(2,las=1)
plot(d.V[rs.index],las=1, xlab="Marker", ylab="Distance (Kb)",xlim=c(1,mm), 
    axes=FALSE,main="b",cex=.5,col=colDisp)
axis(1,at=chr.pos2,labels=rep("",length(chr)+1))
axis(1,at=chr.pos[odd],labels=chr[odd],tick=FALSE)
axis(2,las=1)

r0.hist=hist(r,  plot=FALSE)
r0=r0.hist$counts
r0.demo=ifelse(nchar(max(r0))<=4,1,ifelse(nchar(max(r0))<=8,1000,ifelse(nchar(max(r0))<=12,10000000,100000000000)))
r0.hist$counts=r0/r0.demo
ylab0=ifelse(nchar(max(r0))<=4,1,ifelse(nchar(max(r0))<=8,2,ifelse(nchar(max(r0))<=12,3,4)))
ylab.store=c("Frequency","Frequency (Thousands)","Frequency (Million)","Frequency (Billion)")
# print(nchar(max(r0)))
# print(ylab.store)
plot(r0.hist, xlab="R", las=1,ylab=ylab.store[ylab0], main="c",col="gray")
# hist(r0, xlab="R", las=1,ylab=ylab.store[ylab0], main="d")

d.V.hist=hist(d.V, plot=FALSE)
d.V0=d.V.hist$counts
d.V0.demo=ifelse(nchar(max(d.V0))<=4,1,ifelse(nchar(max(d.V0))<=8,1000,ifelse(nchar(max(d.V0))<=12,10000000,100000000000)))

ylab0=ifelse(nchar(max(d.V0))<=4,1,ifelse(nchar(max(d.V0))<=8,2,ifelse(nchar(max(d.V0))<=12,3,4)))
ylab.store=c("Frequency","Frequency (Thousands)","Frequency (Million)","Frequency (Billion)")
d.V.hist$counts=d.V0/d.V0.demo
plot(d.V.hist, las=1,xlab="Distance (Kb)",col="gray", ylab=ylab.store[ylab0], main="d",cex=.5,xlim=c(0,WS0/Aver.Dis))
plot(d.V[rs.index],r[rs.index], las=1,xlab="Distance (Kb)", ylab="R", main="e",cex=.5,col="gray60",xlim=c(0,WS0/Aver.Dis))
abline(h=0,col="darkred")
# hist(d.V0, las=1,xlab="Distance (Kb)", ylab="Frequency", main="e",cex=.5)
plot(d.V[rs.index],(r^2)[rs.index], las=1,xlab="Distance (Kb)", ylab="R sqaure", main="f",cex=.5,col="gray60",xlim=c(0,WS0/Aver.Dis))

#Moving average
dist2[dist2>WS0]=NA
indOrder=order(dist2)
ma=cbind(as.data.frame(dist2[indOrder]),as.data.frame(r[indOrder]))
indRM=ma[,1]==0
maPure=ma[!indRM,]
maPure=maPure[!is.na(maPure[,1]),]
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
lines(loc[,1]/Aver.Dis,loc[,2],col="darkred",xlim=c(0,WS0/Aver.Dis))

grDevices::dev.off()

H=1-abs(X-1)
het.ind=apply(H,1,mean)
het.snp=apply(H,2,mean)

ss=apply(X,2,sum)
maf=apply(cbind(.5*ss/(nrow(X)),1-.5*ss/(nrow(X))),1,min)
# theCol=array(0,m-1)
# for (i in 2:(m-1)){
#   if(myGI[i,2]==myGI[i-1,2])theCol[i]=theCol[i-1]
#   else theCol[i]=abs(theCol[i-1]-1)
# }


grDevices::pdf("GAPIT.Genotype.MAF_Heterozosity.pdf", width =10, height = 6)

#Display
layout.matrix <- matrix(c(1,2,3), nrow = 3, ncol = 1)
layout(mat = layout.matrix,
       heights = c(100,80,120), # Heights of the two rows
       widths = c(2, 3)) # Widths of the two columns
par(mar = c(1, 5, 1, 1))
plot(het.snp[rs.index],  las=1,ylab="Heterozygosity", xlim=c(1,mm),axes=FALSE,
    cex=.5,col=colDisp,xaxt='n')
# axis(1,at=chr.pos,labels=rep("",length(chr)))
# axis(1,at=chr.pos[odd],labels=chr[odd],tick=FALSE)
axis(2,las=1)
par(mar = c(1, 5, 0, 1))
plot(maf[rs.index], las=1,xlab="Marker", ylab="MAF",xlim=c(1,mm),axes=FALSE,
    cex=.5,col=colDisp,xaxt='n')
# axis(1,at=chr.pos,labels=rep("",length(chr)))
# axis(1,at=chr.pos[odd],labels=chr[odd],tick=FALSE)
axis(2,las=1)
par(mar = c(5, 5, 0, 1))
plot((r^2)[rs.index],  las=1,ylab="R Sqaure", xlab="Marker", xlim=c(1,mm),axes=FALSE,
    cex=.5,col=colDisp)
axis(1,at=chr.pos2,labels=rep("",length(chr)+1))
axis(1,at=chr.pos,labels=chr,tick=FALSE)
axis(2,las=1)
grDevices::dev.off()


#Display Het and MAF distribution
grDevices::pdf("GAPIT.Genotype.Distribution.pdf", width =10, height = 3.5)
layout.matrix <- matrix(c(1,2,3), nrow = 1, ncol = 3)
layout(mat = layout.matrix,
       heights = c(100,80,120), # Heights of the two rows
       widths = c(2, 2,2)) # Widths of the two columns
par(mar = c(5, 5, 2, 0))
hist(as.numeric(het.ind), las=1,xlab="Individual heterozygosity",freq=FALSE,ylab="Frequency", cex=.5,main="a")
par(mar = c(5, 4, 2, 1))
hist(het.snp[rs.index], las=1,xlab="Marker heterozygosity", freq=FALSE,ylab="Frequency",cex=.5,main="b")
par(mar = c(5, 4, 2, 1))
hist(maf[rs.index],  las=1,ylab="Frequency", xlab="MAF",freq=FALSE, cex=.5,main="c")

grDevices::dev.off()




print(paste("GAPIT.Genotype.View ", ".Three pdfs generate.","successfully!" ,sep = ""))

#GAPIT.Genotype.View
}
#=============================================================================================
