`GAPIT.Power.compare.plink` <-function(myG=NULL,myGD=NULL,myGM=NULL,myKI=NULL,myY=NULL,myCV=NULL,rel=NULL,h2=NULL,NQTN=NULL){
# Object: compare to Power against FDR for GLM,MLM,CMLM,ECMLM,SUPER,PLINK
# rel:repetition times
# Authors: You Tang 
# Last update: January 23, 2015
############################################################################################## 

if(is.null(myG)||is.null(myGD)||is.null(myGM)||is.null(myKI)){stop("Read data Invalid. Please select read valid flies !")}
if(is.null(rel))
	rel=100
if(is.null(h2))
	h2=0.85
if(is.null(NQTN))
	NQTN=5
X<-myGD[,-1]
taxa<-myGD[,1]
taxa<-as.character(taxa)

##simulation phyenotype
##-------------------------##
n=nrow(X)
m=ncol(X)

####handle plink tped output work direct####
G<-myG[-1,]
GD<-t(X)
v3<-matrix(0,nrow(G),1)
kk<-cbind(data.frame(G[,3]),data.frame(G[,1]),data.frame(v3),data.frame(G[,4]),GD)
b1<-nrow(kk)
b2<-ncol(kk)
for(i in 1:b1){
	for(j in 5:b2){
##imput number 1
	if(kk[i,j]==0)
		kk[i,j]=1
	}
}
kk4<-cbind(kk[,5],kk[,5])
for(j in 6:b2){
	kk4<-cbind(data.frame(kk4),data.frame(kk[,j]),data.frame(kk[,j]))
}
kk6<-cbind(kk[,1:4],kk4)
##output plink deal with tped
utils::write.table(data.frame(kk6),"mdp_numeric.tped",row.names = FALSE,col.names = FALSE,sep="\t",quote=FALSE)

################----------end tped for pinlk-------------##########

rep.power.GLM<-data.frame(matrix(0,100,6))
rep.FDR.GLM<-data.frame(matrix(0,100,6))
rep.Power.Alpha.GLM<-data.frame(matrix(0,12,6))

rep.power.MLM<-data.frame(matrix(0,100,6))
rep.FDR.MLM<-data.frame(matrix(0,100,6))
rep.Power.Alpha.MLM<-data.frame(matrix(0,12,6))


rep.power.SUPER<-data.frame(matrix(0,100,6))
rep.FDR.SUPER<-data.frame(matrix(0,100,6))
rep.Power.Alpha.SUPER<-data.frame(matrix(0,12,6))

rep.power.CMLM<-data.frame(matrix(0,100,6))
rep.FDR.CMLM<-data.frame(matrix(0,100,6))
rep.Power.Alpha.CMLM<-data.frame(matrix(0,12,6))

rep.power.ECMLM<-data.frame(matrix(0,100,6))
rep.FDR.ECMLM<-data.frame(matrix(0,100,6))
rep.Power.Alpha.ECMLM<-data.frame(matrix(0,12,6))

#####------handle tfam for plink----###
rep.power.plink<-data.frame(matrix(0,100,6))
rep.FDR.plink<-data.frame(matrix(0,100,6))
rep.Power.Alpha.plink<-data.frame(matrix(0,12,6))

k1<-matrix(1,n,1)
k2<-matrix(-9,n,2)
##------step 1 end-----------##
WS=c(1e0,1e3,1e4,1e5,1e6,1e7)
alpha=c(.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1)
maxOut=100

##PCA
##---------------------##

PCA <- stats::prcomp(X)
PCVar<-PCA$sdev^2
myPC<-PCA$x[,1:3]
m1<-as.data.frame(myPC)

myCV<-cbind(taxa,m1)
myCV<-as.data.frame(myCV)

##-----end step 2  for tfam---###
kcv1<-matrix(1,nrow(myCV),1)
kcv<-cbind(data.frame(kcv1),myCV)
utils::write.table(kcv,"pca.txt",row.names = FALSE,col.names = FALSE,sep="\t",quote=FALSE)

for(i in 1:rel)
{
addm<-matrix(stats::rnorm(NQTN,0,1),NQTN,1)
QTN.position<-sample(1:m,NQTN,replace=FALSE)

SNPQ<-as.matrix(X[,QTN.position])
ge<-SNPQ%*%addm

vg <- stats::var(ge)
ve<-vg*(1-h2)/h2
SDE<-sqrt(ve)
res <- stats::rnorm(n,0,SDE)

y=as.data.frame(ge+res)
myY<-cbind(taxa,y)
myY<-as.data.frame(myY)

##-----output tfam for plink----##
k3<-cbind(data.frame(k1),data.frame(taxa),data.frame(k2),data.frame(k1),data.frame(myY[,2]))
utils::write.table(k3,paste("mdp_numeric",i,".tfam",sep=""),row.names = FALSE,col.names = FALSE,sep="\t",quote=FALSE)
 
##-----end step 2  for tfam---###

max.groups=nrow(y)
print(paste("*****************","GWAS by GAPIT...GLM model",i," totle:",rel,sep=""))
#--------------------------
myGAPIT_GLM=GAPIT(
Y=myY,
GD=myGD,
GM=myGM,
PCA.total=3,
file.output=FALSE,
group.from=0,
group.to=0,
group.by=0,
memo="GLM",
QTN.position=QTN.position,
threshold.output=0.001,
iteration.output=TRUE,
) 

print(paste("*****************","GWAS by GAPIT...MLM model",i," totle:",rel,sep=""))
#--------------------------------#
myGAPIT_MLM=GAPIT(
Y=myY,
GD=myGD,
GM=myGM,
KI=myKI,
CV=myCV,

file.output=FALSE,
group.from=max.groups,
group.to=max.groups,
group.by=10,
memo="MLM",
QTN.position=QTN.position,
threshold.output=0.001,
iteration.output=TRUE,
) 

print(paste("*****************","GWAS by GAPIT...SUPER model",i," totle:",rel,sep=""))
##--------------------------------#
myGAPIT_SUPER <- GAPIT(
Y=myY,
GD=myGD,
GM=myGM,
KI=myKI,
CV=myCV,
#PCA.total=3,
sangwich.top="MLM", #options are GLM,MLM,CMLM, FaST and SUPER
sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER
LD=0.1,
QTN.position=QTN.position,
threshold.output=0.001,
iteration.output=TRUE,
file.output=FALSE,
)

print(paste("$$$$$$$$$$$$$$$","GWAS by GAPIT...CMLM model",i," totle:",rel,sep=""))
#--------------------------------#
myGAPIT_CMLM=GAPIT(
Y=myY,
GD=myGD,
GM=myGM,
KI=myKI,
CV=myCV,

file.output=FALSE,
group.from=0,
group.to=max.groups,
group.by=10,
memo="CMLM",
QTN.position=QTN.position,
threshold.output=0.001,
iteration.output=TRUE,
) 

print(paste("-------------------","GWAS by GAPIT...ECMLM model",i," totle:",rel,sep=""))
#--------------------------------#
myGAPIT_ECMLM=GAPIT(
Y=myY,
G=myG,
#GD=myGD,
#GM=myGM,
#KI=myKI,
#CV=myCV,

PCA.total=3,
kinship.cluster=c("average", "complete", "ward"),
kinship.group=c("Mean", "Max"),
file.output=FALSE,
group.from=0,
group.to=max.groups,
group.by=10,
memo="ECMLM",
QTN.position=QTN.position,
threshold.output=0.001,
iteration.output=TRUE,
) 

##ecmlm power
power_ecmlm<-GAPIT.Power(WS=WS, alpha=alpha, maxOut=maxOut,seqQTN=QTN.position,GM=myGM,GWAS=myGAPIT_ECMLM$GWAS)

##-------------GAWS for plink-----##
##---output gwas.mdp_numericx.qassoc by plink.exe,so must be copy it to work path!----## 
system(paste('"plink.exe"', paste('--tped mdp_numeric.tped --tfam mdp_numeric',i,'.tfam --assoc --out gwas.mdp_numeric',i,sep='')), wait = TRUE)
##-------------GAWS for plink-----##
##---output gwas.mdp_numericx.qassoc by plink.exe,so must be copy it to work path!----## 
system(paste('"plink.exe"', paste('--tped mdp_numeric.tped --tfam mdp_numeric',i,'.tfam --covar pca.txt --linear --hide-covar  --out gwas.mdp_numeric',i,sep='')), wait = TRUE)

plinkGWAS <- utils::read.table(paste("gwas.mdp_numeric",i,".assoc.linear",sep=""),header=T)

Format_GWAS=cbind(myGM,plinkGWAS[,9],rep(NA,nrow(myGM)),rep(NA,nrow(myGM)),rep(NA,nrow(myGM))) 
names(Format_GWAS)<-c("SNP","Chromosome","Position","P.value","maf","nobs","FDR_Adjusted_P-values")
power_plink<-GAPIT.Power(WS=WS, alpha=alpha, maxOut=maxOut,seqQTN=QTN.position,GM=myGM,GWAS=Format_GWAS)
##---end powe_plink-----###
##----end step 3 for plink----###

#power #FDR #Power.Alpha
rep.power.GLM<-rep.power.GLM+myGAPIT_GLM$Power
rep.FDR.GLM<-rep.FDR.GLM+myGAPIT_GLM$FDR
rep.Power.Alpha.GLM<-rep.Power.Alpha.GLM+myGAPIT_GLM$Power.Alpha

rep.power.MLM<-rep.power.MLM+myGAPIT_MLM$Power
rep.FDR.MLM<-rep.FDR.MLM+myGAPIT_MLM$FDR
rep.Power.Alpha.MLM<-rep.Power.Alpha.MLM+myGAPIT_MLM$Power.Alpha

rep.power.SUPER<-rep.power.SUPER+myGAPIT_SUPER$Power
rep.FDR.SUPER<-rep.FDR.SUPER+myGAPIT_SUPER$FDR
rep.Power.Alpha.SUPER<-rep.Power.Alpha.SUPER+myGAPIT_SUPER$Power.Alpha


rep.power.CMLM<-rep.power.CMLM+myGAPIT_CMLM$Power
rep.FDR.CMLM<-rep.FDR.CMLM+myGAPIT_CMLM$FDR
rep.Power.Alpha.CMLM<-rep.Power.Alpha.CMLM+myGAPIT_CMLM$Power.Alpha

rep.power.ECMLM<-rep.power.ECMLM+power_ecmlm$Power
rep.FDR.ECMLM<-rep.FDR.ECMLM+power_ecmlm$FDR
rep.Power.Alpha.ECMLM<-rep.Power.Alpha.ECMLM+power_ecmlm$Power.Alpha

##----power-fdr save for mean of plink---##
rep.power.plink<-rep.power.plink+power_plink$Power
rep.FDR.plink<-rep.FDR.plink+power_plink$FDR
rep.Power.Alpha.plink<-rep.Power.Alpha.plink+power_plink$Power.Alpha
##---end sum for power of plink ----##

gc()
}
#mean
rep.power.GLM<-rep.power.GLM/rel
rep.FDR.GLM<-rep.FDR.GLM/rel
rep.Power.Alpha.GLM<-rep.Power.Alpha.GLM/rel

rep.power.MLM<-rep.power.MLM/rel
rep.FDR.MLM<-rep.FDR.MLM/rel
rep.Power.Alpha.MLM<-rep.Power.Alpha.MLM/rel


rep.power.SUPER<-rep.power.SUPER/rel
rep.FDR.SUPER<-rep.FDR.SUPER/rel
rep.Power.Alpha.SUPER<-rep.Power.Alpha.SUPER/rel

rep.power.CMLM<-rep.power.CMLM/rel
rep.FDR.CMLM<-rep.FDR.CMLM/rel
rep.Power.Alpha.CMLM<-rep.Power.Alpha.CMLM/rel

rep.power.ECMLM<-rep.power.ECMLM/rel
rep.FDR.ECMLM<-rep.FDR.ECMLM/rel
rep.Power.Alpha.ECMLM<-rep.Power.Alpha.ECMLM/rel

rep.power.plink<-rep.power.plink/rel
rep.FDR.plink<-rep.FDR.plink/rel
rep.Power.Alpha.plink<-rep.Power.Alpha.plink/rel

#ouput files power FDR for GLM,MLM,SUPER

myWS=c(1e0,1e3,1e4,1e5,1e6,1e7)
myalpha=c(.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1)

colnames(rep.FDR.GLM)=  paste("FDR(",myWS,")",sep="")
colnames(rep.power.GLM)=paste("Power(",myWS,")",sep="")
colnames(rep.Power.Alpha.GLM)=paste("Power(",myWS,")",sep="")

colnames(rep.FDR.MLM)=  paste("FDR(",myWS,")",sep="")
colnames(rep.power.MLM)=paste("Power(",myWS,")",sep="")
colnames(rep.Power.Alpha.MLM)=paste("Power(",myWS,")",sep="")

colnames(rep.FDR.SUPER)=  paste("FDR(",myWS,")",sep="")
colnames(rep.power.SUPER)=paste("Power(",myWS,")",sep="")
colnames(rep.Power.Alpha.SUPER)=paste("Power(",myWS,")",sep="")

colnames(rep.FDR.CMLM)=  paste("FDR(",myWS,")",sep="")
colnames(rep.power.CMLM)=paste("Power(",myWS,")",sep="")
colnames(rep.Power.Alpha.CMLM)=paste("Power(",myWS,")",sep="")

colnames(rep.FDR.ECMLM)=  paste("FDR(",myWS,")",sep="")
colnames(rep.power.ECMLM)=paste("Power(",myWS,")",sep="")
colnames(rep.Power.Alpha.ECMLM)=paste("Power(",myWS,")",sep="")

colnames(rep.FDR.plink)=  paste("FDR(",myWS,")",sep="")
colnames(rep.power.plink)=paste("Power(",myWS,")",sep="")
colnames(rep.Power.Alpha.plink)=paste("Power(",myWS,")",sep="")

utils::write.csv(cbind(rep.FDR.GLM,rep.power.GLM),paste(h2,"_",NQTN,".Power.by.FDR.GLM",rel,".csv",sep=""))
utils::write.csv(cbind(myalpha,rep.Power.Alpha.GLM),paste(h2,"_",NQTN,".Power.by.TypeI.GLM",".csv",sep=""))

utils::write.csv(cbind(rep.FDR.MLM,rep.power.MLM),paste(h2,"_",NQTN,".Power.by.FDR.MLM",rel,".csv",sep=""))
utils::write.csv(cbind(myalpha,rep.Power.Alpha.MLM),paste(h2,"_",NQTN,".Power.by.TypeI.MLM",".csv",sep=""))

utils::write.csv(cbind(rep.FDR.SUPER,rep.power.SUPER),paste(h2,"_",NQTN,".Power.by.FDR.SUPER",rel,".csv",sep=""))
utils::write.csv(cbind(myalpha,rep.Power.Alpha.SUPER),paste(h2,"_",NQTN,".Power.by.TypeI.SUPER",".csv",sep=""))

utils::write.csv(cbind(rep.FDR.CMLM,rep.power.CMLM),paste(h2,"_",NQTN,".Power.by.FDR.CMLM",rel,".csv",sep=""))
utils::write.csv(cbind(myalpha,rep.Power.Alpha.CMLM),paste(h2,"_",NQTN,".Power.by.TypeI.CMLM",".csv",sep=""))

utils::write.csv(cbind(rep.FDR.ECMLM,rep.power.ECMLM),paste(h2,"_",NQTN,".Power.by.FDR.ECMLM",rel,".csv",sep=""))
utils::write.csv(cbind(myalpha,rep.Power.Alpha.ECMLM),paste(h2,"_",NQTN,".Power.by.TypeI.ECMLM",".csv",sep=""))

utils::write.csv(cbind(rep.FDR.plink,rep.power.plink),paste(h2,"_",NQTN,".Power.by.FDR.plink",rel,".csv",sep=""))
utils::write.csv(cbind(myalpha,rep.Power.Alpha.plink),paste(h2,"_",NQTN,".Power.by.TypeI.plink",".csv",sep=""))

utils::write.csv(cbind(rep.FDR.GLM[,6],rep.power.GLM[,6],rep.FDR.MLM[,6],rep.power.MLM[,6],rep.FDR.CMLM[,6],rep.power.CMLM[,6],rep.FDR.ECMLM[,6],rep.power.ECMLM[,6],rep.FDR.SUPER[,6],rep.power.SUPER[,6],rep.FDR.plink[,6],rep.power.plink[,6]),paste(h2,"_",NQTN,".Power.by.FDR.GLM.MLM.SUPER.plink",rel,".csv",sep=""))
	name.of.trait=noquote(names(myY)[2])

grDevices::pdf(paste("GAPIT.Power ", name.of.trait,".compare to GLM,MLM,CMLM,ECMLM,SUPER.plink", ".pdf", sep = ""), width = 4.5, height = 4,pointsize=9)
graphics::par(mar = c(5,6,5,3))
	#win.graph(width=6, height=4, pointsize=9)
	grDevices::palette(c("green4","red","blue","brown4","orange","black",grDevices::rainbow(6)))
	plot(rep.FDR.SUPER[,6],rep.power.SUPER[,6],bg="lightgray",xlab="FDR",ylab="Power",ylim=c(0,1),xlim=c(0,1),main="Power against FDR",type="o",pch=20,col=1,cex=1.0,cex.lab=1.3, cex.axis=1, lwd=2,las=1)
        graphics::lines(rep.power.ECMLM[,6]~rep.FDR.ECMLM[,6], lwd=2,type="o",pch=20,col=2)
	graphics::lines(rep.power.CMLM[,6]~rep.FDR.CMLM[,6], lwd=2,type="o",pch=20,col=3)
	graphics::lines(rep.power.MLM[,6]~rep.FDR.MLM[,6], lwd=2,type="o",pch=20,col=4)
	graphics::lines(rep.power.GLM[,6]~rep.FDR.GLM[,6], lwd=2,type="o",pch=20,col=5)
	graphics::lines(rep.power.plink[,6]~rep.FDR.plink[,6], lwd=2,type="o",pch=20,col=6,lty =1)
	graphics::legend("bottomright",c("SUPER","ECMLM","CMLM","MLM","GLM","PLINK"), pch =c(20,20,20,20,20,20), lty =c(1,1,1,1,1,2),col=c(1:6),lwd=2,cex=1.0,bty="n")
	#

grDevices::dev.off()


###add type I error and power###

kkt<-cbind(rep.Power.Alpha.SUPER[,1],rep.Power.Alpha.ECMLM[,1],rep.Power.Alpha.CMLM[,1],rep.Power.Alpha.MLM[,1],rep.Power.Alpha.GLM[,1],rep.Power.Alpha.plink[,1])
utils::write.csv(cbind(myalpha,rep.Power.Alpha.SUPER[,1],rep.Power.Alpha.ECMLM[,1],rep.Power.Alpha.CMLM[,1],rep.Power.Alpha.MLM[,1],rep.Power.Alpha.GLM[,1],rep.Power.Alpha.plink[,1]),paste(h2,"_",NQTN,".Type I error.Power.by.FDR.GLM.MLM.SUPER",rel,".csv",sep=""))

myalpha1<-myalpha/10

grDevices::pdf(paste("GAPIT.Type I error_Power ", name.of.trait,".compare to GLM,MLM,CMLM,ECMLM,SUPER,plink", ".pdf", sep = ""), width = 6, height = 4.5,pointsize=9)
graphics::par(mar = c(5,6,5,3))
	
	grDevices::palette(c("green4","red","blue","brown4","orange","black",grDevices::rainbow(6)))
	plot(myalpha1,rep.Power.Alpha.SUPER[,1],log="x",bg="lightgray",xlab="Type I error",ylab="Power",main="Power against FDR",type="o",pch=20,col=1,cex=1.0,cex.lab=1.3, cex.axis=1, lwd=2,las=1,ylim=c(min(kkt),max(kkt)))
        graphics::lines(rep.Power.Alpha.ECMLM[,1]~myalpha1, lwd=2,type="o",pch=20,col=2)
	graphics::lines(rep.Power.Alpha.CMLM[,1]~myalpha1, lwd=2,type="o",pch=20,col=3)
	graphics::lines(rep.Power.Alpha.MLM[,1]~myalpha1, lwd=2,type="o",pch=20,col=4)
	graphics::lines(rep.Power.Alpha.GLM[,1]~myalpha1, lwd=2,type="o",pch=20,col=5)
	graphics::lines(rep.Power.Alpha.plink[,1]~myalpha1, lwd=2,type="o",pch=20,col=6,lty =1)
	graphics::legend("bottomright",c("SUPER","ECMLM","CMLM","MLM","GLM","PLINK"), pch =c(20,20,20,20,20,20), lty =c(1,1,1,1,1,2),col=c(1:6),lwd=2,cex=1.0,bty="n")

grDevices::dev.off()

print(paste("GAPIT.Power ", name.of.trait,".compare to GLM,MLM,CMLM,ECMLM,SUPER,PLINK.","successfully!" ,sep = ""))
#return(list(inf_Y_all,ref_Y_all))
}#end compare to GLM,MLM,SUPER
#=============================================================================================

