`GAPIT.Power.compare` <-function(myG=NUll,myGD=NULL,myGM=NULL,myKI=NULL,myY=NULL,myCV=NULL,rep=NULL,h2=NULL,NQTN=NULL){
# Object: compare to Power against FDR for GLM,MLM,CMLM,ECMLM,SUPER
# rep:repetition times
# Authors: You Tang & Jiabo Wang
# Last update: Feb 1, 2020 
############################################################################################## 
if(is.null(myG)||is.null(myGD)||is.null(myGM)||is.null(myKI)){stop("Read data Invalid. Please select read valid flies !")}

if(is.null(rep))
	rep=100
if(is.null(h2))
	h2=0.85
if(is.null(NQTN))
	NQTN=5

X<-myGD[,-1]
taxa<-as.character(myGD[,1])

##simulation phyenotype
##-------------------------##
n=nrow(X)
m=ncol(X)

rep.power.GLM<-data.frame(matrix(0,rep,6))
rep.FDR.GLM<-data.frame(matrix(0,rep,6))
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
##PCA
##---------------------##

PCA<-prcomp(X)
PCVar<-PCA$sdev^2
myPC<-PCA$x[,1:3]
m1<-as.data.frame(myPC)

myCV<-cbind(taxa,m1)
myCV<-as.data.frame(myCV)

##-----end step 2  for tfam---###
kcv1<-matrix(1,nrow(myCV),1)
kcv<-cbind(data.frame(kcv1),myCV)
write.table(kcv,"pca.txt",row.names = FALSE,col.names = FALSE,sep="\t",quote=FALSE)

for(i in 1:rep)
{
addm<-matrix(rnorm(NQTN,0,1),NQTN,1)
QTN.position<-sample(1:m,NQTN,replace=FALSE)

SNPQ<-as.matrix(X[,QTN.position])
ge<-SNPQ%*%addm

vg<-var(ge)
ve<-vg*(1-h2)/h2
SDE<-sqrt(ve)
res<-rnorm(n,0,SDE)

y=as.data.frame(ge+res)
myY<-cbind(taxa,y)
myY<-as.data.frame(myY)

max.groups=nrow(myY)
print(paste("*****************","GWAS by GAPIT...GLM model",i," totle:",rep,sep=""))
#--------------------------
myGAPIT_GLM=GAPIT(
Y=myY,
GD=myGD,
GM=myGM,
PCA.total=3,
file.output=FALSE,
model="GLM",
memo="GLM",
QTN.position=QTN.position,
threshold.output=0.001,
iteration.output=TRUE,
) 

print(paste("*****************","GWAS by GAPIT...MLM model",i," totle:",rep,sep=""))
#--------------------------------#
myGAPIT_MLM=GAPIT(
Y=myY,
GD=myGD,
GM=myGM,
KI=myKI,
CV=myCV,
file.output=FALSE,
model="MLM",
memo="MLM",
QTN.position=QTN.position,
threshold.output=0.001,
iteration.output=TRUE,
) 

print(paste("*****************","GWAS by GAPIT...SUPER model",i," totle:",rep,sep=""))
##--------------------------------#
myGAPIT_SUPER <- GAPIT(
Y=myY,
GD=myGD,
GM=myGM,
KI=myKI,
CV=myCV,
#PCA.total=3,
model="SUPER",
QTN.position=QTN.position,
threshold.output=0.001,
iteration.output=TRUE,
file.output=FALSE,
)

print(paste("$$$$$$$$$$$$$$$","GWAS by GAPIT...CMLM model",i," totle:",rep,sep=""))
#--------------------------------#
myGAPIT_CMLM=GAPIT(
Y=myY,
GD=myGD,
GM=myGM,
KI=myKI,
CV=myCV,
file.output=FALSE,
model="CMLM",
memo="CMLM",
QTN.position=QTN.position,
threshold.output=0.001,
iteration.output=TRUE,
) 

print(paste("-------------------","GWAS by GAPIT...ECMLM model",i," totle:",rep,sep=""))
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
model="ECMLM",
memo="ECMLM",
QTN.position=QTN.position,
threshold.output=0.001,
iteration.output=TRUE,
) 
power_ecmlm<-GAPIT.Power(WS=c(1e0,1e3,1e4,1e5,1e6,1e7), alpha=c(.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1), maxOut=100,seqQTN=QTN.position,GM=myGM,GWAS=myGAPIT_ECMLM$GWAS)

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
gc()
}
#mean
rep.power.GLM<-rep.power.GLM/rep
rep.FDR.GLM<-rep.FDR.GLM/rep
rep.Power.Alpha.GLM<-rep.Power.Alpha.GLM/rep

rep.power.MLM<-rep.power.MLM/rep
rep.FDR.MLM<-rep.FDR.MLM/rep
rep.Power.Alpha.MLM<-rep.Power.Alpha.MLM/rep

rep.power.SUPER<-rep.power.SUPER/rep
rep.FDR.SUPER<-rep.FDR.SUPER/rep
rep.Power.Alpha.SUPER<-rep.Power.Alpha.SUPER/rep

rep.power.CMLM<-rep.power.CMLM/rep
rep.FDR.CMLM<-rep.FDR.CMLM/rep
rep.Power.Alpha.CMLM<-rep.Power.Alpha.CMLM/rep

rep.power.ECMLM<-rep.power.ECMLM/rep
rep.FDR.ECMLM<-rep.FDR.ECMLM/rep
rep.Power.Alpha.ECMLM<-rep.Power.Alpha.ECMLM/rep

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

write.csv(cbind(rep.FDR.GLM,rep.power.GLM),paste(h2,"_",NQTN,".Power.by.FDR.GLM",rep,".csv",sep=""))
write.csv(cbind(myalpha,rep.Power.Alpha.GLM),paste(h2,"_",NQTN,".Power.by.TypeI.GLM",".csv",sep=""))

write.csv(cbind(rep.FDR.MLM,rep.power.MLM),paste(h2,"_",NQTN,".Power.by.FDR.MLM",rep,".csv",sep=""))
write.csv(cbind(myalpha,rep.Power.Alpha.MLM),paste(h2,"_",NQTN,".Power.by.TypeI.MLM",".csv",sep=""))

write.csv(cbind(rep.FDR.SUPER,rep.power.SUPER),paste(h2,"_",NQTN,".Power.by.FDR.SUPER",rep,".csv",sep=""))
write.csv(cbind(myalpha,rep.Power.Alpha.SUPER),paste(h2,"_",NQTN,".Power.by.TypeI.SUPER",".csv",sep=""))

write.csv(cbind(rep.FDR.CMLM,rep.power.CMLM),paste(h2,"_",NQTN,".Power.by.FDR.CMLM",rep,".csv",sep=""))
write.csv(cbind(myalpha,rep.Power.Alpha.CMLM),paste(h2,"_",NQTN,".Power.by.TypeI.CMLM",".csv",sep=""))

write.csv(cbind(rep.FDR.ECMLM,rep.power.ECMLM),paste(h2,"_",NQTN,".Power.by.FDR.ECMLM",rep,".csv",sep=""))
write.csv(cbind(myalpha,rep.Power.Alpha.ECMLM),paste(h2,"_",NQTN,".Power.by.TypeI.ECMLM",".csv",sep=""))

write.csv(cbind(rep.FDR.GLM[,6],rep.power.GLM[,6],rep.FDR.MLM[,6],rep.power.MLM[,6],rep.FDR.CMLM[,6],rep.power.CMLM[,6],rep.FDR.ECMLM[,6],rep.power.ECMLM[,6],rep.FDR.SUPER[,6],rep.power.SUPER[,6]),paste(h2,"_",NQTN,".Power.by.FDR.GLM.MLM.SUPER",rep,".csv",sep=""))
	name.of.trait=noquote(names(myY)[2])


pdf(paste("GAPIT.Power ", name.of.trait,".compare to GLM,MLM,CMLM,ECMLM,SUPER.", ".pdf", sep = ""), width = 4.5, height = 4.5,pointsize=9)
par(mar = c(5,6,5,3))
	#win.graph(width=6, height=4, pointsize=9)
	#palette(c("blue","red","green4","brown4","orange",rainbow(5)))
	palette(c("green4","red","blue","brown4","orange",rainbow(5)))
	plot(rep.FDR.SUPER[,6],rep.power.SUPER[,6],bg="lightgray",xlab="FDR",ylab="Power",ylim=c(0,1),xlim=c(0,1),main="Power against FDR",type="o",pch=20,col=1,cex=1.0,cex.lab=1.3, cex.axis=1, lwd=2,las=1)
        lines(rep.power.ECMLM[,6]~rep.FDR.ECMLM[,6], lwd=2,type="o",pch=20,col=2)
	lines(rep.power.CMLM[,6]~rep.FDR.CMLM[,6], lwd=2,type="o",pch=20,col=3)
	lines(rep.power.MLM[,6]~rep.FDR.MLM[,6], lwd=2,type="o",pch=20,col=4)
	lines(rep.power.GLM[,6]~rep.FDR.GLM[,6], lwd=2,type="o",pch=20,col=5)
	legend("bottomright",c("SUPER","ECMLM","CMLM","MLM","GLM"), pch = 20, lty =1,col=c(1:5),lwd=2,cex=1.0,bty="n")
	#

dev.off()

###add type I error and power###

kkt<-cbind(rep.Power.Alpha.SUPER[,1],rep.Power.Alpha.ECMLM[,1],rep.Power.Alpha.CMLM[,1],rep.Power.Alpha.MLM[,1],rep.Power.Alpha.GLM[,1])
write.csv(cbind(myalpha,rep.Power.Alpha.SUPER[,1],rep.Power.Alpha.ECMLM[,1],rep.Power.Alpha.CMLM[,1],rep.Power.Alpha.MLM[,1],rep.Power.Alpha.GLM[,1]),paste(h2,"_",NQTN,".Type I error.Power.by.FDR.GLM.MLM.SUPER",rep,".csv",sep=""))

myalpha1<-myalpha/10

pdf(paste("GAPIT.Type I error_Power ", name.of.trait,".compare to GLM,MLM,CMLM,ECMLM,SUPER.", ".pdf", sep = ""), width = 6, height = 4.5,pointsize=9)
par(mar = c(5,6,5,3))
	#win.graph(width=6, height=4, pointsize=9)
	#palette(c("blue","red","green4","brown4","orange",rainbow(5)))
	palette(c("green4","red","blue","brown4","orange",rainbow(5)))
	plot(myalpha1,rep.Power.Alpha.SUPER[,1],log="x",bg="lightgray",xlab="Type I error",ylab="Power",main="Power against FDR",type="o",pch=20,col=1,cex=1.0,cex.lab=1.3, cex.axis=1, lwd=2,las=1,ylim=c(min(kkt),max(kkt)))
	#plot(myalpha1,rep.Power.Alpha.SUPER[,1],bg="lightgray",xlab="Type I error",ylab="Power",ylim=c(0,1),xlim=c(0,1),main="Power against FDR",type="o",pch=20,col=1,cex=1.0,cex.lab=1.3, cex.axis=1, lwd=2,las=1)
        lines(rep.Power.Alpha.ECMLM[,1]~myalpha1, lwd=2,type="o",pch=20,col=2)
	lines(rep.Power.Alpha.CMLM[,1]~myalpha1, lwd=2,type="o",pch=20,col=3)
	lines(rep.Power.Alpha.MLM[,1]~myalpha1, lwd=2,type="o",pch=20,col=4)
	lines(rep.Power.Alpha.GLM[,1]~myalpha1, lwd=2,type="o",pch=20,col=5)
	legend("bottomright",c("SUPER","ECMLM","CMLM","MLM","GLM"), pch = 20, lty =1,col=c(1:5),lwd=2,cex=1.0,bty="n")
	#

dev.off()


print(paste("GAPIT.Power ", name.of.trait,".compare to GLM,MLM,CMLM,ECMLM,SUPER.","successfully!" ,sep = ""))
#return(list(inf_Y_all,ref_Y_all))
}#end compare to GLM,MLM,CMLM,ECMLM,SUPER
#=============================================================================================

