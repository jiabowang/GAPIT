`GAPIT.Power.compare` <-function(G=NULL,GD=NULL,GM=NULL,PCA.total=3,KI=NULL,CV=NULL,nrep=10,h2=0.85,NQTN=5,maxOut=100,model=c("GLM","FarmCPU"),seed=99163){
# Object: compare to Power against FDR for GLM,MLM,CMLM,ECMLM,SUPER
# nrep:repetition times
# Authors: You Tang & Jiabo Wang
# Last update: Mar 1, 2022
############################################################################################## 
if(is.null(GD)&is.null(GM)&is.null(G)){stop("Read data Invalid. Please select read valid flies !")}

myWS=c(1e0,1e3,1e4,1e5,1e6,1e7)
myalpha=c(.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1)
all.method=model
##simulation phyenotype
##-------------------------##
set.seed(seed)
Para=list(h2=h2,NQTN=NQTN)
rep.power.store<-list()
rep.FDR.store<-list()
rep.Power.Alpha.store<-list()
for(j in 1:length(all.method))
{
	rep.power.store[[j]]=data.frame(matrix(0,maxOut,length(myWS)))
	rep.FDR.store[[j]]<-data.frame(matrix(0,maxOut,length(myWS)))
    rep.Power.Alpha.store[[j]]<-data.frame(matrix(0,length(myalpha),length(myWS)))
}

for(i in 1:nrep)
{

mysimulation<-GAPIT(Para=Para,GD=GD,GM=GM,PCA.total=PCA.total,file.output=FALSE)
QTN.position=mysimulation$QTN.position
Y=mysimulation$Y
colnames(Y)=c("Taxa","Simu")

# max.groups=nrow(myY)
    for(j in 1:length(all.method))
    {
       myGAPIT=GAPIT(
       Y=Y,
       GD=GD,
       GM=GM,
       PCA.total=PCA.total,
       CV=CV,
       file.output=FALSE,
       model=all.method[j],
       memo=all.method[j],
       QTN.position=QTN.position
       ) 
       power_store<-GAPIT.Power(WS=myWS, alpha=myalpha, maxOut=maxOut,seqQTN=QTN.position,GM=myGM,GWAS=myGAPIT$GWAS)
       rep.power.store[[j]]<-rep.power.store[[j]]+power_store$Power
       rep.FDR.store[[j]]<-rep.FDR.store[[j]]+power_store$FDR
       rep.Power.Alpha.store[[j]]<-rep.Power.Alpha.store[[j]]+power_store$Power.Alpha
    }
gc()
}

#mean
for(j in 1:length(all.method))
{
rep.power.store[[j]]<-rep.power.store[[j]]/nrep
rep.FDR.store[[j]]<-rep.FDR.store[[j]]/nrep
rep.Power.Alpha.store[[j]]<-rep.Power.Alpha.store[[j]]/nrep

colnames(rep.FDR.store[[j]])=  paste("FDR(",myWS,")",sep="")
colnames(rep.power.store[[j]])=paste("Power(",myWS,")",sep="")
colnames(rep.Power.Alpha.store[[j]])=paste("Power(",myWS,")",sep="")

utils::write.csv(cbind(rep.FDR.store[[j]],rep.power.store[[j]]),paste(h2,"_",NQTN,".Power.by.FDR.",all.method[j],".",nrep,".csv",sep=""))
utils::write.csv(cbind(myalpha,rep.Power.Alpha.store[[j]]),paste(h2,"_",NQTN,".Power.by.TypeI.",all.method[j],".",nrep,".csv",sep=""))

}

grDevices::pdf(paste("GAPIT.Power.compare to multiple models", ".pdf", sep = ""), width = 4.5, height = 4.5,pointsize=9)
graphics::par(mar = c(5,6,5,3))
	#win.graph(width=6, height=4, pointsize=9)
	#palette(c("blue","red","green4","brown4","orange",rainbow(5)))
	plot.color=rainbow(length(all.method))
	# grDevices::palette(c(plot.color,grDevices::rainbow(length(all.method))))
kkt=NULL
for(j in 1:length(all.method))
{	
	if(j==1)plot(as.numeric(rep.FDR.store[[j]][,6]),as.numeric(rep.power.store[[j]][,6]),bg="lightgray",xlab="FDR",ylab="Power",ylim=c(0,1),xlim=c(0,1),main="Power against FDR",type="o",pch=20,col=plot.color[j],cex=1.0,cex.lab=1.3, cex.axis=1, lwd=2,las=1)
   if(j!=1) graphics::lines(rep.power.store[[j]][,6]~rep.FDR.store[[j]][,6], lwd=2,type="o",pch=20,col=plot.color[j])
	kkt<-cbind(kkt,rep.Power.Alpha.store[[j]][,1])

}
	graphics::legend("bottomright",c(all.method), pch = 20, lty =1,col=plot.color,lwd=2,cex=1.0,bty="n")
grDevices::dev.off()

###add type I error and power###
utils::write.csv(cbind(myalpha,kkt),paste(h2,"_",NQTN,".Type I error.Power.by.multiple models.",nrep,".csv",sep=""))

myalpha1<-myalpha/10
grDevices::pdf(paste("GAPIT.Type I error_Power.compare to multiple models", ".pdf", sep = ""), width = 6, height = 4.5,pointsize=9)
graphics::par(mar = c(5,6,5,3))
for(j in 1:length(all.method))
{	

	 if(j==1)plot(as.numeric(myalpha1),as.numeric(rep.Power.Alpha.store[[j]][,1]),log="x",bg="lightgray",xlab="Type I error",ylab="Power",main="Power against Type I error",type="o",pch=20,col=plot.color[j],cex=1.0,cex.lab=1.3, cex.axis=1, lwd=2,las=1,ylim=c(min(kkt),max(kkt)))
   if(j!=1) graphics::lines(rep.Power.Alpha.store[[j]][,1]~myalpha1, lwd=2,type="o",pch=20,col=plot.color[j])
	 kkt<-cbind(kkt,rep.Power.Alpha.store[[j]][,1])

}
graphics::legend("bottomright",c(all.method), pch = 20, lty =1,col=plot.color,lwd=2,cex=1.0,bty="n")
grDevices::dev.off()

}#end compare to GLM,MLM,CMLM,ECMLM,SUPER
#=============================================================================================

