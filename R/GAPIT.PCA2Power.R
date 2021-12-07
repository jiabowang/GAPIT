`GAPIT.PCA2Power` <-function(myGD=NULL,myGM=NULL,method="MLM",myPCA=NULL,rep=NULL,h2=NULL,NQTN=NULL,seed=123){
# Object: compare to Power against FDR for GLM,MLM,CMLM,ECMLM,SUPER
# Output: find the optimum number of PCA in model
# Authors: Jiabo Wang
# Last update: Feb 1, 2020 
############################################################################################## 
if(is.null(myGD)||is.null(myGM)){stop("Read data Invalid. Please select read valid flies !")}

if(is.null(rep))
	rep=50
if(is.null(h2))
	h2=0.7
if(is.null(NQTN))
	NQTN=20

X<-myGD[,-1]
taxa<-as.character(myGD[,1])



myGAPIT <- GAPIT(
       #Y=myY[,c(1,2)],
       GD=myGD,
       GM=myGM,
       #model=method[j],
       #memo="simu",
       PCA.total=5,
       file.output=F
       )
myPCA=myGAPIT$PC
##simulation phyenotype
##-------------------------##
n=nrow(X)
m=ncol(X)
npc=ncol(myPCA)-1
legend_text=paste("NUM of PCA~",1:npc,sep="")
nm=length(method)

if(!is.null(seed))set.seed(seed)
  
  power_npca=NULL
  fdr_npca=NULL
  Para=list(h2=h2,NQTN=NQTN)

j=1
for(k in 1:npc)
{
	wholepower=NULL
    wholefdr=NULL
    for(i in 1:rep)
    {
       mysimulation<-GAPIT(Para=Para,GD=myGD,GM=myGM)
       posi=mysimulation$QTN.position
       myY=mysimulation$Y
  
       print(paste("*****************","GWAS by GAPIT...",method[j]," model ",i,sep=""))

       myGAPIT <- GAPIT(
       Y=myY[,c(1,2)],
       GD=myGD,
       GM=myGM,
       model=method[j],
       memo="simu",
       Multi_iter=F,
       file.output=F
       )
       mypower<-GAPIT.Power(WS=c(1), maxOut=m,seqQTN=posi,GM=myGM,GWAS=myGAPIT$GWAS)
       wholepower=cbind(wholepower,mypower$Power)
       wholefdr=cbind(wholefdr,mypower$FDR)
       gc()
    }


    power_rep=apply(wholepower,1,mean)
    fdr_rep=apply(wholefdr,1,mean)
    power_npca=cbind(power_npca,power_rep)
    fdr_npca=cbind(fdr_npca,fdr_rep)

} # end of npca


utils::write.csv(cbind(power_npca,fdr_npca),paste(h2,"_",NQTN,"_",method[j],".Power.by.FDR_rep_",rep,".csv",sep=""))
# write.csv(power_rep,paste(h2,"_",NQTN,"_",method[j],".Power.by.FDR_rep_",rep,".csv",sep=""))

    grDevices::pdf(paste("GAPIT.Power_",h2,"_",NQTN,"_" ,"compare in ",method[j], ".pdf", sep = ""), width = 4.5, height = 4.5,pointsize=9)
    graphics::par(mar = c(5,6,5,3))
	#win.graph(width=6, height=4, pointsize=9)
	#palette(c("blue","red","green4","brown4","orange",rainbow(5)))
	ncol=grDevices::rainbow(npc)
	grDevices::palette(c("green4","red","blue","brown4","orange",grDevices::rainbow(npc)))
	plot(power_npca[,1]~fdr_npca[,1],bg="lightgray",xlab="FDR",ylab="Power",ylim=c(0,1),xlim=c(0,1),main="Power against FDR",type="o",pch=20,col=ncol[1],cex=1,cex.lab=1.3, cex.axis=1, lwd=1,las=1)
    for(i in 2:npc){
    graphics::lines(power_npca[,i]~fdr_npca[,i], lwd=1,type="o",pch=20,col=ncol[i])
	}
	# lines(rep.power.CMLM[,6]~rep.FDR.CMLM[,6], lwd=2,type="o",pch=20,col=3)
	# lines(rep.power.MLM[,6]~rep.FDR.MLM[,6], lwd=2,type="o",pch=20,col=4)
	# lines(rep.power.GLM[,6]~rep.FDR.GLM[,6], lwd=2,type="o",pch=20,col=5)
	graphics::legend("bottomright",legend_text, pch = 20, lty =1,col=ncol,lwd=1,cex=1.0,bty="n")
	#

grDevices::dev.off()

rm(myGAPIT)
} #end of whole function




