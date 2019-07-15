`GAPIT.Power` <-
function(WS=c(1e0,1e3,1e4,1e5,1e6,1e7), GM=NULL,seqQTN=NULL,GWAS=NULL,maxOut=100,
alpha=c(.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),MaxBP=1e10){
#Object: To evaluate power and FDR for the top (maxOut) positive interval defined by WS
#Input: WS- window size 
#Input: GM - m by 3  matrix for SNP name, chromosome and BP
#Input: seqQTN - s by 1 vecter for index of QTN on GM (+1 for GDP column wise)
#Input: GWAS- SNP,CHR,BP,P,MAF
#maxOut: maximum number of rows to report
#Requirement: None
#Output: Table and Plots
#Authors: Zhiwu Zhang
# Date  start: April 2, 2013
# Last update: April 2, 2013
##############################################################################################
#print("GAPIT.Power Started")
if(is.null(seqQTN) | is.null(GM) | is.null(GWAS)) return(list(FDR=NULL,Power=NULL,Power.Alpha=NULL,alpha=NULL))

#-----------------FDR and Power analysis-------------------------
#Information needed: myGAPIT$GWAS,myGM and QTN(r)
nWin=matrix(NA,length(WS),1)

format_GWAS=cbind(GWAS[,1:4],NA,NA,NA) 

names(format_GWAS)<-c("SNP","Chromosome","Position","P.value","maf","nobs","FDR_Adjusted_P-values")
myGM=GM

#loop window size here

theWS=1
for (theWS in 1:length(WS)){

ws=WS[theWS]
#Label QTN intervals
#Restore original order
#QTNList=r-1
QTNList=seqQTN
myGM2=cbind(myGM,rep(0,nrow(myGM)),1:nrow(myGM),NA) #Initial QTN status as 0


#Extract QTN positions
myGM2[,6]=floor((as.numeric(as.character(myGM2[,2]))*MaxBP+as.numeric(as.character(myGM2[,3])))/ws) #Label QTN as 1

QTNInterval=myGM2[QTNList,6]
thePosition=myGM2[,6] %in% QTNInterval

myGM2[thePosition,4]=1 #Label QTN as 1
names(myGM2) <- c("SNP","Chromosome","Position", "QTN","Seq") 

#Merge to P vlaues
#GWAS<- merge(myGAPIT$GWAS[,1:7],myGM2[,c(1,4,5)],by="SNP")
    GWAS<- merge(format_GWAS[,1:7],myGM2[,c(1,4,5)],by="SNP")#xiaoalei changed

#checking
#zw=GWAS[order(GWAS[,4],decreasing = FALSE),]
#zw=GWAS[order(GWAS[,8],decreasing = TRUE),]
#head(zw)

#Creat windows
myQTN=GAPIT.Specify(GI=GWAS[,1:3],GP=GWAS,bin.size=ws,MaxBP=MaxBP)
QTN=GWAS[myQTN$index,]

#Calculate alpha
qtnLoc=which(QTN[,8]==1) #get the position of QTN
P.QTN=QTN[qtnLoc,4] #p value of QTN
P.marker=QTN[-qtnLoc,4] #p value of non qtn (marker)
cutOff=matrix(quantile(P.marker, alpha,na.rm=TRUE),ncol=1)#xiaoalei changed
myPower.Alpha=apply(cutOff,1,function(x){
  Power=length(which(P.QTN<x))/length(P.QTN)
})

      
#Sort on P
#QTN=QTN[order(as.numeric(as.character(QTN[,3])),decreasing = FALSE),]
#QTN=QTN[order(as.numeric(as.character(QTN[,2])),decreasing = FALSE),]
QTN=QTN[order(as.numeric(as.character(QTN[,4])),decreasing = FALSE),]
names(QTN) <- c("SNP","Chromosome","Position", "P","FDR","Power","Order","QTN","Seq") 

#calculate power
QTN[,7]=1:nrow(QTN)
QTN[,5]=cumsum(1-QTN[,8])/QTN[,7]   #FDR
QTN[,6]=cumsum(QTN[,8]) /sum(QTN[,8] ) #Power

#Save results 
if (theWS==1){
nWin=matrix(NA,length(WS),1)
FDR=array(NA,dim=c(nrow(QTN),length(WS)))
Power=array(NA,dim=c(nrow(QTN),length(WS)))
Power.Alpha=array(NA,dim=c(length(alpha),length(WS)))
}

nWin[theWS]=nrow(QTN)
FDR[1:nWin[theWS],theWS]=QTN[,5]
Power[1:nWin[theWS],theWS]=QTN[,6]
Power.Alpha[,theWS]=myPower.Alpha

}#end of window size loop 
nOut=min(maxOut,max(nWin))
index=1:nOut
return(list(FDR=FDR[index,],Power=Power[index,],Power.Alpha=Power.Alpha,alpha=alpha))
}#end of GAPIT.Power
#=============================================================================================


