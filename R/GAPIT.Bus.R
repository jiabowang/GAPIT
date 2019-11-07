`GAPIT.Bus`<-
function(Y=NULL,CV=NULL,Z=NULL,GT=NULL,KI=NULL,GK=NULL,GD=NULL,GM=NULL,
         WS=c(1e0,1e3,1e4,1e5,1e6,1e7),alpha=c(.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),
         method=NULL,delta=NULL,vg=NULL,ve=NULL,LD=0.01,GTindex=NULL,
         cutOff=0.01,Multi_iter=FASLE,windowsize=5e10,
         p.threshold=NA,QTN.threshold=0.01,maf.threshold=0.03,
         method.GLM="FarmCPU.LM",method.sub="reward",method.sub.final="reward",method.bin="static",
         DPP=1000000,bin.size=c(5e5,5e6,5e7),bin.selection=seq(10,100,10),
		 file.output=TRUE,opt="extBIC"){
#Object: To license data by method
#Output: Coresponding numerical value
# This function is used to run multiple method, Thanks MLMM FarmCPU Blink to share program and code.
#Authors: Zhiwu Zhang
#Writen by Jiabo Wang
#Last update: Novenber 3, 2016
##############################################################################################
GR=NULL
if(method=="GLM"){
#print("---------------screening by GLM----------------------------------")

  myGAPIT <- GAPIT(
  Y=Y,			
  CV=CV,
  Z=Z,
  KI=KI,
  GD=GD,
  GM=GM,
  group.from=0,			
  group.to=0,
  QC=FALSE,
  GTindex=GTindex,
  file.output=F				
  )
  GWAS=myGAPIT$GWAS 
  GPS=myGAPIT$GPS 
  REMLs=myGAPIT$REMLs  
  delta=myGAPIT$ve/myGAPIT$va
  vg=myGAPIT$vg
  ve=myGAPIT$ve
}

#Performing first screening with MLM
if(method=="MLM"){
#print("---------------screening by MLM----------------------------------")

  myGAPIT <- GAPIT(
  Y=Y,			
  CV=CV,
  Z=Z,
  KI=KI,
  GD=GD,
  GM=GM,
  group.from=nrow(Y),			
  group.to=nrow(Y),
  QC=FALSE,
  GTindex=GTindex,
  file.output=F				
  )
  GWAS=myGAPIT$GWAS 
  GPS=myGAPIT$GPS 
  REMLs=myGAPIT$REMLs  
  delta=myGAPIT$ve/myGAPIT$va
  vg=myGAPIT$vg
  ve=myGAPIT$ve
}

#Performing first screening with Compressed MLM
if(method=="CMLM"){
#print("---------------screening by CMLM----------------------------------")
  myGAPIT <- GAPIT(
  Y=Y,			
  CV=CV,
  Z=Z,
  KI=KI,
  GD=GD,
  GM=GM,
  group.from=1,			
  group.to=nrow(Y),
  QC=FALSE,
  GTindex=GTindex,
  file.output=F				
  )
  GWAS=myGAPIT$GWAS 
  GPS=myGAPIT$GPS 
  REMLs=myGAPIT$REMLs  
  delta=myGAPIT$ve/myGAPIT$va
  vg=myGAPIT$vg
  ve=myGAPIT$ve
}

#Performing first screening with FaST-LMM
if(method=="FaST" | method=="SUPER"| method=="DC")
{
  GWAS=NULL
  GPS=NULL
  if(!is.null(vg) & !is.null(vg) & is.null(delta)) delta=ve/vg
  if(is.null(vg) & is.null(ve))
  {
    #print("!!!!!!!!!!!!!!!!")
    myFaSTREML=GAPIT.get.LL(pheno=matrix(Y[,-1],nrow(Y),1),geno=NULL,snp.pool=as.matrix(GK[,-1]),X0=as.matrix(cbind(matrix(1,nrow(CV),1),CV[,-1])))
    #print(myFaSTREML)
#print("Transfer data...")    
    REMLs=-2*myFaSTREML$LL  
    delta=myFaSTREML$delta
    vg=myFaSTREML$vg
    ve=myFaSTREML$ve
    #GPS=myFaSTREML$GPS
  }

mySUPERFaST=GAPIT.SUPER.FastMLM(ys=matrix(Y[,-1],nrow(Y),1),X0=as.matrix(cbind(matrix(1,nrow(CV),1),CV[,-1])),snp.pool=as.matrix(GK[-1]), xs=as.matrix(GD[GTindex,-1]),vg=vg,delta=delta,LD=LD,method=method)
GWAS=cbind(GM,mySUPERFaST$ps,mySUPERFaST$stats,mySUPERFaST$dfs,mySUPERFaST$effect)
}#End of if(method=="FaST" | method=="SUPER")


if(method=="FarmCPU")
{
  if(!require(bigmemory)) install.packages("bigmemory")
  if(!require(biganalytics)) install.packages("biganalytics")
library(bigmemory)  #for FARM-CPU
library(biganalytics) #for FARM-CPU
if(!exists('FarmCPU', mode='function'))source("http://www.zzlab.net/FarmCPU/FarmCPU_functions.txt")#web source code
colnames(GM)[1]="SNP"

#print(GTindex)
if(!is.null(CV))
{       farmcpuCV=CV[,2:ncol(CV)]
  }else{
        farmcpuCV=NULL
}

myFarmCPU=FarmCPU(
        Y=Y,#Phenotype
        GD=GD,#Genotype
        GM=GM,#Map information
        CV=farmcpuCV,
        cutOff=cutOff,p.threshold=p.threshold,QTN.threshold=QTN.threshold,
        maf.threshold=maf.threshold,method.GLM=method.GLM,method.sub=method.sub,
        method.sub.final=method.sub.final,method.bin=method.bin,bin.size=c(5e5,5e6,5e7),bin.selection=seq(10,100,10),
        file.output=T
        )
GWAS=myFarmCPU$GWAS
 xs=t(GD[,-1])
# gene_taxa=colnames(GD)[-1]
 ss=apply(xs,1,sum)
 ns=nrow(GD)
# storage=cbind(.5*ss/ns,1-.5*ss/ns)
# #maf=as.data.frame(cbind(gene_taxa,apply(cbind(.5*ss/ns,1-.5*ss/ns),1,min)))
# colnames(maf)=c("SNP","maf")
 nobs=ns
# #GWAS=merge(GWAS[,-5],maf, by.x = "SNP", by.y = "SNP")
# GWAS=merge(GWAS[,-5],maf, by.x = "SNP", by.y = "SNP")  #Jiabo modified at 2019.3.25
 GWAS=cbind(GWAS,nobs)

maf=apply(cbind(.5*ss/ns,1-.5*ss/ns),1,min)
GWAS$maf=maf
#print(head(GWAS))
GR=GAPIT.RandomModel(Y=Y,X=GD[,-1],GWAS=GWAS,CV=cbind(Y[,1],farmcpuCV),cutOff=cutOff)

#print("!!!!!!!!!!!!!")
#print(Multi_iter)
if(Multi_iter)
{
sig=GWAS[GWAS[,4]<(cutOff/(nrow(GWAS))),1:5]
sig=sig[!is.na(sig[,4]),]
#windowsize=500000000
sig_position=as.numeric(as.matrix(sig[,1:3])[,2])*10^10+as.numeric(as.matrix(sig[,1:3])[,3])
sig=sig[order(sig_position),]
sig_position=sig_position[order(sig_position)]
sig_diff=abs(sig_position-c(sig_position[-1],0))
sig_diff_index=sig_diff<windowsize
GWAS0=GWAS
#####################
n=nrow(sig)
print("The number of significant markers is")
print(n)

 if(n>0)
 {
  for(i in 1:n)
  {
    aim_marker=sig[i,]
    #print(aim_marker)
    aim_order=as.numeric(rownames(aim_marker))
    aim_chro=as.character(aim_marker[,2])
    aim_position=as.numeric(as.character(aim_marker[,3]))
    position=as.numeric(as.matrix(GM)[,3])
    aim_area=GM[,2]==aim_chro&position<(aim_position+windowsize)&position>(aim_position-windowsize)    
    aim_matrix=as.matrix(table(aim_area))
    #print(aim_area)
    #print(setequal(aim_area,logical(0)))
    if(setequal(aim_area,logical(0))) next
    if(aim_matrix[rownames(aim_matrix)=="TRUE",1]<10) next
        aim_area[GM[,1]==aim_marker[,1]]=FALSE      
        secondGD=GD[,c(TRUE,aim_area)]
        secondGM=GM[aim_area,]
        myGAPIT_Second <- FarmCPU(
                        Y=Y,
                        GD=secondGD,
                        GM=secondGM,
                        CV=farmcpuCV,
                        file.output=T
                        )
        Second_GWAS= myGAPIT_Second$GWAS [,1:4]
        Second_GWAS[is.na(Second_GWAS[,4]),4]=1
        orignal_GWAS=GWAS[aim_area,]
        GWAS_index=match(Second_GWAS[,1],GWAS[,1])
        #test_GWAS=GWAS
        GWAS[GWAS_index,4]=Second_GWAS[,4]
   }
 }
}

#GWAS=GWAS[order(GWAS$P.value),]
#colnames(GWAS)=c("SNP","Chromosome","Position","mp","mc","maf","nobs")
#print(head(GWAS))
GWAS[,2]=as.numeric(as.character(GWAS[,2]))
GWAS[,3]=as.numeric(as.character(GWAS[,3]))
#rint(head(GWAS))
nobs=ns
# #GWAS=merge(GWAS[,-5],maf, by.x = "SNP", by.y = "SNP")
# GWAS=merge(GWAS[,-5],maf, by.x = "SNP", by.y = "SNP")  #Jiabo modified at 2019.3.25
#GWAS=cbind(GWAS,nobs)

#print(head(GWAS))
GWAS=GWAS[,c(1:5,7,6)]

GPS=myFarmCPU$Pred
#colnames(GPS)[3]=c("Prediction")

h2=NULL
vg=NULL
ve=NULL
delta=NULL
REMLs=NULL
#print(dim(GWAS))
#print(head(GWAS))
print("FarmCPU has been done succeedly!!")
}
if(method=="BlinkC")
{
blink_GD=t(GD[,-1])
blink_GM=GM
blink_Y=Y
blink_Y[is.na(blink_Y)]="NaN"
colnames(blink_Y)=c("taxa","trait1")
blink_CV=CV
write.table(blink_GD,"myData.dat",quote=F,col.names=F,row.names=F)
write.table(blink_GM,"myData.map",quote=F,col.names=T,row.names=F)
write.table(blink_Y,"myData.txt",quote=F,col.names=T,row.names=F)
if(!is.null(CV))
{
  write.table(blink_CV,"myData.cov",quote=F,col.names=T,row.names=F)
}else{
  system("rm myData.cov")
}
system("./blink --gwas --file myData --numeric")
result=read.table("trait1_GWAS_result.txt",head=T)
result=result[,c(1,2,3,5,4)]
xs=t(GD[,-1])
#print(dim(xs))
gene_taxa=colnames(GD)[-1]
ss=apply(xs,1,sum)
ns=nrow(GD)
storage=cbind(.5*ss/ns,1-.5*ss/ns)
maf=result[,5]
#colnames(maf)=c("SNP","maf")
nobs=ns
effect=rep(NA,length(nobs))
#myFarmCPU$GWAS=merge(myFarmCPU$GWAS[,-5],maf, by.x = "SNP", by.y = "SNP")
GWAS=cbind(result[,1:4],effect)
GWAS=cbind(GWAS,maf)
GWAS=cbind(GWAS,nobs)
GWAS[,2]=as.numeric(as.character(GWAS[,2]))
GWAS[,3]=as.numeric(as.character(GWAS[,3]))
#print(dim(GWAS))
#GWAS=GWAS[order(GWAS$P.value),]
colnames(GWAS)=c("SNP","Chromosome","Position","P.value","effec","maf","nobs")

GPS=NULL
#colnames(GPS)[3]=c("Prediction")

h2=NULL
vg=NULL
ve=NULL
delta=NULL
REMLs=NULL
}
if(method=="Blink")
{
  if(!require(devtools))  install.packages("devtools")
  if(!require(BLINK)) devtools::install_github("YaoZhou89/BLINK")
  #source("http://zzlab.net/GAPIT/gapit_functions.txt")
  source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")
  blink_GD=t(GD[,-1])
  blink_GM=GM
  blink_Y=Y
  blink_CV=NULL
  if(!is.null(CV))blink_CV=CV[,-1]

  #print(head(blink_CV))
  library(BLINK)

  myBlink=Blink(Y=blink_Y,GD=blink_GD,GM=blink_GM,CV=blink_CV,maxLoop=10,time.cal=T)
  #print(head(myBlink$GWAS))
  GWAS=myBlink$GWAS[,1:4]
  #print(str(myBlink))
  ns=nrow(GD)
  nobs=ns
  effect=rep(NA,length(nobs))
  xs=t(GD[,-1])
  ss=apply(xs,1,sum)
  #storage=cbind(.5*ss/ns,1-.5*ss/ns)
  maf=apply(cbind(.5*ss/ns,1-.5*ss/ns),1,min)

  GWAS=cbind(GWAS,maf)#, by.x = "SNP", by.y = "SNP")  #Jiabo modified at 2019.3.25
  GWAS=cbind(GWAS,effect)
  GWAS=cbind(GWAS,nobs)
  GR=GAPIT.RandomModel(Y=blink_Y,X=GD[,-1],GWAS=GWAS,CV=CV,cutOff=cutOff)

#   gene_taxa=as.character(blink_GM[,1])
#   ss=apply(blink_GD,1,sum)
#   ns=nrow(GD)
#   nobs=ns
#   storage=cbind(.5*ss/ns,1-.5*ss/ns)
#   maf=as.data.frame(cbind(gene_taxa,apply(cbind(.5*ss/ns,1-.5*ss/ns),1,min)))
#   colnames(maf)=c("SNP","maf")
#   effect=rep(NA,length(nobs))
#   GWAS=cbind(GWAS,effect)
#   GWAS=cbind(GWAS,maf)
#   GWAS=cbind(GWAS,nobs)

# GWAS[,2]=as.numeric(as.character(GWAS[,2]))
# GWAS[,3]=as.numeric(as.character(GWAS[,3]))

# GPS=myBlink$Pred
#colnames(GPS)[3]=c("Prediction")

if(Multi_iter)
{
sig=GWAS[GWAS[,4]<(cutOff/(nrow(GWAS))),1:5]
sig=sig[!is.na(sig[,4]),]
#windowsize=500000000
sig_position=as.numeric(as.matrix(sig[,1:3])[,2])*10^10+as.numeric(as.matrix(sig[,1:3])[,3])
sig=sig[order(sig_position),]
sig_position=sig_position[order(sig_position)]
sig_diff=abs(sig_position-c(sig_position[-1],0))
sig_diff_index=sig_diff<windowsize
GWAS0=GWAS
#####################
n=nrow(sig)
print("The number of significant markers is")
print(n)

if(n>0)
{
  for(i in 1:n)
  {
    aim_marker=sig[i,]
    #print(aim_marker)
    aim_order=as.numeric(rownames(aim_marker))
    aim_chro=as.character(aim_marker[,2])
    aim_position=as.numeric(as.character(aim_marker[,3]))
    position=as.numeric(as.matrix(GM)[,3])
    aim_area=GM[,2]==aim_chro&position<(aim_position+windowsize)&position>(aim_position-windowsize)    
    aim_matrix=as.matrix(table(aim_area))
    if(setequal(aim_area,logical(0))) next

    if(aim_matrix[rownames(aim_matrix)=="TRUE",1]<10) next
        aim_area[GM[,1]==aim_marker[,1]]=FALSE      
        secondGD=GD[,c(TRUE,aim_area)]
        secondGM=GM[aim_area,]
        myGAPIT_Second =Blink(Y=Y,GD=secondGD,GM=secondGM,CV=blink_CV,maxLoop=10,time.cal=T)
        #print(head(myBlink$GWAS))
        #GWAS=myBlink$GWAS[,1:4]
        Second_GWAS= myGAPIT_Second$GWAS [,1:4]
        Second_GWAS[is.na(Second_GWAS[,4]),4]=1
        orignal_GWAS=GWAS[aim_area,]
        GWAS_index=match(Second_GWAS[,1],GWAS[,1])
        #test_GWAS=GWAS
        GWAS[GWAS_index,4]=Second_GWAS[,4]
  }
}



}

#GWAS=GWAS[order(GWAS$P.value),]
#colnames(GWAS)=c("SNP","Chromosome","Position","mp","mc","maf","nobs")
#print(head(GWAS))
GWAS[,2]=as.numeric(as.character(GWAS[,2]))
GWAS[,3]=as.numeric(as.character(GWAS[,3]))
#rint(head(GWAS))

GPS=myBlink$Pred
#colnames(GPS)[3]=c("Prediction")
#print(head(GWAS))
GWAS=GWAS[,c(1:5,7,6)]


h2=NULL
vg=NULL
ve=NULL
delta=NULL
REMLs=NULL
  #print(dim(GWAS))
  #print(head(GWAS))
  print("Bink R is done !!!!!")
}
if(method=="MLMM")
{
print("GWAS by MLMM method !!")
Y=Y[!is.na(Y[,2]),]
taxa_Y=as.character(Y[,1])
taxa_GD=as.character(GD[,1])
taxa_CV=as.character(CV[,1])
GD=GD[taxa_GD%in%taxa_Y,]
CV=CV[taxa_CV%in%taxa_Y,]

#print(dim(Y))
#print(dim(GD))
if(is.null(KI))
{
KI= GAPIT.kinship.VanRaden(snps=as.matrix(GD[,-1]))
colnames(KI)=as.character(GD[,1])
}else{
print("The Kinship is provided by user !!")
KI=KI[,-1] 
#colnames(KI)=as.character(GD[,1])
taxa_KI=as.character(colnames(KI))
KI=KI[taxa_KI%in%as.character(GD[,1]),taxa_KI%in%as.character(GD[,1])]
}
if(ncol(KI)!=nrow(GD)) print("Please make sure dim of K equal number of GD !!")

# print(dim(KI))
# print(dim(GD))
# print(KI[1:5,1:5])

if(is.null(CV))
{
mymlmm=mlmm(
Y=Y[,2],#Phenotype
X=as.matrix(GD[,-1]),#Genotype
K=as.matrix(KI),
#cofs=CV[,2:ncol(CV)],
nbchunks = 2, maxsteps = 10, thresh = 1.2 * 10^-5)

}else{
mymlmm=mlmm_cof(
Y=Y[,2],#Phenotype
X=as.matrix(GD[,-1]),#Genotype
K=as.matrix(KI),
cofs=as.matrix(CV[,2:ncol(CV)]),
nbchunks = 2, maxsteps = 10, thresh = 1.2 * 10^-5)
}
if(opt=='extBIC'){
GWAS_result=mymlmm$opt_extBIC$out
}
if(opt=='mbonf'){
GWAS_result=mymlmm$opt_mbonf$out
}
if(opt=='thresh'){
GWAS_result=mymlmm$opt_thresh$out
}
colnames(GWAS_result)=c("SNP","P.value")
xs=t(GD[,-1])
#print(dim(xs))
gene_taxa=colnames(GD)[-1]
colnames(GM)=c("SNP","Chromosome","position")
ss=apply(xs,1,sum)
ns=nrow(GD)
storage=cbind(.5*ss/ns,1-.5*ss/ns)
maf=as.data.frame(cbind(gene_taxa,apply(cbind(.5*ss/ns,1-.5*ss/ns),1,min)))
colnames(maf)=c("SNP","maf")
nobs=ns
GWAS_GM=merge(GM,GWAS_result, by.x = "SNP", by.y = "SNP")
mc=matrix(NA,nrow(GWAS_GM),1)
GWAS_GM=cbind(GWAS_GM,mc)
#print(head(GWAS_GM))
#print(head(maf))
#maf=NULL
GWAS_GM_maf=merge(GWAS_GM,maf, by.x = "SNP", by.y = "SNP")

GWAS=cbind(GWAS_GM_maf,nobs)
#print(head(GWAS))
GWAS=GWAS[order(GWAS$P.value),]
GWAS[,2]=as.numeric(as.character(GWAS[,2]))
GWAS[,3]=as.numeric(as.character(GWAS[,3]))
GPS=NULL
#h2=mymlmm$step_table$h2[length(mymlmm$step_table$h2)]
h2=NULL
vg=NULL
ve=NULL
delta=NULL
REMLs=NULL
colnames(GWAS)=c("SNP","Chromosome","Position","P.value","effect","maf","nobs")

}
#print("GAPIT.Bus succeed!")  
return (list(GWAS=GWAS, GPS=GPS,REMLs=REMLs,vg=vg,ve=ve,delta=delta,GVs=GR$GVs))
} #end of GAPIT.Bus
#=============================================================================================









