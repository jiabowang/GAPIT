`GAPIT.Bus`<-
function(Y=NULL,CV=NULL,Z=NULL,GT=NULL,KI=NULL,GK=NULL,GD=NULL,GM=NULL,
         WS=c(1e0,1e3,1e4,1e5,1e6,1e7),alpha=c(.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),
         method=NULL,delta=NULL,vg=NULL,ve=NULL,LD=0.01,GTindex=NULL,
         cutOff=0.01,Multi_iter=FALSE,num_regwas=10,Random.model=FALSE,FDRcut=FALSE,N.sig=NULL,
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
seqQTN=NULL

#print(head(CV))
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
#  if(!require(bigmemory)) install.packages("bigmemory")
#  if(!require(biganalytics)) install.packages("biganalytics")
#library(bigmemory)  #for FARM-CPU
#library(biganalytics) #for FARM-CPU
#if(!exists('FarmCPU', mode='function'))source("http://www.zzlab.net/FarmCPU/FarmCPU_functions.txt")#web source code

colnames(GM)[1]="SNP"

#print(GTindex)
if(!is.null(CV))
{       farmcpuCV=CV[,2:ncol(CV)]
  }else{
        farmcpuCV=NULL
}
#print(head(farmcpuCV))
# print(dim(GD))
# print(dim(farmcpuCV))
#print(Y)
# colnames(GD)[-1]=as.character(GM[,1])

myFarmCPU=FarmCPU(
        Y=Y,#Phenotype
        GD=GD,#Genotype
        GM=GM,#Map information
        CV=farmcpuCV,
        cutOff=cutOff,p.threshold=p.threshold,QTN.threshold=QTN.threshold,
        maf.threshold=maf.threshold,method.GLM=method.GLM,method.sub=method.sub,
        method.sub.final=method.sub.final,method.bin=method.bin,bin.size=c(5e5,5e6,5e7),bin.selection=seq(10,100,10),
        file.output=FALSE
        )
# print(head(myFarmCPU$GWAS))
seqQTN=myFarmCPU$seqQTN
seq_farm=myFarmCPU$seqQTN
# print(length(seq_farm))
taxa=names(Y)[2]
#print(taxa)
GWAS=myFarmCPU$GWAS
#print(head(GWAS))
 X=GD[,-1]
 ss=apply(X,2,sum)
 ns=nrow(GD)
 nobs=ns
 GWAS=cbind(GWAS,nobs)

maf=apply(cbind(.5*ss/ns,1-.5*ss/ns),1,min)
GWAS$maf=maf
#print(head(GWAS))
GWAS[is.na(GWAS[,4]),4]=1
GWAS2=GWAS
sig=GWAS[GWAS[,4]<(cutOff/(nrow(GWAS))),1:5]
sig_pass=TRUE
if(nrow(sig)==0)sig_pass=FALSE
# print(Multi_iter&sig_pass)
# print(Multi_iter)
# print(sig_pass)
if(Multi_iter&sig_pass)
{

sig=GWAS[GWAS[,4]<(cutOff/(nrow(GWAS))),1:5]
sig=sig[!is.na(sig[,4]),]
sig_position=as.numeric(as.matrix(sig[,2]))*10^(1+round(log10(max(as.numeric(GWAS[,3]))),0))+as.numeric(as.matrix(sig[,3]))
sig=sig[order(sig_position),]
sig_order=as.numeric(rownames(sig))
#if(setequal(sig_order,numeric(0))) break

n=nrow(sig)
if(n!=1){
  diff_order=abs(sig_order[-n]-sig_order[-1])
  diff_index=diff_order<num_regwas
  count=0
  diff_index2=count
  for(i in 1:length(diff_index))
  {
    if(!diff_index[i]) count=count+1
    diff_index2=append(diff_index2,count)
  }
}else{
  diff_order=0
  diff_index2=0
}

sig_bins=rle(diff_index2)$lengths
num_bins=length(sig_bins)

# sig_diff_index=sig_diff<windowsize
# print(sig_bins)
#GWAS0=GWAS
#####################
print("The number of significant markers is")
print(n)
if(n!=num_bins)
{
  print("The  number of significant bins is")
  print(num_bins)
}
# print(windowsize)
 if(num_bins>0)
 {
  for(i in 1:num_bins)
  { 
    n_sig=sig_bins[i]
    if(i==1)
    {  j=1:n_sig
      }else{
       j=(sum(sig_bins[1:(i-1)])+1):sum(sig_bins[1:i])
      }
    aim_marker=sig[j,]
    # print(aim_marker)
    aim_order=match(as.character(aim_marker[,1]),as.character(GM[,1]))
    aim_area=rep(FALSE,(nrow(GWAS)))
    # print(nrow(GWAS))
    

    #aim_area[c((aim_order-num_regwas):(aim_order-1),(aim_order+1):(aim_order+num_regwas))]=TRUE
    if(min(aim_order)<num_regwas)
    {
      aim_area[c(1:(max(aim_order)+num_regwas))]=TRUE
    }else{
      aim_area[c((min(aim_order)-num_regwas):(max(aim_order)+num_regwas))]=TRUE
    }
    # Next code can control with or without core marker in seconde model
    aim_area[aim_order]=FALSE  # without
    # print(table(aim_area))
    # print(length(seq_farm))
    # print(seq_farm)
    # print(seq_farm[!seq_farm%in%aim_order])
    if(!is.null(farmcpuCV))
    {
      secondCV=cbind(farmcpuCV,X[,seq_farm[!seq_farm%in%aim_order]])
    }else{
      secondCV=cbind(GD[,1],X[,seq_farm[!seq_farm%in%aim_order]])
      # secondCV=cbind(GD[,1],X[,aim_order])

    }
    # print(table(aim_area))
    # print(dim(GD))
    # aim_area=aim_area[1:(nrow(GWAS))]
    if(setequal(aim_area,logical(0))) next
        # this is used to set with sig marker in second model
        # aim_area[GM[,1]==aim_marker[,1]]=FALSE 
        # secondCV=NULL
        secondGD=GD[,c(TRUE,aim_area)]
        # print(dim(secondGD))
        secondGM=GM[aim_area,]
        print("Now that is multiple iteration for new farmcpu !!!")
        myGAPIT_Second <- FarmCPU(
                        Y=Y,
                        GD=secondGD,
                        GM=secondGM,
                        CV=secondCV,
                        file.output=F
                        )
        Second_GWAS= myGAPIT_Second$GWAS [,1:4]
        Second_GWAS[is.na(Second_GWAS[,4]),4]=1
        orignal_GWAS=GWAS[aim_area,]
        # write.csv(cbind(orignal_GWAS,Second_GWAS),paste("TEST_",i,".csv",sep=""),quote=F)

        # GWAS_index=match(Second_GWAS[,1],GWAS[,1])
        #test_GWAS=GWAS
        for(kk in 1:nrow(Second_GWAS))
        {
          GWAS_index=match(as.character(Second_GWAS[kk,1]),as.character(GWAS[,1]))
          GWAS[GWAS_index,4]=Second_GWAS[kk,4]

        }
        # GWAS[GWAS_index,4]=Second_GWAS[,4]
   }
 }
}

GWAS[,2]=as.numeric(as.character(GWAS[,2]))
GWAS[,3]=as.numeric(as.character(GWAS[,3]))
#rint(head(GWAS))
nobs=ns
# print(head(GWAS))
GWAS=GWAS[,c(1:5,7,6)]
#print(head(GWAS))
if(Random.model)GR=GAPIT.RandomModel(Y=Y,X=GD[,-1],GWAS=GWAS,CV=cbind(Y[,1],farmcpuCV),cutOff=cutOff,N.sig=N.sig,GT=GT)

GPS=myFarmCPU$Pred
#colnames(GPS)[3]=c("Prediction")

h2=NULL
vg=NULL
ve=NULL
delta=NULL
REMLs=NULL
#print(dim(GWAS))
#print(head(GWAS))
system(paste("rm -f FarmCPU.",taxa,".GWAS.Results.csv",sep=""))
system(paste("rm -f FarmCPU.",taxa,".Manhattan.Plot.Genomewise.pdf",sep=""))
system(paste("rm -f FarmCPU.",taxa,".QQ-Plot.pdf",sep=""))

print("FarmCPU has been done succeedly!!")
}
if(method=="BlinkC")
{
  print("BlinkC will be started !!")
  colnames(GD)[-1]=as.character(GM[,1])

blink_GD=t(GD[,-1])
blink_GM=GM
blink_Y=Y
blink_Y[is.na(blink_Y)]="NaN"
colnames(blink_Y)=c("taxa","trait1")
blink_CV=CV
utils::write.table(blink_GD,"myData.dat",quote=F,col.names=F,row.names=F)
utils::write.table(blink_GM,"myData.map",quote=F,col.names=T,row.names=F)
utils::write.table(blink_Y,"myData.txt",quote=F,col.names=T,row.names=F)
if(!is.null(CV))
{
  utils::write.table(blink_CV,"myData.cov",quote=F,col.names=T,row.names=F)
}else{
  system("rm myData.cov")
}
print("If there is a error without ./blink , please download the blink excute file from ")
print("https://github.com/Menggg/BLINK/blob/master/blink_mac")
print("Name it as blink. ")
print("And put it into workplace and make it executable with 'chmod 777 blink_versions' ")

system("./blink --gwas --file myData --numeric")

result = utils::read.table("trait1_GWAS_result.txt",head=T)
result=result[,c(1,2,3,5,4)]
xs=t(GD[,-1])
#print(dim(xs))
gene_taxa=as.character(GM[,1])
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
# print(dim(GWAS))
# print(head(GWAS))
#GWAS=GWAS[order(GWAS$P.value),]
colnames(GWAS)=c("SNP","Chromosome","Position","P.value","effect","maf","nobs")

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
  # if(!require(devtools))  install.packages("devtools")
  # if(!require(devtools))  system("git config --global http.proxy http://proxyuser:proxypwd@proxy.server.com:8080")
  # if(!require(devtools))  install.packages("devtools")
  #if(!require(BLINK)) devtools::install_github("YaoZhou89/BLINK", host = "api.github.com")
  # if(!require(BLINK)) devtools::install_github("jiabowang/BLINK")
#  if(!require(bigmemory)) install.packages("bigmemory")
#  if(!require(biganalytics)) install.packages("biganalytics")
#library(bigmemory)  #for FARM-CPU
#library(biganalytics) #for FARM-CPU
  #source("http://zzlab.net/GAPIT/gapit_functions.txt")
  #source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")
  colnames(GD)[-1]=as.character(GM[,1])

  blink_GD=t(GD[,-1])
  blink_GM=GM
  blink_Y=Y
  blink_CV=NULL
  if(!is.null(CV))blink_CV=CV[,-1,drop=FALSE] #Thanks for jloat's suggestion in Jul 23 2021

  #print(head(blink_CV))
  # library(BLINK)
  # source("http://zzlab.net/GAPIT/BLINK.R")
  totaltaxa=cbind(blink_Y[,1],GD[,1])
  # print(totaltaxa)
  myBlink=Blink(Y=blink_Y,GD=blink_GD,GM=blink_GM,CV=blink_CV,maxLoop=10,cutOff=cutOff,time.cal=T,FDRcut=FDRcut)
  # print(head(myBlink$GWAS))
  seqQTN=myBlink$seqQTN
  taxa=names(blink_Y)[2]
  GWAS=myBlink$GWAS[,1:4]
  #print(dim(blink_GD))
   X=GD[,-1]
 ss=apply(X,2,sum)
 ns=nrow(GD)
 nobs=ns
 GWAS=cbind(GWAS,nobs)

maf=apply(cbind(.5*ss/ns,1-.5*ss/ns),1,min)
GWAS$maf=maf
#print(head(GWAS))
GWAS[is.na(GWAS[,4]),4]=1

sig=GWAS[GWAS[,4]<(cutOff/(nrow(GWAS))),1:5]
sig_pass=TRUE
if(nrow(sig)==0)sig_pass=FALSE
# print("!!!!")
# print(Multi_iter&sig_pass)
if(Multi_iter&sig_pass)
{

sig=GWAS[GWAS[,4]<(cutOff/(nrow(GWAS))),1:5]
sig=sig[!is.na(sig[,4]),]
sig_position=as.numeric(as.matrix(sig[,1:3])[,2])*10^10+as.numeric(as.matrix(sig[,1:3])[,3])
sig=sig[order(sig_position),]
sig_order=as.numeric(rownames(sig))
#if(setequal(sig_order,numeric(0))) break

n=nrow(sig)
if(length(sig_order)!=1){
  diff_order=abs(sig_order[-length(sig_order)]-sig_order[-1])

  diff_index=diff_order<num_regwas

  count=0
  diff_index2=count
  for(i in 1:length(diff_index))
  {
    if(!diff_index[i]) count=count+1
    diff_index2=append(diff_index2,count)
  }
}else{
  diff_order=0
  diff_index2=0
}

sig_bins=rle(diff_index2)$lengths
num_bins=length(sig_bins)

# sig_diff_index=sig_diff<windowsize
#GWAS0=GWAS
#####################
print("The number of significant markers is")
print(n)
if(n!=num_bins)
{
  print("The  number of significant bins is")
  print(num_bins)
}
# print(windowsize)
 if(num_bins>0)
 {
  for(i in 1:num_bins)
  { 
    n_sig=sig_bins[i]
    if(i==1)
    {  j=1:n_sig
      }else{
       j=(sum(sig_bins[1:(i-1)])+1):sum(sig_bins[1:i])
      }
    aim_marker=sig[j,]
    #print(aim_marker)
    aim_order=as.numeric(rownames(aim_marker))
    aim_area=rep(FALSE,(nrow(GWAS)))
    # print(head(sig))
    # print(aim_order)

    #aim_area[c((aim_order-num_regwas):(aim_order-1),(aim_order+1):(aim_order+num_regwas))]=TRUE
    if(min(aim_order)<num_regwas)
    {
      aim_area[c(1:(max(aim_order)+num_regwas))]=TRUE

    }else{
      aim_area[c((min(aim_order)-num_regwas):(max(aim_order)+num_regwas))]=TRUE
    }
    # Next code can control with or without core marker in seconde model
    aim_area[aim_order]=FALSE  # without
    if(!is.null(blink_CV))
    {
      secondCV=cbind(blink_CV,X[seqQTN[!seqQTN%in%aim_order]])
    }else{
      secondCV=cbind(GD[,1],X[seqQTN[!seqQTN%in%aim_order]])

    }
    aim_area=aim_area[1:(nrow(GWAS))]
    #if(setequal(aim_area,logical(0))) next
        # this is used to set with sig marker in second model
        # aim_area[GM[,1]==aim_marker[,1]]=FALSE 
        
        secondGD=GD[,c(TRUE,aim_area)]
        secondGM=GM[aim_area,]
        print("Now that is multiple iteration for new farmcpu !!!")
        myGAPIT_Second <- Blink(
                        Y=Y,
                        GD=secondGD,
                        GM=secondGM,
                        CV=secondCV,
                        maxLoop=10,time.cal=T
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

GWAS[,2]=as.numeric(as.character(GWAS[,2]))
GWAS[,3]=as.numeric(as.character(GWAS[,3]))
#rint(head(GWAS))

effect=rep(NA,nrow(GWAS))
GWAS=cbind(GWAS,effect)
GPS=myBlink$Pred
# print(head(GWAS))
colnames(GWAS)[1:3]=c("SNP","Chromosome","Position")
GWAS=GWAS[,c(1:4,6,5,7)]
if(Random.model)GR=GAPIT.RandomModel(Y=blink_Y,X=GD[,-1],GWAS=GWAS,CV=CV,cutOff=cutOff,N.sig=N.sig,GT=GT)


h2=NULL
vg=NULL
ve=NULL
delta=NULL
REMLs=NULL

system(paste("rm -f FarmCPU.",taxa,".GWAS.Results.csv",sep=""))
system(paste("rm -f FarmCPU.",taxa,".Manhattan.Plot.Genomewise.pdf",sep=""))
system(paste("rm -f FarmCPU.",taxa,".QQ-Plot.pdf",sep=""))
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
colnames(KI)[-1]=as.character(KI[,1])
rownames(KI)=as.character(KI[,1])

taxa_KI=as.character(KI[,1])
KI=KI[,-1] 
 # print(dim(KI))
if(!is.null(CV)){
  taxa_com=intersect(taxa_KI,intersect(taxa_GD,intersect(taxa_Y,taxa_CV)))
  }else{
  taxa_com=intersect(taxa_KI,intersect(taxa_GD,taxa_Y))    
  }
# print(head(taxa_com))
KI=KI[taxa_KI%in%taxa_com,taxa_KI%in%taxa_com]
GD=GD[taxa_GD%in%taxa_com,]
Y=Y[taxa_Y%in%taxa_com,]
CV=CV[taxa_CV%in%taxa_com,]
}

if(ncol(KI)!=nrow(GD)) print("Please make sure dim of K equal number of GD !!")

# print(dim(KI))
# print(dim(GD))
# print(colnames(GD))
# print(rownames(GD))
# print(dim(CV))
 # print(KI[1:5,1:5])
rownames(GD)=1:nrow(GD)
# colnames(GD)[-1]=as.character(GM[,1])
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

#print(str(mymlmm))
if(opt=='extBIC'){
GWAS_result=mymlmm$opt_extBIC$out
effect=mymlmm$opt_extBIC$coef[-1,]
}
if(opt=='mbonf'){
GWAS_result=mymlmm$opt_mbonf$out
effect=mymlmm$opt_mbonf$coef[-1,]
}
if(opt=='thresh'){
GWAS_result=mymlmm$opt_thresh$out
effect=mymlmm$opt_thresh$coef[-1,]

}
   taxa=names(Y)[2]
   cof_marker=rownames(effect)
   effect=cbind(cof_marker,effect)
   # write.table(effect, paste("GAPIT.", taxa, ".MLMM.effect.csv", sep = ""), quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)

colnames(GWAS_result)=c("SNP","P.value")
xs=t(GD[,-1])
#print(dim(xs))
gene_taxa=as.character(GM[,1])
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
# print(dim(GWAS_GM))
#print(head(maf))
#maf=NULL
GWAS_GM_maf=merge(GWAS_GM,maf, by.x = "SNP", by.y = "SNP")
# print(dim(GWAS_GM_maf))

GWAS=cbind(GWAS_GM_maf,nobs)
# print(dim(GWAS))
GWAS=GWAS[order(GWAS$P.value),]
GWAS[,2]=as.numeric(as.character(GWAS[,2]))
GWAS[,3]=as.numeric(as.character(GWAS[,3]))
seqQTN=mymlmm$seqQTN
GPS=NULL
#h2=mymlmm$step_table$h2[length(mymlmm$step_table$h2)]
h2=NULL
vg=NULL
ve=NULL
delta=NULL
REMLs=NULL
GWAS=GWAS[,c(1:4,6,7,5)]
colnames(GWAS)=c("SNP","Chromosome","Position","P.value","maf","nobs","effect")

}
# print(head(GWAS))
#print("GAPIT.Bus succeed!")  
return (list(GWAS=GWAS, GPS=GPS,REMLs=REMLs,vg=vg,ve=ve,delta=delta,GVs=GR$GVs,seqQTN=seqQTN))
} #end of GAPIT.Bus
#=============================================================================================









