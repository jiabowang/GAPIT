`GAPIT.Bread` <-
function(Y=NULL,CV=NULL,Z=NULL,KI=NULL,GK=NULL,GD=NULL,GM=NULL,
              method=NULL,delta=NULL,vg=NULL,ve=NULL,LD=0.01,GTindex=NULL,
              file.output=TRUE,opt="extBIC"){
#Object: To calculate p-values of SNPs by using method of GLM, MLM, CMLM, FaST, SUPER and DC  
#Straitegy: NA
#Output: GWAS, GPS,REMLs,vg,ve,delta
#intput: 
#Y: phenotype with columns of taxa,Y1,Y2...
#CV: covariate variables with columns of taxa, v1,v2...
#GD: same as GK. This is the genotype to screen, the columns are taxa,SNP1,SNP2,...
#GK: Genotype data in numerical format, taxa goes to row and snp go ti columns. the first column is taxa
#GM: Genotype map with columns of snpID,chromosome and position
#method: Options are GLM, MLM, CMLM, FaST, SUPER ,FARM-CPU and DC 
#Authors: Zhiwu Zhang
#Last update: November 2, 2011
##############################################################################################
#print("GAPIT.SUPER in progress...")

#Performing first screening with GLM
if(method=="GLM"){
#print("---------------screening by GLM----------------------------------")
  #print(GTindex)
  myGAPIT <- GAPIT(
  Y=Y,			
  CV=CV,
  Z=Z,
  KI=KI,
  GD=GD,
  GM=GM,
  model=("GLM"),
  QC=FALSE,
  GTindex=GTindex,
  file.output=file.output				
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
  file.output=file.output				
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
  file.output=file.output				
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

    myFaSTREML=GAPIT.get.LL(pheno=matrix(Y[,-1],nrow(Y),1),geno=NULL,snp.pool=as.matrix(GK[,-1]),X0=as.matrix(cbind(matrix(1,nrow(CV),1),CV[,-1])))
    
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





#FarmCPU
if(method=="FarmCPU")
{
#  if(!require(bigmemory)) install.packages("bigmemory")
#  if(!require(biganalytics)) install.packages("biganalytics")
#library(bigmemory)  #for FARM-CPU
#library(biganalytics) #for FARM-CPU
#if(!exists('FarmCPU', mode='function'))source("http://www.zzlab.net/FarmCPU/FarmCPU_functions.txt")#web source code
colnames(GM)[1]="SNP"

myFarmCPU=FarmCPU(
Y=Y,#Phenotype
GD=GD,#Genotype
GM=GM,#Map information
CV=CV[,2:ncol(CV)],
file.output=T
)


xs=t(GD[,-1])
#print(dim(xs))
gene_taxa=colnames(GD)[-1]
ss=apply(xs,1,sum)
ns=nrow(GD)
storage=cbind(.5*ss/ns,1-.5*ss/ns)
maf=as.data.frame(cbind(gene_taxa,apply(cbind(.5*ss/ns,1-.5*ss/ns),1,min)))
colnames(maf)=c("SNP","maf")
nobs=ns
#print(dim(myFarmCPU$GWAS))
#print(length(maf))
myFarmCPU$GWAS=merge(myFarmCPU$GWAS[,-5],maf, by.x = "SNP", by.y = "SNP")
GWAS=cbind(myFarmCPU$GWAS,nobs)
GWAS=GWAS[order(GWAS$P.value),]
#colnames(GWAS)=c("SNP","Chromosome","Position","mp","mc","maf","nobs")

GPS=myFarmCPU$Pred

h2=NULL
vg=NULL
ve=NULL
delta=NULL
REMLs=NULL
#colnames(GPS)[3]=c("Prediction")
}
#MLMM
if(method=="MLMM")
{
print(" GWAS by MLMM method !!")
Y=Y[!is.na(Y[,2]),]
taxa_Y=as.character(Y[,1])
taxa_GD=as.character(GD[,1])
taxa_CV=as.character(CV[,1])
GD=GD[taxa_GD%in%taxa_Y,]
CV=CV[taxa_CV%in%taxa_Y,]

#print(dim(Y))
#print(dim(GD))
#print(dim(CV))


KI= GAPIT.kinship.VanRaden(snps=as.matrix(GD[,-1]))
colnames(KI)=as.character(GD[,1])
 
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
GWAS_GM_maf=merge(GWAS_GM,maf, by.x = "SNP", by.y = "SNP")

GWAS=cbind(GWAS_GM_maf,nobs)
#print(head(GWAS))
GWAS=GWAS[order(GWAS$P.value),]
GPS=NULL
#h2=mymlmm$step_table$h2[length(mymlmm$step_table$h2)]
h2=NULL
vg=NULL
ve=NULL
delta=NULL
REMLs=NULL
colnames(GWAS)=c("SNP","Chromosome","Position","P.value","effec","maf","nobs")

}
#print("GAPIT.Bread succeed!")  
return (list(GWAS=GWAS, GPS=GPS,REMLs=REMLs,vg=vg,ve=ve,delta=delta))
} #end of GAPIT.Bread
#=============================================================================================

