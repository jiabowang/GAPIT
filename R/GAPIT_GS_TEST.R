
rm(list=ls())

source("/Users/Jiabo/Documents/gapit/gapit_functions.txt")

setwd("/Users/Jiabo/Dropbox/GAPIT/GAPIT_Tutorial_Data")

myY=read.table(file="mdp_traits.txt", head = TRUE)
myG=read.table(file="mdp_genotype_test.hmp.txt", head = FALSE)
myGD=read.table("mdp_numeric.txt",head=T)
myGM=read.table("mdp_SNP_information.txt",head=T)
# setwd("/Users/Jiabo/Documents/gapit/test/Ainong")
setwd("/Users/Jiabo/Documents/gapit/test/Jiabo")

X=myGD[,-1]
simuCV0=GAPIT.PCA(X=X,taxa=myGD[,1],file.output=FALSE,PCA.total=3)
simuCV=simuCV0$PCs[,1:4]


set.seed(99163)
mysimulation<-GAPIT.Phenotype.Simulation(h2=0.75,NQTN=20,GD=myGD,GM=myGM)
posi=mysimulation$QTN.position
myY=mysimulation$Y

index=sample(nrow(myY),50)
# trainGD=myGD[index,]
trainY=myY[,-index]
testY=myY[,index]

myGAPIT<- GAPIT(
Y=trainY,
GD=myGD,#Genotype
GM=myGM,#Map information
PCA.total=0,
CV=simuCV,
# Geno.View.output=FALSE,
# Random.model=FALSE,
# cutOff=0.01,
# QTN.position=mysimulation$QTN.position,
model=c("gBLUP"),
# model=c("MLM"),
file.output=T
)
gblup.pred=merge(testY,myGAPIT$Pred[,c(1,3,5,7,8,10)],by.x="Taxa",by.y="Taxa") 
gblup.r=cor(gblup.pred[,2],gblup.pred[,7])
# source("/Users/Jiabo/Documents/Github/GAPIT/R/GAPIT.R")
# source("/Users/Jiabo/Documents/Github/GAPIT/R/GAPIT.SUPER.GS.R")

myGAPIT2<- GAPIT(
Y=trainY,
GD=myGD,#Genotype
GM=myGM,#Map information
PCA.total=0,
CV=simuCV,
# Geno.View.output=FALSE,
# Random.model=FALSE,
# cutOff=0.01,
# SUPER_GS=TRUE,
# QTN.position=mysimulation$QTN.position,
model=c("sBLUP"),
# model=c("MLM"),
file.output=T
)
sblup.pred=merge(testY,myGAPIT2$Pred[,c(1,3,5,7,8,10)],by.x="Taxa",by.y="Taxa") 
sblup.r=cor(sblup.pred[,2],sblup.pred[,7])
myGAPIT3<- GAPIT(
Y=trainY,
GD=myGD,#Genotype
GM=myGM,#Map information
PCA.total=0,
CV=simuCV,
# Geno.View.output=FALSE,
# Random.model=FALSE,
# cutOff=0.01,
# QTN.position=mysimulation$QTN.position,
model=c("cBLUP"),
# model=c("MLM"),
file.output=T
)
cblup.pred=merge(testY,myGAPIT$Pred[,c(1,3,5,7,8,10)],by.x="Taxa",by.y="Taxa") 
cblup.r=cor(cblup.pred[,2],cblup.pred[,7])





