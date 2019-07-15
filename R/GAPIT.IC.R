`GAPIT.IC` <-
function(DP=NULL,CV=NULL){
#Object: To Intermediate Components 
#Designed by Zhiwu Zhang
#Writen by Jiabo Wang
#Last update: Novenber 3, 2016
##############################################################################################
print("GAPIT.IC in process...")

Y=DP$Y
PC=DP$PC
GD=DP$GD

if(DP$kinship.algorithm%in%c("FarmCPU","Blink","MLMM"))
{
#Y=Y[!is.na(Y[,2]),]
taxa_Y=as.character(Y[,1])
taxa_GD=as.character(GD[,1])
#print(dim(PC))
if(!all(taxa_GD%in%taxa_Y))
     {
     com_GD=GD[taxa_GD%in%taxa_Y,]
     if(!is.null(PC))PC=PC[taxa_GD%in%taxa_Y,]
     Y=Y[taxa_Y%in%taxa_GD,]
     
     }else{com_GD=GD
     }
     GT=as.matrix(as.character(com_GD[,1]))
}else{
 GT=as.matrix(as.character(GD[,1]))
 
}

if(DP$PCA.total>0&!is.null(DP$CV))CV=GAPIT.CVMergePC(DP$CV,PC)
if(DP$PCA.total>0&is.null(DP$CV))CV=PC

KI=DP$KI
#print(dim(CV))
#print(is.null(CV))
if (is.null(CV))
{my_allCV=CV
}else{my_allCV=CV[order(CV[,1]),]}
  
noCV=FALSE
if(is.null(CV)){
noCV=TRUE
CV=Y[,1:2]
CV[,2]=1
colnames(CV)=c("taxa","overall")
}
#print(dim(CV))
PCA=CV
K=KI
my_allGD=GD

  print("GAPIT.IC accomplished successfully for multiple traits. Results are saved")
 if(DP$kinship.algorithm%in%c("FarmCPU","Blink","MLMM")){ 
  return (list(Y=Y,GT=GT,PCA=PCA,K=K,GD=com_GD,GM=DP$GM,my_allCV=my_allCV,my_allGD=my_allGD))
}else{
  return (list(Y=Y,GT=GT,PCA=PCA,K=K,GD=DP$GD,GM=DP$GM,my_allCV=my_allCV,my_allGD=my_allGD))
}
}  #end of GAPIT IC function
#=============================================================================================

