`GAPIT.IC` <-
function(DP=NULL){
#Object: To Intermediate Components 
#Designed by Zhiwu Zhang
#Writen by Jiabo Wang
#Last update: Novenber 3, 2016
##############################################################################################
print("GAPIT.IC in process...")

     Y=DP$Y
     PC=DP$PC
     CV=DP$CV
     GD=DP$GD

     noCV=FALSE
     if(is.null(CV)){
     noCV=TRUE
     CV=Y[,1:2]
     CV[,2]=1
     colnames(CV)=c("taxa","overall")
     print(paste("There is 0 Covarinces.",sep=""))

     }

     taxa_Y=as.character(Y[,1])
     taxa_GD=as.character(GD[,1])

     if(DP$PCA.total>0&!is.null(DP$CV))CV=GAPIT.CVMergePC(DP$CV,PC)
     if(DP$PCA.total>0&is.null(DP$CV))CV=PC

     taxa_comGD=as.character(GD[,1])
     taxa_comY=as.character(Y[,1])
     taxa_CV=as.character(CV[,1])
     taxa_comall=intersect(intersect(taxa_comGD,taxa_comY),taxa_CV)
     comCV=CV[taxa_CV%in%taxa_comall,]
     comGD=GD[taxa_comGD%in%taxa_comall,]
     comY=Y[taxa_comY%in%taxa_comall,]

     GT=as.matrix(as.character(taxa_comall))
     print(paste("There are ",length(GT)," common individuals in genotype , phenotype and CV files.",sep=""))

     if(nrow(comCV)!=length(GT))stop ("GAPIT says: The number of individuals in CV does not match to the number of individuals in genotype files.")

     print("The dimension of total CV is ")
     print(dim(comCV))

     print("GAPIT.IC accomplished successfully for multiple traits. Results are saved")
     if(DP$kinship.algorithm%in%c("FarmCPU","Blink","MLMM")){ 
        return (list(Y=comY,GT=GT,PCA=comCV,K=DP$KI,GD=comGD,GM=DP$GM,myallCV=CV,myallGD=GD))
     }else{
        return (list(Y=comY,GT=GT,PCA=comCV,K=DP$KI,GD=comGD,GM=DP$GM,myallCV=CV,myallGD=GD,myallY=Y))
     }
}  #end of GAPIT IC function
#=============================================================================================

