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
     if(is.null(GD))
     {
       CV=Y[,1:2]
          }else{
       CV=GD[,1:2]
     }
     CV[,2]=1
     colnames(CV)=c("taxa","overall")
     print(paste("There is 0 Covarinces.",sep=""))
     }
     Y=Y[!is.na(Y[,2]),]
     taxa_Y=as.character(Y[,1])
     # print(head(Y))
     if(DP$PCA.total>0&!is.null(DP$CV))CV=GAPIT.CVMergePC(DP$CV,PC)
     if(DP$PCA.total>0&is.null(DP$CV))CV=PC
     if(is.null(GD)&!is.null(DP$KI))
     {
     taxa_KI=as.character(DP$KI[,1])
     taxa_CV=as.character(CV[,1])
     taxa_comall=intersect(intersect(taxa_KI,taxa_Y),taxa_CV)
     # print(length(taxa_comall))
     comCV=CV[taxa_CV%in%taxa_comall,]
     comCV <- comCV[match(taxa_comall,as.character(comCV[,1])),]
     comY=Y[taxa_Y%in%taxa_comall,]
     comY <- comY[match(taxa_comall,as.character(comY[,1])),]
    
     comGD=NULL

     }else{
     # print("@@@@")
     taxa_GD=as.character(GD[,1])
     taxa_comGD=as.character(GD[,1])
     taxa_CV=as.character(CV[,1])
     taxa_comall=intersect(intersect(taxa_GD,taxa_Y),taxa_CV)
     comCV=CV[taxa_CV%in%taxa_comall,]
     comCV <- comCV[match(taxa_comall,as.character(comCV[,1])),]
     
     comGD=GD[taxa_GD%in%taxa_comall,]
     comGD <- comGD[match(taxa_comall,as.character(comGD[,1])),]

     comY=Y[taxa_Y%in%taxa_comall,]
     comY <- comY[match(taxa_comall,as.character(comY[,1])),]
     }


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

