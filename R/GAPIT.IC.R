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
     # print(dim(CV))
     if(is.null(CV))
     {
       noCV=TRUE
       if(is.null(GD))
       {
         CV=Y[,1:2]
       }else{
         CV=GD[,1:2]
       }
       CV[,2]=1
       colnames(CV)=c("taxa","overall")
       print(paste("There is 0 Covariances.",sep=""))
     }
     Y=Y[!is.na(Y[,2]),]
     taxa_Y=as.character(Y[,1])
     # print(head(PC))
     # print(head(DP$CV))
     if(DP$PCA.total>0&!is.null(DP$CV))CV=GAPIT.CVMergePC(DP$CV,PC)
     if(DP$PCA.total>0&is.null(DP$CV))CV=PC
     if(!is.null(GD))
     {
       if(!is.null(DP$KI))
       {
         taxa_GD=as.character(GD[,1])
         taxa_KI=as.character(DP$KI[,1])
         taxa_CV=as.character(CV[,1])
         # print(Y)
         # print(tail(taxa_Y))
         # print(tail(taxa_GD))
         # print(tail(taxa_CV))

         taxa_comall=intersect(intersect(intersect(taxa_KI,taxa_GD),taxa_Y),taxa_CV)
         taxa_g_cv=intersect(intersect(taxa_KI,taxa_GD),taxa_CV)
     # print(length(taxa_comall))
     # print(length(taxa_g_cv))
         comCV=CV[taxa_CV%in%taxa_comall,]
         comCV <- comCV[match(taxa_comall,as.character(comCV[,1])),]
         comY=Y[taxa_Y%in%taxa_comall,]
         comY <- comY[match(taxa_comall,as.character(comY[,1])),]
         comGD=GD[taxa_GD%in%taxa_comall,]
         comGD <- comGD[match(taxa_comall,as.character(comGD[,1])),]# comGD=NULL
       }else{
     # print("@@@@")
         taxa_GD=as.character(GD[,1])
         taxa_comGD=as.character(GD[,1])
         taxa_CV=as.character(CV[,1])
         taxa_g_cv=intersect(taxa_GD,taxa_CV)
         taxa_comall=intersect(intersect(taxa_GD,taxa_Y),taxa_CV)
         comCV=CV[taxa_CV%in%taxa_comall,]
         comCV <- comCV[match(taxa_comall,as.character(comCV[,1])),]
         comGD=GD[taxa_GD%in%taxa_comall,]
         comGD <- comGD[match(taxa_comall,as.character(comGD[,1])),]
         comY=Y[taxa_Y%in%taxa_comall,]
         comY <- comY[match(taxa_comall,as.character(comY[,1])),]
       }
       GD=GD[taxa_GD%in%taxa_g_cv,]
       GD=GD[match(taxa_g_cv,as.character(GD[,1])),]
     }else{
       # taxa_GD=as.character(GD[,1])
       if(!is.null(DP$KI))
       {
        taxa_KI=as.character(DP$KI[,1])
        taxa_CV=as.character(CV[,1])
        taxa_comall=intersect(intersect(taxa_KI,taxa_Y),taxa_CV)
        taxa_g_cv=intersect(taxa_KI,taxa_CV)
     # print(length(taxa_comall))
        comCV=CV[taxa_CV%in%taxa_comall,]
        comCV <- comCV[match(taxa_comall,as.character(comCV[,1])),]
        comY=Y[taxa_Y%in%taxa_comall,]
        comY <- comY[match(taxa_comall,as.character(comY[,1])),]
        comGD=NULL
        }else{
        # taxa_KI=as.character(DP$KI[,1])
        taxa_CV=as.character(CV[,1])
        taxa_g_cv=taxa_CV
        taxa_comall=intersect(taxa_Y,taxa_CV)
     # print(length(taxa_comall))
        comCV=CV[taxa_CV%in%taxa_comall,]
        comCV <- comCV[match(taxa_comall,as.character(comCV[,1])),]
        comY=Y[taxa_Y%in%taxa_comall,]
        comY <- comY[match(taxa_comall,as.character(comY[,1])),]
        DP$KI=cbind(as.character(taxa_comall),as.data.frame(matrix(rnorm(length(taxa_comall)^2),length(taxa_comall),length(taxa_comall))))
        colnames(DP$KI)=c("taxa",as.character(taxa_comall)[-1])
        comGD=NULL
        }#end of K
     }# end of GD
     # print(DP$KI[1:5,1:5])
         # print(tail(comY[,1]))
         # print(tail(comGD[,1]))
         # print(tail(comCV[,1]))
         # print(tail(GD[,1]))
     GT=as.matrix(as.character(taxa_comall))
     print(paste("There are ",length(GT)," common individuals in genotype , phenotype and CV files.",sep=""))
     if(nrow(comCV)!=length(GT))stop ("GAPIT says: The number of individuals in CV does not match to the number of individuals in genotype files.")
     print("The dimension of total CV is ")
     print(dim(comCV))
     print(dim(CV))
     CV=CV[taxa_CV%in%taxa_g_cv,]
     CV=CV[match(taxa_g_cv,as.character(CV[,1])),]
     # print(head(CV))
     print("GAPIT.IC accomplished successfully for multiple traits. Results are saved")
     if(DP$kinship.algorithm%in%c("FarmCPU","BLINK","MLMM")){ 
        return (list(Y=comY,GT=GT,PCA=comCV,KI=DP$KI,GD=comGD,GM=DP$GM,myallCV=CV,myallGD=GD))
     }else{
        return (list(Y=comY,GT=GT,PCA=comCV,KI=DP$KI,GD=comGD,GM=DP$GM,myallCV=CV,myallGD=GD,myallY=Y))
     }
}  #end of GAPIT IC function
#=============================================================================================

