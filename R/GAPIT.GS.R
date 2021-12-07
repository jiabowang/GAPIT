`GAPIT.GS` <-
function(KW,KO,KWO,GAU,UW){
#Object: to derive BLUP for the individuals without phenotype
#UW:BLUP and PEV of ID with phenotyp
#Output: BLUP
#Authors: Zhiwu Zhang 
# Last update: Oct 22, 2015  by Jiabo Wang
##############################################################################################
#print(dim(UW))
UO=try(t(KWO)%*%solve(KW)%*%UW,silent=TRUE)
#print(dim(KWO)) #kinship without inference
#print(dim(KW))  #kinship within inference
# print(dim(UW))  #BLUP AND PEV of reference

if(inherits(UO, "try-error")) 
{
	GTT=try(t(KWO)%*%MASS::ginv(as.matrix(KW))%*%UW)
	if(inherits(GTT,"try-error"))
	{
        utils::write.csv(KW,"KW.csv",quote=F,row.names=F)
        KW=utils::read.csv("KW.csv",head=T)
        UO=t(KWO)%*%MASS::ginv(as.matrix(KW))%*%UW
        system("rm KW.csv")
    }else{
    	UO=GTT
    }
}
n=ncol(UW) #get number of columns, add additional for individual name

BLUP=data.frame(as.matrix(GAU[,1:4]))
BLUP.W=BLUP[which(GAU[,3]<2),]
W_BLUP=BLUP.W[order(as.numeric(as.matrix(BLUP.W[,4]))),]
UW=UW[which(rownames(UW)==colnames(KW)),] # get phenotype groups order

ID.W=as.numeric(as.matrix(W_BLUP[,4]))
n.W=max(ID.W)
DS.W=diag(n.W)[ID.W,]
# print(dim(DS.W))
# print(dim(UW))
ind.W=DS.W%*%UW

all.W=cbind(W_BLUP,ind.W)
all=all.W

BLUP.O=BLUP[which(GAU[,3]==2),]
O_BLUP=BLUP.O[order(as.numeric(as.matrix(BLUP.O[,4]))),]
#print(dim(O_BLUP))
if(nrow(O_BLUP)>0){

ID.O=as.numeric(as.matrix(O_BLUP[,4]))
n.O=max(ID.O)
DS.O=diag(n.O)[ID.O,]
ind.O=DS.O%*%UO
all.O=cbind(O_BLUP,ind.O)
all=rbind(all.W,all.O)
}

colnames(all)=c("Taxa", "Group", "RefInf","ID","BLUP","PEV")

print("GAPIT.GS accomplished successfully!")
return(list(BLUP=all))
}#The function GAPIT.GS ends here
#=============================================================================================

