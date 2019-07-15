`GAPIT.ZmatrixCompress` <-
function(Z,GAU){
#Object: To assign the fraction of a individual belonging to a group
#Output: Z
#Authors: Zhiwu Zhang
# Last update: April 14, 2011 
##############################################################################################
#Extraction of GAU coresponding to Z, sort GAU rowwise to mach columns of Z, and make design matrix
#print("GAPIT.ZmatrixCompress")
#print(dim(Z))
#print(dim(GAU))

effect.Z=as.matrix(Z[1,-1])
effect.GAU=as.matrix(GAU[,1])
taxa=as.data.frame(Z[-1,1])

GAU0=GAU[effect.GAU%in%effect.Z,]
order.GAU=order(GAU0[,1])
GAU1 <- GAU0[order.GAU,]
#id.1=GAU1[which(GAU1[,3]==1),4]
id.1=GAU1[which(GAU1[,3]<2),4]
n=max(as.numeric(as.vector(id.1)))
x=as.numeric(as.matrix(GAU1[,4]))
DS=diag(n)[x,]

#sort Z column wise
order.Z=order(effect.Z)
Z=Z[-1,-1]
Z <- Z[,order.Z]

#Z matrix from individual to group
#Z1.numeric <- as.numeric(as.matrix(Z))
Z <- matrix(as.numeric(as.matrix(Z)), nrow = nrow(Z), ncol = ncol(Z)) 
Z=Z%*%DS

#Z3=data.frame(cbind(as.character(Z[-1,1]),Z2))
Z=data.frame(cbind(taxa,Z))

#Z=Z3[order(Z3[,1]),]

Z=Z[order(as.matrix(taxa)),]


#print("GAPIT.ZmatrixCompress accomplished successfully!")
return(list(Z=Z))
}#The function GAPIT.ZmatrixCompress ends here
#=============================================================================================

