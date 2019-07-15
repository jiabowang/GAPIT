`GAPIT.kinship.separation` <-
function(PCs=NULL,EV=NULL,nPCs=0 ){
#Object: To calculate kinship from PCS
#       PCs: the principal component as columns and individual as rows, the first column is taxa
#       EV: Eigen values
#       nPCs: the number of front PCs excluded to calculate kinship
#Output: kinship
#Authors: Huihui Li and Zhiwu Zhang
#Last update: April 17, 2012
##############################################################################################
print("Calling GAPIT.kinship.separation")  
  Total.number.PCs=ncol(PCs)
  n=nrow(PCs)
print(Total.number.PCs)
print(n)
  #Choose Total.number.PCs-nPCs PCs and EV to calculate K
  sep.PCs=PCs[, (nPCs+2):(Total.number.PCs)]  #first column is taxa
  sep.EV=EV[(nPCs+1):Total.number.PCs]

  Weighted.sep.EV=sep.EV/sum(sep.EV)
  
  #X=t(t(sep.PCs)*Weighted.sep.EV)  
  X=sep.PCs
   
  XMean= apply(X,2,mean)
  X=as.matrix(X-XMean)
  K=tcrossprod((X), (X))

  #Extract diagonals
  i =1:n
  j=(i-1)*n
  index=i+j
  d=K[index]
  DL=min(d)
  DU=max(d)
  floor=min(K)
  
  K=(K-floor)/(DL-floor)
  MD=(DU-floor)/(DL-floor)
     
  if(is.na(K[1,1])) stop ("GAPIT says: Missing data is not allowed for numerical genotype data")
  if(MD>2)K[index]=K[index]/(MD-1)+1
print("GAPIT.kinship.separation called succesfuly")
  return (K)
}
#=============================================================================================

