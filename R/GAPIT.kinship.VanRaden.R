`GAPIT.kinship.VanRaden` <-
function(snps,hasInbred=TRUE) {
# Object: To calculate the kinship matrix using the method of VanRaden (2009, J. Dairy Sci. 91:4414???C4423)
# Input: snps is n individual rows by m snps columns
# Output: n by n relationship matrix
# Authors: Zhwiu Zhang
# Last update: March 2, 2016 
############################################################################################## 
print("Calculating kinship with VanRaden method...")
#Remove invariants
fa=colSums(snps)/(2*nrow(snps))
index.non=fa>=1| fa<=0
snps=snps[,!index.non]

nSNP=ncol(snps)
nInd=nrow(snps)
n=nInd 

##allele frequency of second allele
p=colSums(snps)/(2*nInd)
P=2*(p-.5) #Difference from .5, multiple by 2
snps=snps-1 #Change from 0/1/2 coding to -1/0/1 coding

print("substracting P...")
Z=t(snps)-P#operation on matrix and vector goes in direction of column
print("Getting X'X...")
#K=tcrossprod((snps), (snps))
K=crossprod((Z), (Z)) #Thanks to Peng Zheng, Meng Huang and Jiafa Chen for finding the problem

print("Adjusting...")
adj=2*sum(p*(1-p))
K=K/adj

print("Calculating kinship with VanRaden method: done")

return(K)
}
#=============================================================================================

