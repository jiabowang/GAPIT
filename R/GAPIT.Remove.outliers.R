`GAPIT.Remove.outliers`=function(x,na.rm=TRUE,pro=0.25,size=1.5,...){
# Remove outliers of phenotype, and set them as max values
#
#    
#
## Input:
#  x: The given phenotype vector
#  
#
## Output:
#  y: The removed phenetype vector
#  idx:  The indices of the removed 
#Authors: Jiabo Wang
#Writer:  Jiabo Wang
# Last update: MAY 12, 2022 
##############################################################################################

	qnt=quantile(x,probs=c(pro,1-pro),na.rm=na.rm,...)
	y=x
	H=size*IQR(y,na.rm=na.rm)
	y[x<=(qnt[1]-H)]=min(y,na.rm=na.rm)
	y[x>=(qnt[2]+H)]=max(y,na.rm=na.rm)
	idx=x<=(qnt[1]-H)|x>=(qnt[2]+H)
	res <- vector("list")
    res$y=y
    res$idx=idx
    return(res)	
}
