`GAPIT.Imputation` <-
function(x,GI=NULL,impute="Middle",byRow=TRUE){
#Object: To impute NA in genome
#Output: Coresponding numerical value
#Authors: Zhiwu Zhang
#Writer:  Jiabo Wang
# Last update: April 13, 2016 
##############################################################################################
n=length(x)
lev=levels(as.factor(x))
lev=setdiff(lev,NA)
#print(lev)
len=length(lev)
count=1:len
for(i in 1:len){
	count[i]=length(x[(x==lev[i])])
}
position=order(count)
#print(position)
if(impute=="Middle") {x[is.na(x)]=1 }

if(len==3){
	if(impute=="Minor")  {x[is.na(x)]=position[1]  -1}
	if(impute=="Major")  {x[is.na(x)]=position[len]-1}

}else{
	if(impute=="Minor")  {x[is.na(x)]=2*(position[1]  -1)}
	if(impute=="Major")  {x[is.na(x)]=2*(position[len]-1)}
}

if(byRow) {
  result=matrix(x,n,1)
}else{
  result=matrix(x,1,n)  
}


return(result)
}#end of GAPIT.Imputation function
#=============================================================================================







