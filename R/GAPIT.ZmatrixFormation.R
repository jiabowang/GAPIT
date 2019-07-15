`GAPIT.ZmatrixFormation` <-
function(Z,Y){
#Object: To expande the proportion Z to final Z
#Output: Z
#Authors: Zhiwu Zhang 
# Last update: April 22, 2011 
##############################################################################################
#split individuals in Y to the ones that are given Z and the one not
taxa.Z=as.matrix(Z[-1,1])
taxa.Y=as.matrix(Y[,1])
taxa.diff=setdiff(taxa.Y,taxa.Z)
taxa.I=as.matrix(taxa.Y[match(taxa.diff,taxa.Y,nomatch = 0)])
taxa.Z.col=as.matrix(Z[1,-1])

#Create final Z with zero block and identity block
Z0=matrix(data=0,nrow=nrow(taxa.Z),ncol=nrow(taxa.I))
Z1=diag(1,nrow(taxa.I))
ZC=as.matrix(rbind(Z0,Z1))

#To label rows and columns
label.row=rbind(as.matrix(Z[,1]),taxa.I)
label.col=t(taxa.I)

#update the zero block by the given Z matrix
position=t(as.matrix(match(taxa.Z.col,taxa.I,nomatch = 0)))
ZC[1:nrow(taxa.Z),position]=as.matrix(Z[-1,-1])

#habdler of parents do not have phenotype (colums of Z are not in taxa.I)
# To do list

#To form final Z matrix
dataPart=rbind(label.col,ZC)
Z=data.frame(cbind(label.row,dataPart))

#print("GAPIT.ZmatrixFormation accomplished successfully!")
return(Z)
}#The function GAPIT.ZmatrixFormation ends here
#=============================================================================================

