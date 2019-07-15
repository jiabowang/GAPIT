`GAPIT.Log` <-
function(Y=Y,KI=KI,Z=Z,CV=CV,SNP.P3D=SNP.P3D,
				group.from = group.from ,group.to =group.to ,group.by = group.by ,kinship.cluster = kinship.cluster, kinship.group= kinship.group,
                      	ngrid = ngrid , llin = llin , ulim = ulim , esp = esp ,name.of.trait = name.of.trait){
#Object: To report model factors
#Output: Text file (GAPIT.Log.txt)
#Authors: Zhiwu Zhang
# Last update: may 16, 2011 
##############################################################################################

#Creat storage
facto <- list(NULL)
value <- list(NULL)

#collecting model factors

facto[[1]]="Trait"
value[[1]]=paste(dim(Y))

facto[[2]]="group.by "
value[[2]]=group.by 

facto[[3]]="Trait name "
value[[3]]=name.of.trait

facto[[4]]="Kinship"
value[[4]]=dim(KI)

facto[[5]]="Z Matrix"
value[[5]]=dim(Z)

facto[[6]]="Covariate"
value[[6]]=dim(CV)

facto[[7]]="SNP.P3D"
value[[7]]=SNP.P3D

facto[[8]]="Clustering algorithms"
value[[8]]=kinship.cluster

facto[[9]]="Group kinship"
value[[9]]=kinship.group

facto[[10]]="group.from "
value[[10]]=group.from 

facto[[11]]="group.to "
value[[11]]=group.to 



theLog=as.matrix(cbind(facto,value))
#theLog=as.character(as.matrix(cbind(facto,value)))
colnames(theLog)=c("Model", "Value")
file=paste("GAPIT.", name.of.trait,".Log.csv" ,sep = "")
write.table(theLog, file, quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)

return (theLog)
}
#=============================================================================================

