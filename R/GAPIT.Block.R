`GAPIT.Block` <-
function(Z,GA,KG){
#Object: To split a group kinship into two blocks containing individuals with and without phenotype
#Output: GAU,KW,KO,KWO
#Authors: Zhiwu Zhang and Alex Lipka 
# Last update: April 14, 2011 
##############################################################################################
# To separate group kiship into two blocks: with and without phenotype.
# A group goes to with phenotype as loog as it has one phenotyped individual.

#find position in group assignment (GA) for the individual associate with phenotype (specified by Z)
#taxa=unique(intersect(as.matrix(Z[1,-1]),GA[,1]))

taxa.Z=as.matrix(Z[1,-1])
taxa.GA=as.matrix(GA[,1])
position=taxa.GA%in%taxa.Z

#Initial block as 2
GAU=cbind(GA,2)

#Assign block as 1 if the individual has phenotype
GAU[position,3]=1

#Modify the non-phenotyped individuals if they in a group with phenotyped individuals
#To find the groups with phenotyped individuals
#update block assignment for all these groups
#get list of group that should be block 1

#grp.12=as.matrix(unique(GAU[,2]))
#grp.1=as.matrix(unique(GAU[which(GAU[,3]==1),2]))
#grp.2= as.matrix(setdiff(grp.12,grp.1))

grp.12=as.matrix(as.vector(unique(GAU[,2])) ) #unique group
grp.1=as.matrix(as.vector(unique(GAU[which(GAU[,3]==1),2])) ) #unique phenotyped group
grp.2= as.matrix(as.vector(setdiff(grp.12,grp.1))) #unique unphenotyped group

numWithout=length(grp.2)

order.1=1:length(grp.1)
order.2=1:length(grp.2)
if(numWithout >0) grpblock=as.matrix(rbind(cbind(grp.1,1,order.1), cbind(grp.2,   2,    order.2)))
if(numWithout==0) grpblock=as.matrix(      cbind(grp.1,1,order.1),                       )

order.block=order(as.matrix(GAU[,3]))
colnames(grpblock)=c("grp","block","ID")

#Indicators: 1-Phenotype, 1.5- unphenotyped but in a group with other phenotyped, 2-rest  (Zhiwu, Dec 7,2012)
#GAU0 <- merge(GAU[order.block,-3], grpblock, by.x = "X2", by.y = "grp")
#GAU=GAU0[,c(2,1,3,4)]
#print(head(GAU))
GAU1 <- merge(GAU[order.block,], grpblock, by.x = "X2", by.y = "grp")
#print(GAU1)
GAU1[,4]=(as.numeric(GAU1[,3])+as.numeric(GAU1[,4]))/2
#print(GAU1)

GAU=GAU1[,c(2,1,4,5)]
KW=KG[grp.1,grp.1]
KO=KG[grp.2,grp.2]
KWO=KG[grp.1,grp.2]

#write.table(GAU, "GAU.txt", quote = FALSE, sep = "\t", row.names = TRUE,col.names = TRUE)

#print("GAPIT.Block accomplished successfully!")

return(list(GAU=GAU,KW=KW,KO=KO,KWO=KWO))
}#The function GAPIT.Block ends here
#=============================================================================================

