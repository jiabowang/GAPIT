`GAPIT.Judge`<-
function(Y=Y,G=NULL,GD=NULL,KI=NULL,GM=NULL,group.to=group.to,group.from=group.from,sangwich.top=sangwich.top,sangwich.bottom=sangwich.bottom,kinship.algorithm=kinship.algorithm,PCA.total=PCA.total,model="MLM",SNP.test=TRUE){
#Object: To judge Pheno and Geno data practicability
#Designed by Zhiwu Zhang
#Writen by Jiabo Wang
#Last update: Novenber 3, 2016
##############################################################################################
print("--------------------Phenotype and Genotype ----------------------------------")
if(ncol(Y)<2)  stop ("Phenotype should have taxa name and one trait at least. Please correct phenotype file!")
print(kinship.algorithm)
if(is.null(KI)&is.null(GD) & kinship.algorithm!="SUPER"&is.null(G)) stop ("GAPIT says: Kinship is required. As genotype is not provided, kinship can not be created.")
if(kinship.algorithm=="FarmCPU"&SNP.test==FALSE)stop("FarmCPU is only for GWAS, plase set: SNP.test= TRUE")
#if((!is.null(GD))&(!is.null(G))) stop("GAPIT Says:Please put in only one type of geno data.")
if(is.null(GD)&is.null(G)&is.null(KI))stop ("GAPIT Says:GAPIT need genotype!!!")
if(!is.null(GD) & is.null(GM) & (is.null(G)) &SNP.test) stop("GAPIT Says: Genotype data and map files should be in pair")
if(is.null(GD) & !is.null(GM) & (is.null(G)) &SNP.test) stop("GAPIT Says: Genotype data and map files should be in pair")

if(!is.null(GD)&!is.null(Y))
{
  if(nrow(GD)>1){
    if (is.null(GD[,1]%in%Y[,1]))stop("GAPIT Says: There are no common taxa between genotype and phenotype")
     }
}
if(!is.null(G)&!is.null(Y))
{
if (is.null(colnames(G)[-c(1:11)]%in%Y[,1]))stop("GAPIT Says: There are no common taxa between genotype and phenotype")
}

if (!is.null(Y)) nY=nrow(Y)
if (!is.null(Y)) ntrait=ncol(Y)-1
print(paste("There are ",ntrait," traits in phenotype data."))
print(paste("There are ",nY," individuals in phenotype data."))
if (!is.null(G)) nG=nrow(G)-11
if (!is.null(GD)) 
{nG=ncol(GD)-1
print(paste("There are ",nG," markers in genotype data."))}
print("Phenotype and Genotype are test OK !!")

print("--------------------GAPIT Logical ----------------------------------")
#if (group.to>nY&is.null(KI))group.to=nY
#if (group.from>group.to&is.null(KI)) group.from=group.to
if(!is.null(sangwich.top) & is.null(sangwich.bottom) ) stop("GAPIT Says: SUPER method need sangwich.top and bottom")
if(is.null(sangwich.top) & !is.null(sangwich.bottom) ) stop("GAPIT Says: SUPER method need sangwich.top and bottom")
 if(kinship.algorithm=="Separation"&PCA.total==0) stop ("GAPIT Says: Separation kinship need PCA.total>0")







return (list(group.to=group.to,group.from=group.from))
}#end of GAPIT.Pheno.Geno.judge function
#=============================================================================================







