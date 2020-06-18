`GAPIT.QC2` <-
function(Y=NULL,KI=NULL,GT=NULL,CV=NULL,Z=NULL,GK=NULL){
#Object: to do data quality control
#Output: Y, KI, GD, CV, Z, flag
#Authors: Zhiwu Zhang and Alex Lipka 
# Last update: April 14, 2011 
##############################################################################################
#Remove duplicates 
print("Removing duplicates...")
#print(dim(CV))
Y=GAPIT.RemoveDuplicate(Y)
CV=GAPIT.RemoveDuplicate(CV)
GK=GAPIT.RemoveDuplicate(GK)
if(!is.null(Z))Z=GAPIT.RemoveDuplicate(Z)

#Remove missing phenotype
print("Removing NaN...")
Y=Y[which(Y[,2]!="NaN"),]
Y=Y[!is.na(Y[,2]),]

# Remove duplicates for GT 
# GT row wise, Z column wise, and KI both direction.
print("Remove duplicates for GT...")
#print(dim(GT))
if(!is.null(GT))
{ 
  if(is.null(dim(GT)))taxa.kept=unique(GT)
  if(!is.null(dim(GT)))taxa.kept=unique(GT[,1])

}else{
  taxa.kept=unique(Y[,1])
}

# Remove duplicates for KI 
print("Remove duplicates for KI...")
# improve speed: remove t() and use cbind
if(!is.null(KI))
{
  taxa.all=KI[,1]
  taxa.uniqe=unique(taxa.all)
  position=match(taxa.uniqe, taxa.all,nomatch = 0)
  position.addition=cbind(1,t(1+position))
  KI=KI[position,position.addition]
}

#Sort KI
# if(!is.null(KI))
# {
#   taxa.all=KI[,1]
#   position=order(taxa.all)
#   position.addition=cbind(1,t(1+position))
#   KI=KI[position,position.addition]
# }

# Remove duplicates for Z rowwise
print("Remove duplicates for Z (column wise)...")
if(!is.null(Z))
{
  taxa.all=as.matrix(Z[1,])
  taxa.uniqe=intersect(taxa.all,taxa.all)
  position=match(taxa.uniqe, taxa.all,nomatch = 0)
  Z=Z[,position]
}


#Remove the columns of Z if they are not in KI/GT. KI/GT are allowed to have individuals not in Z
print("Maching Z with Kinship colwise...")
if(!is.null(KI))
{
  taxa.all=as.character(KI[,1])
  taxa.kinship=unique(taxa.all)
}

if(!is.null(Z) & !is.null(KI))
{
  #get common taxe between KI and Z
  taxa.Z=as.character(as.matrix(Z[1,-1]))
  #taxa.Z=colnames(Z) #This does not work for names starting with numerical or "-"   \
  if(is.null(KI)){
  taxa.Z_K_common=taxa.Z
  }else{
  taxa.Z_K_common=intersect(taxa.kinship,taxa.Z)
  }
  # print(taxa.kinship)
  # print(taxa.Z)
  if(!identical(taxa.kinship,taxa.Z))Z <-cbind(Z[,1], Z[,1+match(taxa.Z_K_common, taxa.Z, nomatch = 0)])
  
  #Remove the rows of Z if all the ellements sum to 0
  #@@@ improve speed: too many Zs
  print("Maching Z without origin...")
  Z1=Z[-1,-1]
  # Z2=data.frame(Z1) # by jiabo
  # Z3=as.matrix(Z2)
  # Z4=as.numeric(Z3) #one dimemtion
  # Z5=matrix(data = Z4, nrow = nrow(Z1), ncol = ncol(Z1))
  Z5=as.matrix(Z1)
  RS=apply(Z5,1,function(one) sum(as.numeric(one)))>0
  #The above process could be simplified!
  Z <- Z[c(TRUE,RS),]
  
  #make individuals the same in Z, Y, GT and CV
  print("Maching GT and CV...")
  if(length(Z)<=1)stop("GAPIT says: there is no place to match IDs!")
}# end of  if(!is.null(Z) & !is.null(K))

# get intersect of all the data
taxa=intersect(Y[,1],Y[,1])
if(!is.null(Z))taxa=intersect(Z[-1,1],taxa)
if(!is.null(GT))taxa=intersect(taxa,taxa.kept)
if(!is.null(CV))taxa=intersect(taxa,CV[,1])
if(!is.null(GK))taxa=intersect(taxa,GK[,1])
if(length(taxa)<=1)stop("GAPIT says: There is no individual ID matched to covariate. Please check!")

print(length(taxa))
# print(taxa)
# if(!is.null(Z))
# {
#   #Remove taxa in Z that are not in others, columnwise
#   t=c(TRUE, taxa.Z%in%taxa)
#   if(length(t)<=2)stop("GAPIT says: There is no individual ID matched among data. Please check!")
#   Z <- Z[t,]
  
#   #Remove the columns of Z if all the ellements sum to 0
#   print("QC final process...")
#   #@@@ improve speed: too many Zs
#   Z1=Z[-1,-1]
#   Z2=data.frame(Z1)
#   Z3=as.matrix(Z2)
#   Z4=as.numeric(Z3) #one dimemtion
#   Z5=matrix(data = Z4, nrow = nrow(Z1), ncol = ncol(Z1))
#   CS=colSums(Z5)>0
#   #The above process could be simplified!
#   Z <- Z[,c(TRUE,CS)]
# }
# print(cbind(Y[,1],CV[,1]))
#Filtering with comman taxa
# print(head(Y))
Y <- Y[Y[,1]%in%taxa,]
Y <- GAPIT.CVMergePC(as.matrix(cbind(taxa,1)),Y)[,-c(2)]
# print(head(Y))
if(!is.null(CV)) 
  {
  CV=CV[CV[,1]%in%taxa,]
  CV <- GAPIT.CVMergePC(cbind(taxa,1),CV)[,-c(2)]
  }
if(!is.null(GK)) 
  {
  GK=GK[GK[,1]%in%taxa,]
  GK <- GAPIT.CVMergePC(cbind(taxa,1),GK)[,-c(2)]
  }
if(!is.null(GT)) taxa.kept=data.frame(taxa.kept[taxa.kept%in%taxa])

#Y <- Y[Y[,1]%in%taxa.kept,]
# print(cbind(as.character(Y[,1]),as.character(CV[,1])))

#To sort Y, GT, CV and Z
# Y=Y[order(Y[,1]),]
# CV=CV[order(CV[,1]),]
# if(!is.null(GK))GK=GK[order(GK[,1]),]
# if(!is.null(Z))Z=Z[c(1,1+order(Z[-1,1])),]

#get position of taxa.kept in GT
#position=match(taxa.kept[,1], GT[,1],nomatch = 0)
if(is.null(dim(GT)))position=match(taxa.kept, GT,nomatch = 0)
if(!is.null(dim(GT)))position=match(taxa.kept[,1], GT[,1],nomatch = 0)


if(is.null(dim(taxa.kept)))order.taxa.kept=order(taxa.kept)
if(!is.null(dim(taxa.kept)))order.taxa.kept=order(taxa.kept[,1])

GTindex=position[order.taxa.kept]
flag=nrow(Y)==nrow(Z)-1&nrow(Y)==nrow(GT)&nrow(Y)==nrow(CV)

print("GAPIT.QC accomplished successfully!")

#print(dim(Y))
#print(dim(CV))
#print(dim(KI))
return(list(Y = Y, KI = KI, GT = GT, CV = CV, Z = Z, GK = GK, GTindex=GTindex, flag=flag))
}#The function GAPIT.QC ends here
#=============================================================================================

`GAPIT.QC` <-
function(Y=NULL,KI=NULL,GT=NULL,CV=NULL,Z=NULL,GK=NULL){
#Object: to do data quality control
#Output: Y, KI, GD, CV, Z, flag
#Authors: Zhiwu Zhang and Alex Lipka 
# Last update: April 14, 2011 
##############################################################################################
#Remove duplicates 
print("Removing duplicates...")
#print(dim(CV))
Y=GAPIT.RemoveDuplicate(Y)
CV=GAPIT.RemoveDuplicate(CV)
GK=GAPIT.RemoveDuplicate(GK)
if(!is.null(Z))Z=GAPIT.RemoveDuplicate(Z)

#Remove missing phenotype
print("Removing NaN...")
Y=Y[which(Y[,2]!="NaN"),]

# Remove duplicates for GT 
# GT row wise, Z column wise, and KI both direction.
print("Remove duplicates for GT...")
#print(dim(GT))
if(!is.null(GT))
{ 
  if(is.null(dim(GT)))taxa.kept=unique(GT)
  if(!is.null(dim(GT)))taxa.kept=unique(GT[,1])

}else{
  taxa.kept=unique(Y[,1])
}

# Remove duplicates for KI 
print("Remove duplicates for KI...")
# improve speed: remove t() and use cbind
if(!is.null(KI))
{
  taxa.all=KI[,1]
  taxa.uniqe=unique(taxa.all)
  position=match(taxa.uniqe, taxa.all,nomatch = 0)
  position.addition=cbind(1,t(1+position))
  KI=KI[position,position.addition]
}

#Sort KI
if(!is.null(KI))
{
  taxa.all=KI[,1]
  position=order(taxa.all)
  position.addition=cbind(1,t(1+position))
  KI=KI[position,position.addition]
}

# Remove duplicates for Z rowwise
print("Remove duplicates for Z (column wise)...")
if(!is.null(Z))
{
  taxa.all=as.matrix(Z[1,])
  taxa.uniqe=intersect(taxa.all,taxa.all)
  position=match(taxa.uniqe, taxa.all,nomatch = 0)
  Z=Z[,position]
}


#Remove the columns of Z if they are not in KI/GT. KI/GT are allowed to have individuals not in Z
print("Maching Z with Kinship colwise...")
if(!is.null(KI))
{
  taxa.all=KI[,1]
  taxa.kinship=unique(taxa.all)
}

if(!is.null(Z) & !is.null(KI))
{
  #get common taxe between KI and Z
  taxa.Z=as.matrix(Z[1,])
  #taxa.Z=colnames(Z) #This does not work for names starting with numerical or "-"   \
  if(is.null(KI)){
  taxa.Z_K_common=taxa.Z
  }else{
  taxa.Z_K_common=intersect(taxa.kinship,taxa.Z)
  }
  Z <-cbind(Z[,1], Z[,match(taxa.Z_K_common, taxa.Z, nomatch = 0)])
  
  #Remove the rows of Z if all the ellements sum to 0
  #@@@ improve speed: too many Zs
  print("Maching Z without origin...")
  Z1=Z[-1,-1]
  Z2=data.frame(Z1)
  Z3=as.matrix(Z2)
  Z4=as.numeric(Z3) #one dimemtion
  Z5=matrix(data = Z4, nrow = nrow(Z1), ncol = ncol(Z1))
  RS=rowSums(Z5)>0
  #The above process could be simplified!
  Z <- Z[c(TRUE,RS),]
  
  #make individuals the same in Z, Y, GT and CV
  print("Maching GT and CV...")
  if(length(Z)<=1)stop("GAPIT says: there is no place to match IDs!")
}# end of  if(!is.null(Z) & !is.null(K))

# get intersect of all the data
taxa=intersect(Y[,1],Y[,1])
if(!is.null(Z))taxa=intersect(Z[-1,1],taxa)
if(!is.null(GT))taxa=intersect(taxa,taxa.kept)
if(!is.null(CV))taxa=intersect(taxa,CV[,1])
if(!is.null(GK))taxa=intersect(taxa,GK[,1])
if(length(taxa)<=1)stop("GAPIT says: There is no individual ID matched to covariate. Please check!")


if(!is.null(Z))
{
  #Remove taxa in Z that are not in others, columnwise
  t=c(TRUE, Z[-1,1]%in%taxa)
  if(length(t)<=2)stop("GAPIT says: There is no individual ID matched among data. Please check!")
  Z <- Z[t,]
  
  #Remove the columns of Z if all the ellements sum to 0
  print("QC final process...")
  #@@@ improve speed: too many Zs
  Z1=Z[-1,-1]
  Z2=data.frame(Z1)
  Z3=as.matrix(Z2)
  Z4=as.numeric(Z3) #one dimemtion
  Z5=matrix(data = Z4, nrow = nrow(Z1), ncol = ncol(Z1))
  CS=colSums(Z5)>0
  #The above process could be simplified!
  Z <- Z[,c(TRUE,CS)]
}

#Filtering with comman taxa
Y <- Y[Y[,1]%in%taxa,]
if(!is.null(CV)) CV=CV[CV[,1]%in%taxa,]
if(!is.null(GK)) GK=GK[GK[,1]%in%taxa,]
if(!is.null(GT)) taxa.kept=data.frame(taxa.kept[taxa.kept%in%taxa])
#Y <- Y[Y[,1]%in%taxa.kept,]

#To sort Y, GT, CV and Z
Y=Y[order(Y[,1]),]
CV=CV[order(CV[,1]),]
if(!is.null(GK))GK=GK[order(GK[,1]),]
if(!is.null(Z))Z=Z[c(1,1+order(Z[-1,1])),]

#get position of taxa.kept in GT
#position=match(taxa.kept[,1], GT[,1],nomatch = 0)
if(is.null(dim(GT)))position=match(taxa.kept, GT,nomatch = 0)
if(!is.null(dim(GT)))position=match(taxa.kept[,1], GT[,1],nomatch = 0)


if(is.null(dim(taxa.kept)))order.taxa.kept=order(taxa.kept)
if(!is.null(dim(taxa.kept)))order.taxa.kept=order(taxa.kept[,1])

GTindex=position[order.taxa.kept]
flag=nrow(Y)==nrow(Z)-1&nrow(Y)==nrow(GT)&nrow(Y)==nrow(CV)

print("GAPIT.QC accomplished successfully!")

#print(dim(Y))
#print(dim(CV))
#print(dim(KI))
return(list(Y = Y, KI = KI, GT = GT, CV = CV, Z = Z, GK = GK, GTindex=GTindex, flag=flag))
}#The function GAPIT.QC ends here
#=============================================================================================
