`GAPIT.Numericalization` <-
function(x,bit=2,effect="Add",impute="None", Create.indicator = FALSE, Major.allele.zero = FALSE, byRow=TRUE){
#Object: To convert character SNP genotpe to numerical
#Output: Coresponding numerical value
#Authors: Feng Tian and Zhiwu Zhang
# Last update: May 30, 2011 
##############################################################################################
if(bit==1)  {
x[x=="X"]="N"
x[x=="-"]="N"
x[x=="+"]="N"
x[x=="/"]="N"
x[x=="K"]="Z" #K (for GT genotype)is replaced by Z to ensure heterozygose has the largest value
}

if(bit==2)  {
x[x=="XX"]="N"
x[x=="--"]="N"
x[x=="++"]="N"
x[x=="//"]="N"
x[x=="NN"]="N"
x[x=="00"]="N"

}

n=length(x)
lev=levels(as.factor(x))
lev=setdiff(lev,"N")
#print(lev)
len=length(lev)
#print(len)
#Jiabo creat this code to convert AT TT to 1 and 2. 2018.5.29
if(bit==2)
{
   inter_store=c("AT","AG","AC","TA","GA","CA","GT","TG","GC","CG","CT","TC")
   inter=intersect(lev,inter_store)
   if(length(inter)>1)
   {
     x[x==inter[2]]=inter[1]
     n=length(x)
     lev=levels(as.factor(x))
     lev=setdiff(lev,"N")
     #print(lev)
     len=length(lev)
   }
   # if(len==2)
   # { #inter=intersect(lev,inter_store)
   #   if(!setequal(character(0),inter))
   #   { 
   #     lev=union(lev,"UU")
   #     len=len+1
   #   }
   # }
   if(len==3&bit==2)
   {
     inter=intersect(lev,inter_store)
   }
}
#print(lev)
#print(len)
#Jiabo code is end here

#Genotype counts
count=1:len
for(i in 1:len){
	count[i]=length(x[(x==lev[i])])
}

#print(count)

if(Major.allele.zero){
  if(len>1 & len<=3){
    #One bit: Make sure that the SNP with the major allele is on the top, and the SNP with the minor allele is on the second position
    if(bit==1){ 
      count.temp = cbind(count, seq(1:len))
      if(len==3) count.temp = count.temp[-3,]
      count.temp <- count.temp[order(count.temp[,1], decreasing = TRUE),]
      if(len==3)order =  c(count.temp[,2],3)else order = count.temp[,2]
    }
    #Two bit: Make sure that the SNP with the major allele is on the top, and the SNP with the minor allele is on the third position
    if(bit==2){ 
      count.temp = cbind(count, seq(1:len))
      if(len==3) count.temp = count.temp[-2,]
      count.temp <- count.temp[order(count.temp[,1], decreasing = TRUE),]
      if(len==3) order =  c(count.temp[1,2],2,count.temp[2,2])else order = count.temp[,2]
    }

    count = count[order]
    # print(count)
    lev = lev[order]
    # print(lev)

  }   #End  if(len<=1 | len> 3)
} #End  if(Major.allele.zero)

#print(x)

#make two  bit order genotype as AA,AT and TT, one bit as A(AA),T(TT) and X(AT)
if(bit==1 & len==3){
	temp=count[2]
	count[2]=count[3]
	count[3]=temp
}

#print(lev)
#print(count)
position=order(count)

#Jiabo creat this code to convert AT TA to 1 and 2.2018.5.29

# lev1=lev
# if(bit==2&len==3) 
# {
# lev1[1]=lev[count==sort(count)[1]]
# lev1[2]=lev[count==sort(count)[2]]
# lev1[3]=lev[count==sort(count)[3]]
# position=c(1:3)
# lev=lev1
# }
#print(lev)
#print(position)
#print(inter)
#Jiabo code is end here
if(bit==1){
  lev0=c("R","Y","S","W","K","M") 
  inter=intersect(lev,lev0)
}

#1status other than 2 or 3
if(len<=1 | len> 3)x=0

#2 status
if(len==2)
{
  
  if(!setequal(character(0),inter))
  {
    x=ifelse(x=="N",NA,ifelse(x==inter,1,0)) 
    }else{
    x=ifelse(x=="N",NA,ifelse(x==lev[1],0,2))     # the most is set 0, the least is set 2
  }
}

#3 status
if(bit==1){
	if(len==3)x=ifelse(x=="N",NA,ifelse(x==lev[1],0,ifelse(x==lev[3],1,2)))
}else{
	if(len==3)x=ifelse(x=="N",NA,ifelse(x==lev[lev!=inter][1],0,ifelse(x==inter,1,2)))
}

#print(paste(lev,len,sep=" "))
#print(position)

#missing data imputation
if(impute=="Middle") {x[is.na(x)]=1 }

if(len==3){
	if(impute=="Minor")  {x[is.na(x)]=position[1]  -1}
	if(impute=="Major")  {x[is.na(x)]=position[len]-1}

}else{
	if(impute=="Minor")  {x[is.na(x)]=2*(position[1]  -1)}
	if(impute=="Major")  {x[is.na(x)]=2*(position[len]-1)}
}

#alternative genetic models
if(effect=="Dom") x=ifelse(x==1,1,0)
if(effect=="Left") x[x==1]=0
if(effect=="Right") x[x==1]=2

if(byRow) {
  result=matrix(x,n,1)
}else{
  result=matrix(x,1,n)  
}

return(result)
}#end of GAPIT.Numericalization function
#=============================================================================================

