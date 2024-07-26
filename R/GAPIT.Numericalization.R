`GAPIT.Numericalization` <-
function(x,bit=2,effect="Add",impute="Middle", Create.indicator = FALSE, Major.allele.zero = FALSE, byRow=TRUE){
#Object: To convert character SNP genotpe to numerical
#Output: Coresponding numerical value
#Authors: Feng Tian and Zhiwu Zhang
# Last update: May 30, 2011 
##############################################################################################
if(bit==1)  {
x[x=="X"]="N"
# x[x=="-"]="N"
x[x=="+"]="N"
x[x=="/"]="N"
x[x=="K"]="Z" #K (for GT genotype)is replaced by Z to ensure heterozygose has the largest value
}

if(bit==2)  {
x[x=="XX"]="N"
# x[x=="--"]="N"
x[x=="++"]="N"
x[x=="//"]="N"
x[x=="NN"]="N"
x[x=="00"]="N"
# x[x=="-"]="N"

}

n=length(x)
lev=levels(as.factor(x))
lev=setdiff(lev,"N")
lev=setdiff(lev,"NN")
#print(lev)
len=length(lev)
#print(len)
#Jiabo creat this code to convert AT TT to 1 and 2. 2018.5.29

   if(bit==2)inter_store=c("AT","AG","AC","TA","GA","CA","GT","TG","GC","CG","CT","TC","A-","-A","C-","-C","G-","-G","G-","-G")
   if(bit==1)inter_store=c("R","Y","S","W","K","M") 
   inter=intersect(lev,inter_store)
#else{
  #inter=lev
#}#
   if(length(inter)==2)
   {
     x[x==inter[2]]=inter[1]
     n=length(x)
     lev=levels(as.factor(x))
     lev=setdiff(lev,"N")
     lev=setdiff(lev,"NN")
     #print(lev)
     len=length(lev)
   }
  
#Genotype counts
count=1:len
for(i in 1:len){
	count[i]=length(x[(x==lev[i])])
}

# print(count)
# print(len)
# print()

if(Major.allele.zero){
  # if(len>1 & len<=3){
    #One bit: Make sure that the SNP with the major allele is on the top, and the SNP with the minor allele is on the second position
   # if(bit==1){ 
      count.temp = cbind(lev,count)
      if(length(inter)!=0&lev[1]!=inter)count.temp = count.temp[-which(lev==inter),,drop=FALSE]
      # if(nrow(count.temp)==0) return()
      # print("!!!")
      order.index=order(as.numeric(count.temp[,2]), decreasing = FALSE)
      count.temp <- count.temp[order.index,]
      # if(len==3)order =  c(count.temp[,2],3)else order = count.temp[,2]
    # }
    #Two bit: Make sure that the SNP with the major allele is on the top, and the SNP with the minor allele is on the third position
    # if(bit==2){ 
    #   count.temp = cbind(count, seq(1:len))
    #   if(len==3) count.temp = count.temp[-2,]
    #   count.temp <- count.temp[order(count.temp[,1], decreasing = FALSE),]
    #   if(len==3) order =  c(count.temp[1,2],2,count.temp[2,2])else order = count.temp[,2]
    # }

    count = count[order.index]
    # print(count)
    lev = lev[order.index]
    # print(lev)

}else{
    # if(bit==1){ 
      count.temp = cbind(lev,count)
      print(count.temp)
      if(length(inter)!=0)
      {
        if(lev[1]!=inter)
        {
          count.temp = count.temp[-which(lev==inter),,drop=FALSE]
          count=count[-which(lev==inter)]
          lev=lev[-which(lev==inter)]
          len=length(lev)
        }
      }
      order.index=order(as.numeric(count.temp[,2]), decreasing = TRUE)
      count.temp <- count.temp[order.index,]
      # if(len==3)order =  c(count.temp[,2],3)else order = count.temp[,2]
    # }
    #Two bit: Make sure that the SNP with the major allele is on the top, and the SNP with the minor allele is on the third position
    # if(bit==2){ 
    #   count.temp = cbind(count, seq(1:len))
    #   if(len==3) count.temp = count.temp[-2,]
    #   count.temp <- count.temp[order(count.temp[,1], decreasing = FALSE),]
    #   if(len==3) order =  c(count.temp[1,2],2,count.temp[2,2])else order = count.temp[,2]
    # }

    count = count[order.index]
    # print(count)
    lev = lev[order.index]
    # print(lev)

} #End  if(Major.allele.zero)

#1status other than 2 or 3
if(len<=1 | len> 3)x=0

#2 status
if(len==2)
{
  if(!setequal(character(0),inter))
  {
    x=ifelse(x=="N",NA,ifelse(x==inter,1,2)) 
    # if(bit==2)x=ifelse(x=="NN",NA,ifelse(x==inter,1,2)) 
  }else{
   x=ifelse(x=="N",NA,ifelse(x==lev[1],2,0))     # the most is set 0, the least is set 2
    # if(bit==2)x=ifelse(x=="NN",NA,ifelse(x==lev[1],2,0))     # the most is set 0, the least is set 2
  }
}

#3 status
# print(table(x))
if(len==3)
{
  x=ifelse(x=="N",NA,ifelse(x==lev[1],2,ifelse(x==inter,1,0)))
  # if(bit==2)x=ifelse(x=="NN",NA,ifelse(x==lev[1],2,ifelse(x==inter,1,0)))
}
# print(table(x))
#print(paste(lev,len,sep=" "))
#print(position)

#missing data imputation

if(impute=="Middle") {x[is.na(x)]=1}

# if(len==3){
	if(impute=="Minor")  {x[is.na(x)]=lev[1]}
	if(impute=="Major")  {x[is.na(x)]=lev[2]}

# }else{
# 	if(impute=="Minor")  {x[is.na(x)]=2*(position[1]  -1)}
# 	if(impute=="Major")  {x[is.na(x)]=2*(position[len]-1)}
# }

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

