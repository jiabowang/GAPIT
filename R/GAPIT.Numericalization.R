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
   if(length(inter)==2)
   {
     x[x==inter[2]]=inter[1]
     n=length(x)
     lev=levels(as.factor(x))
     lev=setdiff(lev,"N")
     lev=setdiff(lev,"NN")
     inter=inter[1]
     #print(lev)
     len=length(lev)
   }
  
#Genotype counts
count=1:len
for(i in 1:len){
	count[i]=length(x[(x==lev[i])])
}
up=0
down=0
if(Major.allele.zero){
    count.temp = cbind(lev,count)
    up=2
    if(length(inter)!=0)
      {
        if(lev[1]!=inter)
        {
          count.temp = count.temp[-which(lev==inter),,drop=FALSE]
          count=count[-which(lev==inter)]
          lev=lev[-which(lev==inter)]
          len=length(lev)
        }
      }      # if(nrow(count.temp)==0) return()
      order.index=order(as.numeric(count.temp[,2]), decreasing = FALSE)
      count.temp <- count.temp[order.index,]
      count = count[order.index]
      lev = lev[order.index]
}else{
      count.temp = cbind(lev,count)
      down=2
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
      count = count[order.index]
      lev = lev[order.index]
    # print(lev)
} #End  if(Major.allele.zero)

#1status other than 2 or 3
if(len<=1 | len> 3)    x=ifelse(x=="N",NA,ifelse(x==inter,1,ifelse(x==lev[1],down,up))) 
#2 status
if(len==2)
{
  if(!setequal(character(0),inter))
  {
    x=ifelse(x=="N",NA,ifelse(x==inter,1,ifelse(x==lev[1],2,0))) 
  }else{
   x=ifelse(x=="N",NA,ifelse(x==lev[1],2,0))     # the most is set 0, the least is set 2
  }
}

#3 status
if(len==3)
{
  x=ifelse(x=="N",NA,ifelse(x==lev[1],2,ifelse(x==inter,1,0)))
  # if(bit==2)x=ifelse(x=="NN",NA,ifelse(x==lev[1],2,ifelse(x==inter,1,0)))
}

#missing data imputation

if(impute=="Middle") {x[is.na(x)]=1}
if(impute=="Minor")  {x[is.na(x)]=lev[1]}
if(impute=="Major")  {x[is.na(x)]=lev[2]}

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

