`GAPIT.HMP.Amplification`<-function(x
	)
{
#Object: To amplificate HapMap file from MNPs
#Output: New HapMap file without MNPs
#Authors: Jiabo Wang
#Writen by Jiabo Wang
#Last update: Oct 1, 2024
##############################################################################################
# setwd("/Users/Jiabo/Documents/Data/Meiai")
# myG=data.table::fread("gene.hmp.txt",
                          # header = F,
                          # na.strings = c("NA", "NaN"),
                          # data.table = F)
# X=myG[-1,]
muta=x[2]
n.muta=nchar(muta)

if(n.muta==3)
  {
     map=as.character(x[c(1:11)])
     genotype=as.character(x[-c(1:11)])
     genotype[genotype=="A"]="AA"
     genotype[genotype=="C"]="CC"
     genotype[genotype=="T"]="TT"
     genotype[genotype=="G"]="GG"
     genotype[genotype=="R"]="AG"
     genotype[genotype=="Y"]="CT"
     genotype[genotype=="S"]="CG"
     genotype[genotype=="W"]="AT"
     genotype[genotype=="K"]="GT"
     genotype[genotype=="M"]="AC"
     genotype[genotype=="N"]="NN"


     new.X=append(map,genotype)
     new.X=matrix(new.X,1,length(new.X))
     return(as.data.frame(new.X))
  }else{
     # muta2=strsplit(muta,"/")[[1]]
     # genotype=X
     map=as.character(x[c(1:11)])
     genotype=as.character(x[-c(1:11)])
     genotype[genotype=="A"]="AA"
     genotype[genotype=="C"]="CC"
     genotype[genotype=="T"]="TT"
     genotype[genotype=="G"]="GG"
     genotype[genotype=="R"]="AG"
     genotype[genotype=="Y"]="CT"
     genotype[genotype=="S"]="CG"
     genotype[genotype=="W"]="AT"
     genotype[genotype=="K"]="GT"
     genotype[genotype=="M"]="AC"
     genotype[genotype=="N"]="NN"
     genotypes=unlist(strsplit(genotype, ""))
     genotypes2=genotypes[genotypes!="N"]
     lev=levels(as.factor(genotypes2))
     len=length(lev)
     count=1:len
     for(i in 1:len){
	     count[i]=length(genotypes2[(genotypes2==lev[i])])
     }
      count.temp = cbind(lev,count)
      order.index=order(as.numeric(count.temp[,2]), decreasing = FALSE)
      count.temp <- count.temp[order.index,]
      count = count[order.index]
      lev = lev[order.index]
     MAF.lev=lev[1]
     heter=c("AT","AG","AC","TA","GA","CA","GT","TG","GC","CG","CT","TC")
     homo=paste(lev[1],lev[1],sep="")
     homo0=NULL
     for(i in 2:len)
     {
     	homo0=c(homo0,paste(lev[i],lev[i],sep=""))
     }
     new.genotype=matrix(homo,len-1,length(genotype))
     for(i in 2:len)
     {
     heter0=c(paste(lev[i],lev[-c(i)],sep=""),paste(lev[-c(i)],lev[i],sep=""))
# heter2=setdiff(heter,heter0)
     homo1=paste(lev[i],lev[i],sep="")
     homo2=setdiff(homo0,homo1)
# homo2=c(homo2,)
# new.genotype[i-1,]=genotype
     new.genotype[i-1,genotype%in%heter0]=heter0[1]
     new.genotype[i-1,genotype%in%homo2]=homo
     new.genotype[i-1,genotype%in%homo1]=homo1
     }
     maps=matrix(NA,len-1,11)
     maps[,1]=paste(map[1],"_",1:(len-1),sep="")
     maps[,2]=paste(MAF.lev,"/",lev[-1],sep="")
     maps[,3]=rep(map[3],len-1)
     maps[,4]=as.numeric(map[4])+(1:(len-1))
     new.X=cbind(maps,new.genotype)
     return(as.data.frame(new.X))
  }

}

