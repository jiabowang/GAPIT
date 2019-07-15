`GAPIT.HMP2Num` <-
function(nLines=n,fileHMP="hmp.txt",fileNum="num.txt",bit=1,SNP.effect="Add",SNP.impute="Middle",heading=TRUE, Create.indicator = FALSE, Major.allele.zero = FALSE){
    
#Object: To convert hmp file to numerical file
#Input: hmp genotype file
#Output: Numerical genotype file
#Authors: Zhiwu Zhang
# Last update: May 23, 2013 
##############################################################################################
#print("GAPIT.HMP2Num start")

setwd("/Users/Zhiwu/Dropbox/Current/paper/BigData/BUS/Robust/MaizeGBS")
fileHMP="NAMs26HM2.c10.imp.hmp.txt"
fileNum="NAMs26HM2.c10s.imp.num.txt"

bit=1
SNP.effect="Add"
SNP.impute="Middle"
Major.allele.zero = FALSE

system.time({
n=2000
fileHMPCon<-file(fileHMP, open="r")
#fileNumCon<-file(fileNum, open="r")
tt<-readLines(fileHMPCon, n=1) #header
for(i in 1:n){
  if(i %% 100 == 0)print(i)
  tt<-readLines(fileHMPCon, n=1) 
  #tt2<-na.omit(as.numeric(unlist(strsplit(tt, "\t")))) 
  tt2<-unlist(strsplit(tt, "\t"))
  #GM
  rs=tt2[1]
  chrom=tt2[3]
  pos=tt2[4]
  #GD
  GD= GAPIT.Numericalization(x=tt2[-c(1:11)],bit=bit,effect=SNP.effect,impute=SNP.impute, Major.allele.zero=Major.allele.zero)
  
  #Output
  #print(i)
  #print(tt2[12:52]) 
  #print(GD[1:41]) 
  #writeLines(tt2, fileNumCon,append=TRUE)
 
}
close.connection(fileHMPCon)
})

#print("GAPIT.HMP2Num accomplished successfully!")
}   #GAPIT.HMP2Num ends here
#=============================================================================================

