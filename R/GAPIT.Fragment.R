`GAPIT.Fragment` <-
function(file.path=NULL,file.from=NULL, file.to=NULL,file.total=NULL,file.G=NULL,
                          file.Ext.G=NULL,seed=123,SNP.fraction=1,SNP.effect="Add",SNP.impute="Middle",
                          genoFormat=NULL, file.GD=NULL, file.Ext.GD=NULL, file.GM=NULL, file.Ext.GM=NULL, file.fragment=NULL,
                          file=1,frag=1,LD.chromosome=NULL,LD.location=NULL,LD.range=NULL, Create.indicator = FALSE, Major.allele.zero = FALSE){
#Object: To load SNPs on a (frag)ment in file (this is to replace sampler)
#Output: genotype data sampled
#Authors: Alex Lipka and Zhiwu Zhang
# Last update: August 18, 2011
##############################################################################################
#print("Fragmental reading...")
genoFormat="hapmap"
if(!is.null(file.GD)&is.null(file.G)) genoFormat="EMMA"
  
if(genoFormat=="hapmap"){
        #Initical G
        #print("Reading file...")
        G=NULL
        if(frag==1){
          skip.1=0
          G <- try(utils::read.delim(paste(file.path,file.G,file, ".",file.Ext.G,sep=""),
                          head = FALSE,skip = skip.1, nrows = file.fragment+1),silent=TRUE)
        }else{
          skip.1 <- (frag-1)*file.fragment +1
          G <- try(utils::read.delim(paste(file.path,file.G,file, ".",file.Ext.G,sep=""),
                          head = FALSE,skip = skip.1, nrows = file.fragment),silent=TRUE )
        }
        
        #print("processing the data...")
        if(inherits(G, "try-error"))  {
          G=NULL
          #print("File end reached for G!!!")
        }

        if(is.null(G)){
        #print("The above error indicating reading after end of file (It is OK).")
        return(list(GD=NULL,GI=NULL,GT=NULL,linesRead=NULL,GLD=NULL,heading=NULL) )
        }

        #print("Calling hapmap...")
        heading=(frag==1)
        
        #Recording number of lineas read
        if(heading){
          n= nrow(G)-1
        }else{
          n= nrow(G)
        } 
       
       linesRead=n
               
        #Sampling
       if(SNP.fraction<1){

          #print("Number of SNP in this pragment:")
          #print(n)
          
          #set.seed(seed+(file*1000)+frag)
          #mySample=sample(1:n,max(2,floor(n*as.numeric(as.vector(SNP.fraction)))))
          mySample=sample(1:n,max(2,floor(n*SNP.fraction)))
          #print("@@@@@@@@@@")
          #print(mySample)
          #print(length(mySample))
          if(heading){
            G=G[c(1,(1+mySample)),]
          }else{
            G=G[mySample,]
          }
        } #end of if(SNP.fraction<1)
        

        print("Call hapmap from fragment")      
        hm=GAPIT.HapMap(G,SNP.effect=SNP.effect,SNP.impute=SNP.impute,heading=heading, Create.indicator = Create.indicator, Major.allele.zero = Major.allele.zero)

        #print("Extracting snps for LD plot...")
        #Extract SNPs for LD plot
        if(!is.null(LD.chromosome) & !is.null(hm$GD)){
          index=(G[,3]==LD.chromosome[1]) & abs((as.numeric(G[,4])-as.numeric(LD.location[1]))<(as.numeric(LD.range[1])/2))   
          GLD=G[index,]
        }else{
          GLD=NULL
        }
        
        #rm(G)
        #gc()
        print("hapmap called successfuly from fragment")

        return(list(GD=hm$GD,GI=hm$GI,GT=hm$GT,linesRead=linesRead,GLD=GLD,heading=heading,G=G))

          print("ERROR: It should not get here!!!")        
} #end of "hapmap"



if(genoFormat=="EMMA"){
#print("The file is a numerical format!")
        #Initial GD
        GD=NULL
        skip.1 <- (frag-1)*file.fragment
        #Skip the remaining columns
        GD.temp <- try(utils::read.table(paste(file.path,file.GD, file, ".", file.Ext.GD,sep=""), head = TRUE, nrows = 1),silent=TRUE)
        num.SNP <- ncol(GD.temp)-1
        rm(GD.temp)
        read.in <- min(file.fragment,(num.SNP-skip.1))
        skip.2 <- max((num.SNP - (skip.1 + read.in)),0)
        print(paste(file.path,file.GD,file, ".",file.Ext.GD,sep=""))

        GD <- try(utils::read.table(paste(file.path,file.GD,file, ".",file.Ext.GD,sep=""), head = TRUE,
                  colClasses = c("factor", rep("NULL", skip.1), rep("numeric", read.in),
                  rep("NULL", skip.2))) ,silent=TRUE)
        GI <- try(utils::read.table(paste(file.path,file.GM,file, ".",file.Ext.GM,sep=""), head = TRUE,
                  skip=skip.1, nrows=file.fragment) ,silent=TRUE)
                  
        if(inherits(GD, "try-error"))  {
          GD=NULL
          print("File end reached for GD!!!")
        }
        if(inherits(GI, "try-error"))  {
          GI=NULL
          print("File end reached for GI!!!")
        }                          
                  
        if(is.null(GD)) return(list(GD=NULL, GI=NULL,GT=NULL,linesRead=NULL,GLD=NULL))
        
        GT=GD[,1]  #Extract infividual names

        GD=GD[,-1] #Remove individual names
#print("Numerical file read sucesfuly from fragment") 
        linesRead=ncol(GD)       
        if(SNP.fraction==1) return(list(GD=GD, GI=GI,GT=GT,linesRead=linesRead,GLD=NULL))
        
        if(SNP.fraction<1){
          n= ncol(GD)
          #set.seed(seed+file)
          sample=sample(1:n,floor(n*SNP.fraction))
          return(list(GD=GD[,sample], GI=GI[sample,],GT=GT,linesRead=linesRead,GLD=NULL))
        }
    } # end of the "EMMA"
#print("fragment ended succesfully!")
}#End of fragment
#=============================================================================================

