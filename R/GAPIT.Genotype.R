`GAPIT.Genotype` <-
function(G=NULL,GD=NULL,GM=NULL,KI=NULL,
  kinship.algorithm="Zhang",SNP.effect="Add",SNP.impute="Middle",PCA.total=0,PCA.col=NULL,PCA.3d=FALSE,seed=123, SNP.fraction =1,
  file.path=NULL,file.from=NULL, file.to=NULL, file.total=NULL, file.fragment = 1000,SNP.test=TRUE,
  file.G =NULL,file.Ext.G =NULL,
  file.GD=NULL,file.Ext.GD=NULL,
  file.GM=NULL,file.Ext.GM=NULL,
  SNP.MAF=0.05,FDR.Rate = 0.05,SNP.FDR=1,
  Timmer=NULL,Memory=NULL,WS0=1e6,ws=200,Aver.Dis=1000,
  LD.chromosome=NULL,LD.location=NULL,LD.range=NULL, SNP.CV=NULL,
  GP = NULL,GK = NULL,GTindex=NULL,  
  bin.size = 1000,inclosure.size = 100,
  sangwich.top=NULL,sangwich.bottom=NULL,PCA.legend=NULL,
  file.output=TRUE,kinship.cluster="average",NJtree.group=NULL,NJtree.type=c("fan","unrooted"),
  Create.indicator = FALSE, Major.allele.zero = FALSE,Geno.View.output=TRUE){
#Object: To unify genotype and calculate kinship and PC if required:
#       1.For G data, convert it to GD and GI
#       2.For GD and GM data, nothing change 
#       3.Samling GD and create KI and PC
#       4.Go through multiple files
#       5.In any case, GD must be returned (for QC)
#Output: GD, GI, GT, KI and PC
#Authors: Zhiwu Zhang
#Last update: August 11, 2011
##############################################################################################

#print("Genotyping: numericalization, sampling kinship, PCs and much more...")



Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Genotype start")
Memory=GAPIT.Memory(Memory=Memory,Infor="Genotype start")
compress_z=NULL
type_col=NULL
#Create logical variables
byData=!is.null(G) | !is.null(GD)
byFile=!is.null(file.G) | !is.null(file.GD)
hasGenotype=(byData | byFile  )
needKinPC=(is.null(KI) | PCA.total>0 | kinship.algorithm=="Separation")

if(!is.null(KI) & !byData & !byFile & !SNP.test &kinship.algorithm!="SUPER") 
  { 
  print("It return unexpected")
  return (list(GD=NULL,GI=NULL,GT=NULL,hasGenotype=FALSE, genoFormat=NULL, KI=KI,PC=NULL,byFile=FALSE,fullGD=TRUE,Timmer=Timmer,Memory=Memory))
  }


#Set indicator for full GD
fullGD=FALSE
if(byData) fullGD=TRUE
if(byFile & SNP.fraction==1 & needKinPC) fullGD=TRUE

#SET GT to NULL in case of no genotype
if(!byData & !byFile & is.null(GK) &kinship.algorithm!="SUPER") 
  {
  if(is.null(KI) & is.null(GP) & is.null(GK)) stop("GAPIT says: Kinship has to be provided or estimated from genotype!!!")
  return (list(GD=NULL,GI=NULL,GT=NULL,hasGenotype=FALSE, genoFormat=NULL, KI=KI,PC=NULL,byFile=FALSE,fullGD=TRUE,Timmer=Timmer,Memory=Memory))
  }

genoFormat="hapmap"
if(is.null(G)&is.null(file.G)) genoFormat="EMMA"

#Multiple genotype files
#In one of the 3 situations, calculate KI with the algorithm specified, otherwise skip cit by setting algorithm to "SUPER"
kinship.algorithm.save=kinship.algorithm
kinship.algorithm="SUPER"
#Normal
if(is.null(sangwich.top) & is.null(sangwich.bottom) ) kinship.algorithm=kinship.algorithm.save
#TOP or Bottom is MLM
pass.top=FALSE
if(!is.null(sangwich.top))   pass.top=!(sangwich.top=="FaST" | sangwich.top=="SUPER" | sangwich.top=="DC")
pass.bottom=FALSE
if(!is.null(sangwich.bottom))   pass.bottom=!(sangwich.bottom=="FaST" | sangwich.bottom=="SUPER" | sangwich.bottom=="DC")
if(pass.top | pass.bottom )kinship.algorithm=kinship.algorithm.save
#Compatibility of input

#agreement among file from, to and total
if(!is.null(file.from) &!is.null(file.to) &!is.null(file.total))
  {
  if(file.total!=(file.to-file.from+1))  stop("GAPIT says: Conflict among file (from, to and total)")
  }
if(!is.null(file.from) &!is.null(file.to)) 
  {
  if(file.to<file.from)  stop("GAPIT says: file.from should smaller than file.to")
  }
#file.from and file.to must be in pair
if(is.null(file.from) &!is.null(file.to) ) stop("GAPIT says: file.from and file.to must be in pair)")
if(!is.null(file.from) &is.null(file.to) ) stop("GAPIT says: file.from and file.to must be in pair)")

#assign file.total
if(!is.null(file.from) &!is.null(file.to) ) file.total=file.to-file.from+1
if(byFile& is.null(file.total)) stop("GAPIT says: file.from and file.to must be provided!)")

if(!is.null(GP) & !is.null(GK) ) stop("GAPIT Says: You can not provide GP and GK at same time")
if(!is.null(GP) & !is.null(KI) ) stop("GAPIT Says: You can not provide GP and KI at same time")
if(!is.null(GK) & !is.null(KI))   stop("GAPIT says: You can not specify GK and KI at same time!!!")

#GP does not allow TOP
if(!is.null(GP) & !is.null(sangwich.top) ) stop("GAPIT Says: You provided GP. You can not spycify sangwich.top")

#Top require a bottom
if(!is.null(sangwich.top) & is.null(sangwich.bottom) ) stop("GAPIT Says: Top require its Bottom")

#naked bottom require GP or GK
if(is.null(sangwich.top) & !is.null(sangwich.bottom) & (is.null(GP) & is.null(GK)) ) stop("GAPIT Says: Uncovered Bottom (without TOP) requires GP or GK")

#Pseudo top (GK or GP) requires a bottom
if(is.null(sangwich.top) & is.null(sangwich.bottom) & (!is.null(GP)|!is.null(GK  ))) stop("GAPIT Says: You have provide GP or GK, you need to provide Bottom")

#if(!is.null(KI) &!is.null(kinship.algorithm))  stop("GAPIT says: You can not specify kinship.algorithm and provide kinship at same time!!!")



if(!needKinPC &SNP.fraction<1)  stop("GAPIT says: You did not require calculate kinship or PCs. SNP.fraction should not be specified!!!")
if(!SNP.test & is.null(KI) & !byData & !byFile)  stop("GAPIT says: For SNP.test optioin, please input either use KI or use genotype")

#if(is.null(file.path) & !byData & byFile) stop("GAPIT Ssays: A path for genotype data should be provided!")
if(is.null(file.total) & !byData & byFile) stop("GAPIT Says: Number of file should be provided: >=1")
if(!is.null(G) & !is.null(GD)) stop("GAPIT Says: Both hapmap and EMMA format exist, choose one only.")

if(!is.null(file.GD) & is.null(file.GM) & (!is.null(GP)|!is.null(GK)) ) stop("GAPIT Ssays: Genotype data and map files should be in pair")
if(is.null(file.GD) & !is.null(file.GM) & (!is.null(GP)|!is.null(GK)) ) stop("GAPIT Ssays: Genotype data and map files should be in pair")

if(!is.null(GD) & is.null(GM) & (is.null(GP)&is.null(GK)) &kinship.algorithm!="SUPER") stop("GAPIT Says: Genotype data and map files should be in pair")
if(is.null(GD) & !is.null(GM) & (is.null(GP)&is.null(GK)) &kinship.algorithm!="SUPER") stop("GAPIT Says: Genotype data and map files should be in pair")


#if(!byData & !byFile) stop("APIT Ssays: Either genotype data or files should be given!")
#if(byData&(!is.null(file.path))) stop ("APIT Ssays: You have provided geotype data. file.path should not be provided!")

#print("Pass compatibility of input")
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Genotype loaded")
Memory=GAPIT.Memory(Memory=Memory,Infor="Genotype loaded")
  
#Inital GLD
GLD=NULL
SNP.QTN=NULL #Intitial
GT=NULL

#Handler of read data in numeric format (EMMA)
#Rename GM as GI
if(!is.null(GM))GI=GM
rm(GM)
gc()
#Extract GD and GT from read data GD
if(!is.null(GD) )
  {
  GT=as.matrix(GD[,1])  #get taxa
  GD=as.matrix(GD[,-1]) #remove taxa column
  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="GT created from GD)")
  Memory=GAPIT.Memory(Memory=Memory,Infor="GT created from GD")
  }

#Hapmap format
if(!is.null(G))
  {
  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Before HapMap")
  Memory=GAPIT.Memory(Memory=Memory,Infor="Before HapMap")
  #Convert HapMap to numerical
  print(paste("Converting genotype...",sep=""))
  hm=GAPIT.HapMap(G,SNP.effect=SNP.effect,SNP.impute=SNP.impute, Create.indicator = Create.indicator, Major.allele.zero = Major.allele.zero)
  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="after HapMap")
  Memory=GAPIT.Memory(Memory=Memory,Infor="after HapMap")
  #Extracting SNP for LD plot
  if(!is.null(LD.chromosome))
    {
  #print("Extracting SNP for LD plot...")
    chromosome=(G[,3]==LD.chromosome[1])
    bp=as.numeric(as.vector(G[,4]))
    deviation=abs(bp-as.numeric(as.vector(LD.location[1])) )
    location=deviation< as.numeric(as.vector(LD.range[1])  )
    index=chromosome&location
    GLD=G[index,]
    }else{
    #print("No data in GLD")
    GLD=NULL
    }
    Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="HapMap")
    Memory=GAPIT.Memory(Memory=Memory,Infor="HapMap")
    print(paste("Converting genotype done.",sep=""))
    #rm(G)
    #gc()
    Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="G removed")
    Memory=GAPIT.Memory(Memory=Memory,Infor="G removed")
    GT=hm$GT
    GD=hm$GD
    GI=hm$GI
#
#print(unique(GI[,2]))
    rm(hm)
    gc()
    Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="hm removed")
    Memory=GAPIT.Memory(Memory=Memory,Infor="hm removed")
  }

#From files
if(!byData & byFile)
  {
  #print("Loading genotype from files...")
  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="byFile")
  Memory=GAPIT.Memory(Memory=Memory,Infor="byFile")
  numFileUsed=file.to
  if(!needKinPC) numFileUsed=file.from
  #Initial GI as storage
  GD=NULL
  GT=NULL
  GI=NULL
  GLD=NULL
  #multiple fragments or files
  for (file in file.from:numFileUsed)
    {
    frag=1
    numSNP=file.fragment
    myFRG=NULL
   #print(paste("numSNP  before while is ",numSNP))
    while(numSNP==file.fragment) 
         {     #this is problematic if the read end at the last line
         print(paste("Reading file: ",file,"Fragment: ",frag))
         Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Before Fragment")
         Memory=GAPIT.Memory(Memory=Memory,Infor="Before Fragment")
         myFRG=GAPIT.Fragment( file.path=file.path,file.from=file.from, file.to=file.to,file.total=file.total,file.G=file.G,file.Ext.G=file.Ext.G,
                            seed=seed,SNP.fraction=SNP.fraction,SNP.effect=SNP.effect,SNP.impute=SNP.impute,genoFormat=genoFormat,
                            file.GD=file.GD,file.Ext.GD=file.Ext.GD,file.GM=file.GM,file.Ext.GM=file.Ext.GM,
                            file.fragment=file.fragment,file=file,frag=frag,
                            LD.chromosome=LD.chromosome,LD.location=LD.location,LD.range=LD.range, Create.indicator = Create.indicator, Major.allele.zero = Major.allele.zero)
         Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="After Fragment")
         Memory=GAPIT.Memory(Memory=Memory,Infor="After Fragment")
 
         if(is.null(GT) & !is.null(myFRG$GT))GT= as.matrix(myFRG$GT)

         if(is.null(GD))
           {
           GD= myFRG$GD
           }else{
           if(!is.null(myFRG$GD))
             {
             GD=cbind(GD,myFRG$GD)
             }
           }
           if(is.null(GI))
             {
             GI= myFRG$GI
             }else{
             if(!is.null(myFRG$GI)) 
               {
               colnames(myFRG$GI)=c("SNP","Chromosome","Position")
               GI=as.data.frame(rbind(as.matrix(GI),as.matrix(myFRG$GI)))
               }
             }

           if(is.null(G))
             {
             G= myFRG$G
             }else{
             if(!is.null(myFRG$G)) 
               {
               G=as.data.frame(rbind(as.matrix(G),as.matrix(myFRG$G[-1,])))
               }
             }
      
           if(is.null(GLD))
             {
             GLD= myFRG$GLD
             }else{
             if(!is.null(myFRG$GLD))
               {
               if(myFRG$heading)
                 {
                 GLD=as.data.frame(rbind(as.matrix(GLD),as.matrix(myFRG$GLD[-1,])))
                 }else{
                 GLD=as.data.frame(rbind(as.matrix(GLD),as.matrix(myFRG$GLD)))
                 }
               }
             }

            if(file==file.from & frag==1)GT=as.matrix(myFRG$GT)
            frag=frag+1
            if(!is.null(myFRG$GI))
              {
              numSNP=myFRG$linesRead[1]
              }else{
              numSNP=0
              }

            if(!needKinPC)numSNP=0  #force to end the while loop
            if(is.null(myFRG))numSNP=0  #force to end the while loop

            Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="END this Fragment")
            Memory=GAPIT.Memory(Memory=Memory,Infor="END this Fragment")



          } #end whileof repeat on fragment
   # print("This file is OK")
    } #end of file loop
  print("All files loaded")
  } #end of if(!byData&byFile)

#GM=as.matrix(GI)
#GI=GM
# GM=GI

# modified by Jiabo in 20190927. sorted number of chrom by numeric and charicter
chor_taxa=as.character(unique(GI[,2]))
chor_taxa=chor_taxa[order(as.numeric(as.character(chor_taxa)))]
letter.index=grep("[A-Z]|[a-z]",chor_taxa)

if(!setequal(integer(0),letter.index))
  {     
  # myGI=as.matrix(myGI)
      if(length(letter.index)!=length(chor_taxa))
        {
          chr.letter=chor_taxa[letter.index]
          chr.taxa=chor_taxa[-letter.index]
        }else{
          chr.letter=chor_taxa
          chr.taxa=NULL
        }
      Chr=as.character(GI[,2])
      for(i in letter.index)
        {
         index=Chr==chor_taxa[i]
         Chr[index]=i 
        }
      GI[,2]=as.data.frame(Chr)
  }

#print(chor_taxa)
#print(head(GI))
#print("@@@@@@@@@@@")
#print(GD[1:5,1:5])
#print(dim(GI))
#Follow the MAF to filter markers
if(!is.null(GD))
  { 
  #maf=apply(as.matrix(GD),2,function(one) abs(1-sum(one)/(2*nrow(GD))))
  #maf[maf>0.5]=1-maf[maf>0.5]
  ss=apply(GD,2,sum)
  maf=apply(cbind(.5*ss/(nrow(GD)),1-.5*ss/(nrow(GD))),1,min)
#print(max(maf))
#print(min(maf))
  maf_index=maf>=SNP.MAF
  print(paste("GAPIT will filter marker with MAF setting !!"))
  print(paste("The markers will be filtered by SNP.MAF: ",SNP.MAF,sep=""))
  print(table(maf_index))

#print(head(maf[!maf_index]))

  GD=GD[,maf_index]
  GI=as.data.frame(GI[maf_index,])
  # GM=as.data.frame(GM[maf_index,])
  #GI=GM
  }
#print("file loaded")

Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Sampling genotype")
Memory=GAPIT.Memory(Memory=Memory,Infor="Sampling genotype")
#print(KI)
#Plot third part kinship
if(!is.null(KI)&file.output)
  {
  # if(KI!=1) 
    # {
    if(nrow(KI)<2000)
      {
      print("Plotting Kinship")
      #print(dim(KI))
      theKin=as.matrix(KI[,-1])
      line.names <- KI[,1]
      colnames(theKin)=KI[,1]
      rownames(theKin)=KI[,1]
      distance.matrix = stats::dist(theKin,upper=TRUE)
      hc = stats::hclust(distance.matrix,method=kinship.cluster)
      hcd = stats::as.dendrogram(hc)
    ##plot NJtree
      if(!is.null(NJtree.group))
        {
        clusMember <- stats::cutree(hc, k = NJtree.group)
        compress_z=table(clusMember,paste(line.names))
        type_col = grDevices::rainbow(NJtree.group)
        Optimum=c(nrow(theKin),kinship.cluster,NJtree.group)
        }
      Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="set kinship")
      Memory=GAPIT.Memory(Memory=Memory,Infor="set kinship")
      if(file.output)
      {
      print("Creating heat map for kinship...")
      grDevices::pdf(paste("GAPIT.Genotype.Kin_thirdPart.pdf",sep=""), width = 12, height = 12)
      graphics::par(mar = c(25,25,25,25))
      Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="prepare heatmap")
      Memory=GAPIT.Memory(Memory=Memory,Infor="prepare heatmap")
      gplots::heatmap.2(theKin,  cexRow =.2, cexCol = 0.2, col=rev(grDevices::heat.colors(256)), scale="none", symkey=FALSE, trace="none")
      grDevices::dev.off()
      print("Kinship heat map PDF created!") 
      Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="plot heatmap")
      Memory=GAPIT.Memory(Memory=Memory,Infor="plot heatmap")
      }
## Jiabo Wang add NJ Tree of kinship at 4.5.2017
      if(!is.null(NJtree.group)&file.output)
        {            
        for(tr in 1:length(NJtree.type))
           {
           print("Creating NJ Tree for kinship...")
           grDevices::pdf(paste("GAPIT.Genotype.Kin_NJtree_",NJtree.type[tr],".pdf",sep=""), width = 12, height = 12)
           graphics::par(mar = c(5,5,5,5))
           Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="prepare NJ TREE")
           Memory=GAPIT.Memory(Memory=Memory,Infor="prepare NJ TREE")
           plot(ape::as.phylo(hc), type = NJtree.type[tr], tip.color =type_col[clusMember],  use.edge.length = TRUE, col = "gray80",cex=0.8)
           graphics::legend("topright",legend=paste(c("Tatal individuals is: ","Cluster method: ","Group number: "), Optimum[c(1:3)], sep=""),lty=0,cex=1.3,bty="n",bg=graphics::par("bg"))
           grDevices::dev.off()
           }
        }
        if(!is.null(compress_z)){
            utils::write.table(compress_z, 
                paste("GAPIT.Genotype.Kin_NJtre_compress_z.txt",sep=""),
                quote=F)
        }
        print("Kinship NJ TREE PDF created!")
 
        Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="plot NJ TREE")
        Memory=GAPIT.Memory(Memory=Memory,Infor="plot NJ TREE")
    #rm(hc,clusMember)
      }#end 
## NJ Tree end    } #end of if(nrow(KI)<1000)
    # } #end of if(KI!=1)
  } #end of if(!is.null(KI))

Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Before SUPER")
Memory=GAPIT.Memory(Memory=Memory,Infor="Before SUPER")

# SUPER
if(!is.null(GP) & kinship.algorithm=="SUPER" & !is.null(bin.size) & !is.null(inclosure.size))
{
  mySpecify=GAPIT.Specify(GI=GI,GP=GP,bin.size=bin.size,inclosure.size=inclosure.size)
  SNP.QTN=mySpecify$index
  if(!is.null(GD))
  {
	  GK=GD[,SNP.QTN]
    SNPVar=apply(as.matrix(GK),2, stats::var)
    GK=GK[,SNPVar>0]
    GK=cbind(as.data.frame(GT),as.data.frame(GK)) #add taxa  
  } 
}
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Before creating kinship")
Memory=GAPIT.Memory(Memory=Memory,Infor="Before creating kinship")

#Create kinship from genotype if not provide
if(is.null(KI) & (!is.null(GD) |!is.null(GK)) & !kinship.algorithm%in%c("FarmCPU","Blink","MLMM"))
{
  print("Calculating kinship...")
  if(!is.null(GK))
  {
    thisGD=GK[,-1]
    myGT=as.matrix(GK[,1])
    print("GK is used to create KI")
  }else{
    thisGD=GD
    myGT=GT
  }
  print(paste("Number of individuals and SNPs are ",nrow(thisGD)," and ",ncol(thisGD)))
  theKin=NULL
  #if(is.null(PCA.col)&!is.null(NJtree.group))PCA.col=rainbow(NJtree.group)[clusMember]
  if(kinship.algorithm=="EMMA")
    {
    half.thisGD = as.matrix(.5*thisGD)
    if(length(which(is.na(half.thisGD))) > 0)
      {
      print("Substituting missing values with heterozygote for kinship matrrix calculation....")
      half.thisGD[which(is.na(half.thisGD))] = 1
      }
      theKin= emma.kinship(snps=t(as.matrix(.5*thisGD)), method="additive", use="all")
    }
  if(kinship.algorithm=="Loiselle")theKin= GAPIT.kinship.loiselle(snps=t(as.matrix(.5*thisGD)), method="additive", use="all")
  if(kinship.algorithm=="VanRaden")theKin= GAPIT.kinship.VanRaden(snps=as.matrix(thisGD)) 
  if(kinship.algorithm=="Zhang")theKin= GAPIT.kinship.Zhang(snps=as.matrix(thisGD)) 
  if(kinship.algorithm=="Separation")
  {
    thePCA=GAPIT.PCA(X = GD, taxa = GT, PC.number = PCA.total,file.output=F,PCA.total=PCA.total,PCA.col=NULL,PCA.3d=F,PCA.legend=PCA.legend)
    PC=thePCA$PCs[,1:(1+PCA.total)]
    theKin= GAPIT.kinship.separation(PCs=thePCA$PCs,EV=thePCA$EV,nPCs=PCA.total)
  }
  if(!is.null(theKin))
    {
    colnames(theKin)=myGT
    rownames(theKin)=myGT
    line.names <- myGT
    if (!is.null(NJtree.group))
      {
      distance.matrix = stats::dist(theKin,upper=TRUE)
      hc = stats::hclust(distance.matrix, method = kinship.cluster)
      hcd = stats::as.dendrogram(hc)
      clusMember <- stats::cutree(hc, k = NJtree.group)
      compress_z=table(clusMember,paste(line.names))
      type_col = grDevices::rainbow(NJtree.group)
      Optimum=c(nrow(theKin),kinship.cluster,NJtree.group)
      }
    print("kinship calculated")
    if(length(GT)<2000 &file.output)
      {
    #Create heat map for kinship
      print("Creating heat map for kinship...")
      grDevices::pdf(paste("GAPIT.Genotype.Kin_",kinship.algorithm,".pdf",sep=""), width = 12, height = 12)
      graphics::par(mar = c(25,25,25,25))
      gplots::heatmap.2(theKin,  cexRow =.2, cexCol = 0.2, col=rev(grDevices::heat.colors(256)), scale="none", symkey=FALSE, trace="none")
      grDevices::dev.off()
      print("Kinship heat map created")
    ## Jiabo Wang add NJ Tree of kinship at 4.5.2017
      if (!is.null(NJtree.group))      
        {
        print("Creating NJ Tree for kinship...")
        for(tr in 1:length(NJtree.type))
           {
           grDevices::pdf(paste("GAPIT.Genotype.Kin_NJtree_",NJtree.type[tr],".pdf",sep=""), width = 12, height = 12)
           graphics::par(mar = c(0,0,0,0))
           Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="prepare NJ TREE")
           Memory=GAPIT.Memory(Memory=Memory,Infor="prepare NJ TREE")   
           plot(ape::as.phylo(hc), type = NJtree.type[tr], tip.color =type_col[clusMember],  use.edge.length = TRUE, col = "gray80",cex=0.6)
    #legend("topright",legend=c(paste("Tatal numerber of individuals is ",),lty=0,cex=1.3,bty="n",bg=par("bg"))
           graphics::legend("topright",legend=paste(c("Tatal individuals is: ","Group method: ","Group number: "), Optimum[c(1:3)], sep=""),lty=0,cex=1.3,bty="n",bg=graphics::par("bg"))
           grDevices::dev.off()
           }
    # print(Optimum)   
        utils::write.table(compress_z,paste("GAPIT.Genotype.Kin_NJtree_compress_z.txt",sep=""),quote=F)
        print("Kinship NJ TREE PDF created!")  
        Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="plot NJ TREE")
        Memory=GAPIT.Memory(Memory=Memory,Infor="plot NJ TREE")
        rm(hc)
        }#end NJtree
      }
    print("Adding IDs to kinship...")
    #Write the kinship into a text file
    KI=cbind(myGT,as.data.frame(theKin)) #This require big memory. Need a way to solve it.
    print("Writing kinship to file...")
    if(file.output) utils::write.table(KI, paste("GAPIT.Genotype.Kin_",kinship.algorithm,".csv",sep=""), quote = FALSE, sep = ",", row.names = FALSE,col.names = FALSE)
    print("Kinship save as file")    
    rm(theKin)
    gc()
    }
  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Estimating kinship")
  Memory=GAPIT.Memory(Memory=Memory,Infor="Estimating kinship")
  print("Kinship created!")
}  #end of if(is.null(KI)&!is.null(GD))

Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="after creating kinship")
Memory=GAPIT.Memory(Memory=Memory,Infor="after creating kinship")


PC=NULL
thePCA=NULL

if(PCA.total>0)
{
  if(is.null(PCA.col)&!is.null(type_col))PCA.col=type_col[clusMember]
  thePCA=GAPIT.PCA(X = GD, taxa = GT, PC.number = PCA.total,file.output=file.output,PCA.total=PCA.total,PCA.col=PCA.col,PCA.3d=PCA.3d)
  PC=thePCA$PCs[,1:(1+PCA.total)]
  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="PCA")
  Memory=GAPIT.Memory(Memory=Memory,Infor="PCA")
  print("PC created")
}

#LD plot
#print("LD section")
if(!is.null(GLD) &file.output)
  {
    # print(dim(GLD))
  if(nrow(GLD)>500)
    {
    GLD=GLD[1,]
    print("WARNING: The number of SNPs requested is beyond limitation. No LD plot created.")
    }
  if(nrow(GLD)>1)
    {
    print("Plot LD...")
    hapmapgeno= data.frame(as.matrix(t(GLD[,-c(1:11)])))
    hapmapgeno[hapmapgeno=="NN"]=NA
    hapmapgeno[hapmapgeno=="XX"]=NA
    hapmapgeno[hapmapgeno=="--"]=NA
    hapmapgeno[hapmapgeno=="++"]=NA
    hapmapgeno[hapmapgeno=="//"]=NA
    # print(GLD)
    LDdist=as.numeric(as.vector(GLD[,4]))
    LDsnpName=GLD[,1]
    colnames(hapmapgeno)=LDsnpName
    sigsnp=match(LD.location,LDdist)
#Prune SNM names
#LDsnpName=LDsnpName[GAPIT.Pruning(LDdist,DPP=7)]
    LDsnpName=LDsnpName[c(sigsnp)] #keep the first and last snp names only
    # print(LDsnpName)
    color.rgb <- grDevices::colorRampPalette(rev(c("snow","red")),space="rgb")
    
    # print(color.rgb)
    print("Getting genotype object")
    LDsnp=genetics::makeGenotypes(hapmapgeno,sep="",method=genetics::as.genotype)   #This need to be converted to genotype object
    print("Caling LDheatmap...")
#pdf(paste("GAPIT.LD.pdf",sep=""), width = 12, height = 12)
    grDevices::pdf(paste("GAPIT.Genotype.LD_chromosom_",LD.chromosome,"(",round(max(0,LD.location-LD.range)/1000000),"_",round((LD.location+LD.range)/1000000),"Mb)",".pdf",sep=""), width = 12, height = 12)
    graphics::par(mar = c(25,25,25,25))

    MyHeatmap <- try(LDheatmap::LDheatmap(LDsnp, LDdist, LDmeasure="r", add.map=TRUE,
    SNP.name=LDsnpName,color=color.rgb(20), name="myLDgrob", add.key=TRUE,geneMapLabelY=0.1) )  
    if(!inherits(MyHeatmap, "try-error")) 
      {
  #Modify the plot
#      library(grid)
      # LDheatmap.highlight(MyHeatmap, i = sigsnp, j=sigsnp+1, col = "red")
      grid::grid.edit(grid::gPath("myLDgrob","heatMap","heatmap"),gp=grid::gpar(col="white",lwd=8))
      grid::grid.edit(grid::gPath("myLDgrob", "Key", "title"), gp=grid::gpar(cex=.5, col="blue"))  #edit key title size and color
      grid::grid.edit(grid::gPath("myLDgrob", "heatMap", "title"), gp=grid::gpar(just=c("center","bottom"), cex=0.8, col="black")) #Edit gene map title
      grid::grid.edit(grid::gPath("myLDgrob", "geneMap","SNPnames"), gp = grid::gpar(cex=0.3,col="black")) #Edit SNP name
      }else{
      print("Warning: error in converting genotype. No LD plot!")
      }
    grDevices::dev.off()
    print("LD heatmap crated")
    }else{ # alternative of if(nrow(GLD)>1)
    print("Warning: There are less than two SNPs on the region you sepcified. No LD plot!")
    } #end of #if(nrow(GLD)>1)
  }#end of if(!is.null(GLD))

Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="after LD plot")
Memory=GAPIT.Memory(Memory=Memory,Infor="after LD plot")


###output Marker density and decade of linkage disequilibrium over distance
if(!is.null(GI) & !is.null(GD) & file.output & Geno.View.output)
{
ViewGenotype<-GAPIT.Genotype.View(
GI=GI,
X=GD,
WS0=WS0,
ws=ws,
Aver.Dis=Aver.Dis
)
}

#print("Genotype successfully acomplished")
return (list(G=G,GD=GD,GI=GI,GT=GT,hasGenotype=hasGenotype, genoFormat=genoFormat, KI=KI,PC=PC,byFile=byFile,fullGD=fullGD,Timmer=Timmer,Memory=Memory,SNP.QTN=SNP.QTN,chor_taxa=chor_taxa))
}
#=============================================================================================

