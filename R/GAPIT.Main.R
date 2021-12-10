#'
#' GAPIT.Main
#' 
#' @description 
#' GAPIT.Main
#'
#' @param Y data.frame of phenotype data, samples in rows, traits in columns
#' @param G data.frame of genotypic data, HAPMAP format
#' @param GD data.frame of genotypic data
#' @param GM genetic map for GD
#' @param KI param
#' @param Z param
#' @param CV param
#' @param CV.Inheritance param
#' @param SNP.P3D param
#' @param GP param
#' @param GK param
#' @param group.from param
#' @param group.to param
#' @param group.by param
#' @param kinship.cluster param
#' @param kinship.group param
#'
#' @param kinship.algorithm param
#' @param DPP param
#' @param ngrid param
#' @param llin param
#' @param ulim param
#' @param esp param
#' @param GAPIT3.output param
#' @param file.path param
#' @param file.from param
#' @param file.to param
#' @param file.total param
#' @param file.fragment param
#' @param file.G param
#' @param file.Ext.G param
#' @param file.GD param
#' @param file.GM param
#' @param file.Ext.GD param
#' @param file.Ext.GM param
#' @param SNP.MAF param
#' @param FDR.Rate param
#' @param SNP.FDR param
#' @param SNP.effect param
#' @param SNP.impute param
#' @param PCA.total param
#' @param GAPIT.Version param
#' @param name.of.trait param
#' @param GT param
#' @param SNP.fraction param
#' @param seed param
#' @param BINS param
#' @param SNP.test param
#' @param SNP.robust param
#' @param LD.chromosome param
#' @param LD.location param
#' @param LD.range param
#' @param model param
#' @param bin.from param
#' @param bin.to param
#' @param bin.by param
#' @param inclosure.from param
#' @param inclosure.to param
#' @param inclosure.by param
#' @param SNP.permutation param
#' @param SNP.CV param
#' @param NJtree.group param
#' @param NJtree.type param
#' @param plot.bin param
#' @param genoFormat param
#' @param hasGenotype param
#' @param byFile param
#' @param fullGD param
#' @param PC param
#' @param GI param
#' @param Timmer param
#' @param Memory param
#' @param sangwich.top param
#' @param sangwich.bottom param
#' @param QC param
#' @param GTindex param
#' @param LD param
#' @param file.output param
#' @param cutOff param
#' @param Model.selection param
#' @param Create.indicator param
#' @param QTN param
#' @param QTN.round param
#' @param QTN.limit param
#' @param QTN.update param
#' @param QTN.method param
#' @param Major.allele.zero param
#' @param QTN.position param
#' @param SUPER_GD param
#' @param SUPER_GS param
#' @param plot.style param
#' @param CG param
#' @param chor_taxa param
#'
#' @return 
#' A list
#'
#' @author Zhiwu Zhang and Jiabo Wang
#'
#' @export
`GAPIT.Main` <-
function(Y,
         G=NULL,
         GD=NULL,
         GM=NULL,
         KI=NULL,
         Z=NULL,
         CV=NULL,
         CV.Inheritance=NULL,
         SNP.P3D=TRUE,
         GP=NULL,
         GK=NULL,
         group.from=1000000,
         group.to=1,
         group.by=10,
         kinship.cluster="average",
         kinship.group='Mean',
         kinship.algorithm=NULL,
         DPP=50000,
         ngrid = 100, 
         llin = -10, 
         ulim = 10, 
         esp = 1e-10,
         GAPIT3.output=TRUE,
         file.path=NULL,
         file.from=NULL, 
         file.to=NULL, 
         file.total=NULL, 
         file.fragment = 512, 
         file.G=NULL, 
         file.Ext.G=NULL,
         file.GD=NULL, 
         file.GM=NULL, 
         file.Ext.GD=NULL,
         file.Ext.GM=NULL,
         SNP.MAF=0,
         FDR.Rate=1,
         SNP.FDR=1,
         SNP.effect="Add",
         SNP.impute="Middle",
         PCA.total=0,  
         GAPIT.Version=GAPIT.Version,
         name.of.trait, 
         GT = NULL, 
         SNP.fraction = 1, 
         seed = 123, 
         BINS = 20,
         SNP.test=TRUE,
         SNP.robust="FaST",
         LD.chromosome=NULL,
         LD.location=NULL,
         LD.range=NULL,
         model=model,
         bin.from=10000,
         bin.to=5000000,
         bin.by=1000,
         inclosure.from=10,
         inclosure.to=1000,
         inclosure.by=10,
         SNP.permutation=FALSE,
         SNP.CV=NULL,
         NJtree.group=NJtree.group,
         NJtree.type=NJtree.type,
         plot.bin=plot.bin,
         genoFormat=NULL,
         hasGenotype=NULL,
         byFile=NULL,
         fullGD=NULL,
         PC=NULL,
         GI=NULL, 
         Timmer = NULL, 
         Memory = NULL,
         sangwich.top=NULL,
         sangwich.bottom=NULL,
         QC=TRUE,
         GTindex=NULL,
         LD=0.05,
         file.output=TRUE,
         cutOff=0.05, 
         Model.selection = FALSE, 
         Create.indicator = FALSE,
				 QTN=NULL, 
				 QTN.round=1,
				 QTN.limit=0, 
				 QTN.update=TRUE, 
				 QTN.method="Penalty", 
				 Major.allele.zero = FALSE,
         QTN.position=NULL,
				 SUPER_GD=NULL,
				 SUPER_GS=SUPER_GS,
				 plot.style="Beach",
				 CG=CG,
				 chor_taxa=chor_taxa){
#Object: To perform GWAS and GPS (Genomic Prediction or Selection)
#Output: GWAS table (text file), QQ plot (PDF), Manhattan plot (PDF), genomic prediction (text file), and
#        genetic and residual variance components
#Authors: Zhiwu Zhang
# Last update: Oct 23, 2015  by Jiabo Wang add REML threshold and SUPER GD KI
##############################################################################################

#Initial p3d and h2.opt temporaryly
  h2.opt=NULL
  p3d=list(
    ps=NULL,
    REMLs=NULL,
    stats=NULL,
    effect.est=NULL,
    rsquare_base=NULL,
    rsquare=NULL,
    dfs=NULL,
    df=NULL,
    tvalue=NULL,
    stderr=NULL,
    maf=NULL,
    nobs=NULL,
    Timmer=NULL,
    Memory=NULL,
    vgs=NULL,
    ves=NULL,
    BLUP=NULL,
    BLUP_Plus_Mean=NULL,
    PEV=NULL,
    BLUE=NULL,
    logLM=NULL,
    effect.snp=NULL,
    effect.cv=NULL
  )
  
  
  if (SUPER_GS){
    Compression=NULL
    kinship.optimum=NULL
    kinship=NULL
    PC=PC
    REMLs=NULL
    GWAS=NULL
    QTN=NULL
    Timmer=GAPIT.Timmer(Infor="GAPIT.SUPER.GS")
    Memory=GAPIT.Memory(Infor="GAPIT.SUPER.GS")
    #print(model)
    SUPER_GS_GAPIT = GAPIT.SUPER.GS(Y=Y,
                                    GD=GD,
                                    GM=GM,
                                    KI=KI,
                                    Z=Z,
                                    CV=CV,
                                    GK=GK,
                                    kinship.algorithm=kinship.algorithm,
                                    bin.from=bin.from,
                                    bin.to=bin.to,
                                    bin.by=bin.by,
                                    inclosure.from=inclosure.from,
                                    inclosure.to=inclosure.to,
                                    inclosure.by=inclosure.by,
                                    group.from=group.from,
                                    group.to=group.to,
                                    group.by=group.by,
                                    kinship.cluster=kinship.cluster,
                                    kinship.group=kinship.group,
                                    PCA.total=PCA.total,
                                    GT=GT,
                                    PC=PC,
                                    GI=GI,
                                    Timmer = Timmer, 
                                    Memory = Memory,
                                    model=model,
                                    sangwich.top=sangwich.top,
                                    sangwich.bottom=sangwich.bottom,
                                    QC=QC,
                                    GTindex=GTindex,
                                    LD=LD,
                                    file.output=GAPIT3.output,
                                    cutOff=cutOff
                        )
# Compression=as.matrix(SUPER_GS_GAPIT$Compression)
# opt=
	  print("SUPER_GS_GAPIT FUNCTION DONE")	
	  return (list(Compression=SUPER_GS_GAPIT$Compression,
	               kinship.optimum=SUPER_GS_GAPIT$SUPER_kinship,
	               kinship=SUPER_GS_GAPIT$kinship, 
	               PC=SUPER_GS_GAPIT$PC,
	               GWAS=GWAS, 
                 GPS=SUPER_GS_GAPIT$GPS,
	               Pred=SUPER_GS_GAPIT$Pred,
	               Timmer=Timmer,
	               Memory=Memory,
	               h2=SUPER_GS_GAPIT$h2,
	               SUPER_GD=SUPER_GS_GAPIT$SUPER_GD,
	               GWAS=NULL,
	               QTN=NULL)
	          )
					
  }else{
  #print("@@@@@@@")
  #print(group.from)

#Handler of SNP.test=F
#Iniciate with two by seven NA matrix
#The seventh is for p values of SNP
    DTS=rbind(rep(NA,7),rep(NA,7) )
  
  
#End imediatly in one of these situtiona
    shortcut=FALSE
    LL.save=1e10
#In case of null Y and null GP, sent back genotype only  
    thisY=Y[,2]
    thisY=thisY[!is.na(thisY)]
    if(length(thisY) <3){
      shortcut=TRUE
    }else{
      if(stats::var(thisY) ==0) shortcut=TRUE
    }
        
    if(shortcut){
      print(paste("Y is empty. No GWAS/GS performed for ",name.of.trait,sep=""))
      return (list(compression=NULL,
                   kinship.optimum=NULL, 
                   kinship=KI,
                   PC=PC,
                   GWAS=NULL, 
                   GPS=NULL,
                   Pred=NULL, 
                   REMLs=NULL,
                   Timmer=Timmer,
                   Memory=Memory,
                   h2=NULL))
    }

#QC
    print("------------Examining data (QC)------------------------------------------")
    if(is.null(Y)) stop ("GAPIT says: Phenotypes must exist.")
    if(is.null(KI)&missing(GD) & kinship.algorithm!="SUPER") stop ("GAPIT says: Kinship is required. As genotype is not provided, kinship can not be created.")

#When GT and GD are missing, force to have fake ones (creating them from Y),GI is not required in this case
    if(is.null(GD) & is.null(GT)) {
	    GT=as.matrix(Y[,1])
	    GD=matrix(1,nrow(Y),1)	
      GI=as.data.frame(matrix(0,1,3) )
      colnames(GI)=c("SNP","Chromosome","Position")
    }

    if(is.null(GT)) {
      GT=as.character(GD[,1])
    }
#print("@@@@@@@@")
#print(GD)
#merge CV with PC: Put CV infront of PC
    if(PCA.total>0&!is.null(CV))CV=GAPIT.CVMergePC(CV,PC)
    if(PCA.total>0&is.null(CV))CV=PC
    #for GS merge CV with GD name
    if (is.null(CV)){
      my_allCV=CV
    }else{
      taxa_GD=rownames(GD)
      my_allCV=CV[order(CV[,1]),]
      my_allCV=my_allCV[my_allCV[,1]%in%taxa_GD,]
    #print(dim(my_allCV))
    }

    #Handler of CV.Inheritance
    if(is.null(CV) & !is.null(CV.Inheritance)){
      stop ("GAPIT says: CV.Inheritance is more than avaiable.")
    }

    if(!is.null(CV)& !is.null(CV.Inheritance)){  
      if(CV.Inheritance>(ncol(CV)-1)){
        stop ("GAPIT says: CV.Inheritance is more than avaiable.")
      }
    }

    #Create Z as identity matrix from Y if it is not provided
    if(kinship.algorithm!="None" & kinship.algorithm!="SUPER" & is.null(Z)){
      taxa=as.character(Y[,1]) #this part will make GS without CV not present all prediction
      Z=as.data.frame(diag(1,nrow(Y)))
#taxa=as.character(KI[,1])
#Z=as.data.frame(diag(1,nrow(KI)))
      Z=rbind(taxa,Z)
      taxa=c('Taxa',as.character(taxa))
      Z=cbind(taxa,Z)
    }
    ZI=Z

    #Add the part of non proportion in Z matrix
    if(kinship.algorithm!="None" & kinship.algorithm!="SUPER" & !is.null(Z))
    {
      if(nrow(Z)-1<nrow(Y)) Z=GAPIT.ZmatrixFormation(Z=Z,Y=Y)
    }

    #Create CV with all 1's if it is not provided
    noCV=FALSE
    if(is.null(CV)){
      noCV=TRUE
      CV=Y[,1:2]
      CV[,2]=1
      colnames(CV)=c("taxa","overall")
    }

    #Remove duplicat and integragation of data
    print("QC is in process...")

    CVI <- CV

# print(dim(Z))
# print("!!!!!")
# if(QC)
# {
#   qc <- GAPIT.QC(Y=Y,KI=KI, GT=GT,CV=CV,Z=Z,GK=GK)
#   GTindex=qc$GTindex
#   Y=qc$Y  # here make twice qc and chaos with numeric taxa, Thanks for Dennis. 20210913
#   KI=qc$KI
#   CV=qc$CV
#   Z=qc$Z
#   GK=qc$GK
#   # if(noCV)CVI=qc$CV #this part will make GS without CV not present all prediction
# }

# print(length(GK))
    GTindex=match(as.character(KI[,1]),as.character(Y[,1]))
    GTindex=GTindex[!is.na(GTindex)]
# KI=KI[GTindex,GTindex+1]
    my_taxa=as.character(KI[,1])
    CV=CV[as.character(CV[,1])%in%as.character(Y[,1]),]
# print(dim(KI))
# print(dim(CV))

    #Output phenotype
    colnames(Y)=c("Taxa",name.of.trait)
    if(file.output){
      try(utils::write.table(Y, paste("GAPIT.", name.of.trait,".phenotype.csv" ,sep = ""),
                             quote = FALSE, sep = ",", 
                             row.names = FALSE,
                             col.names = TRUE))
    }
    
    # Default kinship.algorithm = "VanRaden".
    # This if() may be seldom used.
    #TDP
    if( kinship.algorithm =="None" ){
      if(min(CV[,2])==max(CV[,2])) CV=NULL
	
      # GAPIT.TDP() does not appear to exist.
      # theTDP = GAPIT.TDP(Y=Y,
      #                    CV=CV,
      #                    SNP = as.data.frame(cbind(GT[GTindex],as.matrix(as.data.frame(GD[GTindex,])))),
      #                    QTN=QTN,
      #                    Round=QTN.round,
      #                    QTN.limit=QTN.limit, 
      #                    QTN.update=QTN.update, 
      #                    Method=QTN.method
      #                    )
#print(dim(GM))
#print(length(theTDP$p))

      #theGWAS=cbind(GM,theTDP$p,NA,NA,NA)	
      theGWAS=cbind(GM,NA,NA,NA,NA)	
      
      return (list(Compression = NULL,
                   kinship.optimum = NULL, 
                   kinship = NULL,
                   PC = NULL,
                   GWAS = theGWAS, 
                   GPS = NULL,
                   Pred = NULL,
                   REMLs = NULL,
                   #QTN = theTDP$QTN,
                   QTN = NULL,
                   Timmer = Timmer,
                   Memory = Memory,
                   h2 = NULL))
    }

#rm(qc)
    gc()

    Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="QC")
    Memory=GAPIT.Memory(Memory=Memory,Infor="QC")

    #Get indicator of sangwich top and bottom
    byPass.top=FALSE
    byPass=FALSE
    NOBLUP=FALSE
    if(group.from<2&group.to<2) NOBLUP=TRUE
    #if(!is.null(sangwich.bottom)) byPass=((sangwich.bottom=="FaST" | sangwich.bottom=="SUPER" | sangwich.bottom=="DC" )& is.null(GP)   )
    if(!is.null(sangwich.top)){
      byPass.top=((sangwich.top=="FaST" | sangwich.top=="SUPER" | sangwich.top=="DC" ) )
    }
    if(!is.null(sangwich.bottom)){
      byPass=((sangwich.bottom=="FaST" | sangwich.bottom=="SUPER" | sangwich.bottom=="DC" ) )
    }
    
    print("Try to group from and to were set to 1")

    if(byPass){
      print("group from and to were set to 1")
      group.from=1
      group.to=1
    }

    print("------------Examining data (QC) done-------------------------------------")

    #Sagnwich top bun: To gep GP if it is not provided
    if(!is.null(sangwich.top) & is.null(GP)){
      print("-------------------Sandwich top bun-----------------------------------")
#print(dim(GD))
#print(GD[1:5,1:5])

      #Create GK if not provided
      if(is.null(GK)){
#    set.seed(1)
        nY=floor(nrow(Y)*.9)
        nG=ncol(GD)
        if(nG>nY){
          snpsam=sample(1:nG,nY)
        }else{
            snpsam=1:nG
        }
        GK=GD[GTindex,snpsam]
    # print(dim(GK))
    # print(GK[270:279,1:5])
        SNPVar=apply(as.matrix(GK), 2, stats::var)
    # print(SNPVar)
        GK=GK[,SNPVar>0]
        GK=cbind(as.data.frame(GT[GTindex]),as.data.frame(GK)) #add taxa
      }
  
  #myGD=cbind(as.data.frame(GT),as.data.frame(GD)) 
      file.output.temp=file.output
      file.output=FALSE
  #print(sangwich.top)[GTindex,c(1,GTindex+1)]
      GP=GAPIT.Bread(Y=Y,CV=CV,Z=Z,KI=KI,GK=GK,GD=cbind(as.data.frame(GT[GTindex]),as.data.frame(GD[GTindex,])),GM=GI,method=sangwich.top,GTindex=GTindex,LD=LD,file.output=file.output)$GWAS
      file.output=file.output.temp
  
  
      GK=NULL
  
      print("-------------------Sagnwich top bun: done-----------------------------")  

    } 

    Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="SagnwichTop")
    Memory=GAPIT.Memory(Memory=Memory,Infor="SagnwichTop")

    #Sandwich burger and dressing
    print("-------------------Sandwich burger and dressing------------------------")

    #Handler of group boundry
    if(group.from>group.to) stop("GAPIT says: group.to should  be larger than group.from. Please correct them!")

    if(is.null(CV) | (!is.null(CV) & group.to<(ncol(CV)+1))) {
      #The minimum of group is 1 + number of columns in CV
      group.from=1
      group.to=1
      warning("The upper bound of groups (group.to) is not sufficient. both boundries were set to a and GLM is performed!")
    }

    if(!is.null(CV)& group.from<1) {
      group.from=1 #minimum of group is number of columns in CV
      warning("The lower bound of groups should be 1 at least. It was set to 1!")
    }
 
    nk=1000000000
    if(!is.null(KI)) nk=min(nk,nrow(KI))
    if(!is.null(GK)) nk=min(nk,nrow(GK))

    if(!is.null(KI)){
      if(group.to>nk) {
        #group.to=min(nrow(KI),length(GTindex)) #maximum of group is number of rows in KI
        group.to=nk #maximum of group is number of rows in KI
        warning("The upper bound of groups is too high. It was set to the size of kinship!") 
      }
      if(group.from>nk){ 
        group.from=nk
        warning("The lower bound of groups is too high. It was set to the size of kinship!") 
      } 
    }

    if(!is.null(CV)){
      if(group.to<=ncol(CV)+1) {
        #The minimum of group is number of columns in CV
        #group.from=ncol(CV)+2
        #group.to=ncol(CV)+2
        warning("The upper bound of groups (group.to) is not sufficient. both boundries were set to their minimum and GLM is performed!")
      }
    }

#bin.fold=ceiling(log2(bin.to/bin.from))
#bin.seq=0:bin.fold
#bin.level=bin.from*2^bin.seq

    #Set upper bound for inclosure.to
    if(inclosure.to>nrow(Y))inclosure.to=nrow(Y)-1

    #set inclosure loop levels
    bin.level=seq(bin.from,bin.to,by=bin.by)
    inclosure=seq(inclosure.from,inclosure.to,by=inclosure.by)

    #Optimization for group number, cluster algorithm and kinship type
    GROUP=seq(group.to,group.from,by=-group.by)#The reverse order is to make sure to include full model
    if(missing("kinship.cluster")) kinship.cluster=c("ward", "single", "complete", "average", "mcquitty", "median", "centroid")
    if(missing("kinship.group")) kinship.group=c("Mean", "Max", "Min", "Median")
    numSetting=length(GROUP)*length(kinship.cluster)*length(kinship.group)*length(bin.level)*length(inclosure)

    #Reform Y, GD and CV into EMMA format
    ys=as.matrix(Y[,2])
    X0=as.matrix(CV[,-1])
    CV.taxa=CVI[,1]
    #print(length(ys))
    #Initial
    count=0
    Compression=matrix(,numSetting,6)
    colnames(Compression)=c("Type","Cluster","Group","REML","VA","VE")

    #add indicator of overall mean
    if(min(X0[,1])!=max(X0[,1])) X0 <- cbind(1, X0) #do not add overall mean if X0 has it already at first column


    Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="DataProcessing")
    Memory=GAPIT.Memory(Memory=Memory,Infor="DataProcessing")

    print("-------------------------Iteration in process--------------------------")
    print(paste("Total iterations: ",numSetting,sep=""))

    #Loop to optimize cluster algorithm, group number and kinship type
    for (bin in bin.level){
      for (inc in inclosure){

        #Grill: update KI if GK or GP is provided
        if(!byPass & (!is.null(GK) | !is.null(GP))){  
          print("Grilling KI...")

          myGenotype<-GAPIT.Genotype(G=NULL,
                                     GD=cbind(as.data.frame(GT),as.data.frame(GD)),
                                     GM=GI,
                                     KI=NULL,
                                     kinship.algorithm=kinship.algorithm,
                                     PCA.total=0,
                                     SNP.fraction=SNP.fraction,
                                     SNP.test=SNP.test,
                                     file.path=file.path,
                                     file.from=file.from, 
                                     file.to=file.to,
                                     file.total=file.total, 
                                     file.fragment = file.fragment, 
                                     file.G=file.G,
                                     file.Ext.G=file.Ext.G,
                                     file.GD=file.GD, 
                                     file.GM=file.GM, 
                                     file.Ext.GD=file.Ext.GD,
                                     file.Ext.GM=file.Ext.GM,
                                     SNP.MAF=SNP.MAF,
                                     FDR.Rate = FDR.Rate,
                                     SNP.FDR=SNP.FDR,
                                     SNP.effect=SNP.effect,
                                     SNP.impute=SNP.impute,
                                     kinship.cluster=kinship.cluster,
                                     NJtree.group=NJtree.group,
                                     NJtree.type=NJtree.type,
                                     LD.chromosome=LD.chromosome,
                                     LD.location=LD.location,
                                     LD.range=LD.range,
                                     GP=GP,
                                     GK=GK,
                                     bin.size=bin,
                                     inclosure.size=inc,
                                     SNP.CV=SNP.CV,
                                     Timmer = Timmer, 
                                     Memory = Memory,
                                     GTindex=GTindex,
                                     sangwich.top=NULL,
                                     sangwich.bottom=sangwich.bottom,
                                     file.output=file.output,
                                     Create.indicator = Create.indicator,
                                     Major.allele.zero = Major.allele.zero)
   
          Timmer=myGenotype$Timmer
          Memory=myGenotype$Memory

          KI=myGenotype$KI
          #update group set by new KI
          nk=nrow(KI)
          GROUP=GROUP[GROUP<=nk]
        }
        
        
        for (ca in kinship.cluster){
          for (group in GROUP){
            for (kt in kinship.group){
              #Do not screen SNP unless existing genotype and one combination
              if(numSetting==1 & hasGenotype){
                optOnly=FALSE
              }else{
                optOnly=TRUE
              }
              if(!SNP.test) optOnly=TRUE

              if(optOnly | Model.selection){
                colInclude=1
                optOnly = TRUE
              }else{
                colInclude=c(1:ncol(GD))
              }

              if(!optOnly) {print("Compressing and Genome screening..." )}
              count=count+1

#Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="PreP3D 1")
#Memory=GAPIT.Memory(Memory=Memory,Infor="PreP3D 1")

              if(!byPass){
                if(count==1)print("-------Mixed model with Kinship-----------------------------")
                if(group<ncol(X0)+1) group=1 # the emma function (emma.delta.REML.dLL.w.Z) does not allow K has dim less then CV. turn to GLM (group=1)

                cp <- GAPIT.Compress(KI=KI,kinship.cluster=ca,kinship.group=kt,GN=group,Timmer=Timmer,Memory=Memory)
                Timmer=cp$Timmer
                Memory=cp$Memory

                Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="PreP3D 2_cp")
                Memory=GAPIT.Memory(Memory=Memory,Infor="PreP3D 2_cp")

#print("BK...")

                bk <- GAPIT.Block(Z=Z,GA=cp$GA,KG=cp$KG)

                Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="PreP3D 2_bk")
                Memory=GAPIT.Memory(Memory=Memory,Infor="PreP3D 2 bk")

#print("ZC...")
                zc <- GAPIT.ZmatrixCompress(Z=Z,GAU =bk$GA)

                Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="PreP3D 2_zc")
                Memory=GAPIT.Memory(Memory=Memory,Infor="PreP3D 2 zc")

#print("wraping...")
#Reform KW and Z into EMMA format

                zrow=nrow(zc$Z)
                zcol=ncol(zc$Z)-1
#Z1=matrix(as.numeric(as.matrix(zc$Z[,-1])),nrow=zrow,ncol=zcol)

                Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Prio PreP3D")
                Memory=GAPIT.Memory(Memory=Memory,Infor="Prio PreP3D")

                p3d <- GAPIT.EMMAxP3D(ys=ys,
                                      xs=as.matrix(as.data.frame(GD[GTindex,colInclude])),
                                      K = as.matrix(bk$KW),
                                      Z=matrix(as.numeric(as.matrix(zc$Z[,-1])),nrow=zrow,ncol=zcol),
                                      X0=X0,
                                      CVI=CVI,
                                      CV.Inheritance=CV.Inheritance,
                                      GI=GI,
                                      SNP.P3D=SNP.P3D,
                                      Timmer=Timmer,
                                      Memory=Memory,
                                      fullGD=fullGD,
                                      SNP.permutation=SNP.permutation,
                                      GP=GP,
                                      SNP.fraction=SNP.fraction,
			                                file.path=file.path,
			                                file.from=file.from,
			                                file.to=file.to,
			                                file.total=file.total, 
			                                file.fragment = file.fragment, 
			                                byFile=byFile, 
			                                file.G=file.G,
			                                file.Ext.G=file.Ext.G,
			                                file.GD=file.GD, 
			                                file.GM=file.GM, 
			                                file.Ext.GD=file.Ext.GD,
			                                file.Ext.GM=file.Ext.GM,
			                                GTindex=GTindex,
			                                genoFormat=genoFormat,
			                                optOnly=optOnly,
			                                SNP.effect=SNP.effect,
			                                SNP.impute=SNP.impute,
			                                name.of.trait=name.of.trait, 
			                                Create.indicator = Create.indicator, 
			                                Major.allele.zero = Major.allele.zero
			                                )

                Timmer=p3d$Timmer
                Memory=p3d$Memory

                Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Post PreP3D")
                Memory=GAPIT.Memory(Memory=Memory,Infor="Post PreP3D")

#print("Cluster algorithm, kinship type, groups, VG, Ve and REML:")
                print(paste(count, "of",numSetting,"--","Vg=",round(p3d$vgs,4), "VE=",round(p3d$ves,4),"-2LL=",round(p3d$REMLs,2), "  Clustering=",ca,"  Group number=", group ,"  Group kinship=",kt,sep = " "))
#print(table(GTindex))

                #Recoding the optimum KI
                if(count==1){
                  KI.save=KI
                  LL.save=p3d$REMLs
                }else{
                  if(p3d$REMLs<LL.save){
                    KI.save=KI
                    LL.save=p3d$REMLs
                  }
                }

#print(paste("CA is ",ca))
#print(paste("group is ",group))
#print(paste("kt is ",kt))

                #recording Compression profile on array
                Compression[count,1]=kt
                Compression[count,2]=ca
                Compression[count,3]=group
                Compression[count,4]=p3d$REMLs
                Compression[count,5]=p3d$vgs
                Compression[count,6]=p3d$ves
#print("result saved")

              }else{# end of if(!byPass)

                #Set QTNs
                if(count==1)print("-------The burger is SNP-----------------------------------")
  #bin.size=bin
  #inclosure.size=inc

 
                #@@@This section is not useful
                if(!is.null(GP)){
  #print("Being specific...")

                  myGenotype <- GAPIT.Genotype(
                    G=NULL,
                    GD=NULL,
                    GM=GI,
                    KI=NULL,
                    kinship.algorithm="SUPER",
                    PCA.total=0,
                    SNP.fraction=SNP.fraction,
                    SNP.test=SNP.test,
                    file.path=file.path,
                    file.from=file.from, 
                    file.to=file.to, 
                    file.total=file.total, 
                    file.fragment = file.fragment, 
                    file.G=file.G, 
                    file.Ext.G=file.Ext.G,
                    file.GD=file.GD, 
                    file.GM=file.GM, 
                    file.Ext.GD=file.Ext.GD,
                    file.Ext.GM=file.Ext.GM,
                    SNP.MAF=SNP.MAF,
                    FDR.Rate = FDR.Rate,
                    SNP.FDR=SNP.FDR,
                    SNP.effect=SNP.effect,
                    SNP.impute=SNP.impute,
                    LD.chromosome=LD.chromosome,
                    LD.location=LD.location,
                    LD.range=LD.range,
                    kinship.cluster=kinship.cluster,#NJtree.group=NJtree.group,NJtree.type=NJtree.type,
                    GP=GP,
                    GK=NULL,
                    bin.size=bin,
                    inclosure.size=inc,
                    SNP.CV=SNP.CV,
                    GTindex=GTindex,
                    sangwich.top=NULL,
                    sangwich.bottom=sangwich.bottom,
                    Timmer = Timmer, 
                    Memory = Memory,
                    file.output=file.output, 
                    Create.indicator = Create.indicator, 
                    Major.allele.zero = Major.allele.zero
                    )
    
                  Timmer=myGenotype$Timmer
                  Memory=myGenotype$Memory
  
                  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Genotype for burger")
                  Memory=GAPIT.Memory(Memory=Memory,Infor="Genotype for burger")
  
                  print(paste("bin---",bin,"---inc---",inc,sep=""))
                  GK=GD[GTindex,myGenotype$SNP.QTN]
                  SUPER_GD=GD[,myGenotype$SNP.QTN]
                  SNPVar=apply(as.matrix(GK), 2, stats::var)
  
                  GK=GK[,SNPVar>0]
                  SUPER_GD=SUPER_GD[,SNPVar>0]
                  GK=cbind(as.data.frame(GT[GTindex]),as.data.frame(GK)) #add taxa
  # print(length(GT))
  # print(dim(SUPER_GD))
                  SUPER_GD=cbind(as.data.frame(GT[GTindex]),as.data.frame(SUPER_GD)) #add taxa
# print(dim(GK))
  #GP=NULL
                }# end of if(is.null(GK)) 


if(!is.null(GK) & numSetting>1)
{
print("-------Calculating likelihood-----------------------------------")
 # myBurger=GAPIT.Burger(Y=Y,CV=CV,GK=GK)
    myBurger=GAPIT.Burger(Y=Y,CV=NULL,GK=GK)   #########modified by Jiabo Wang

  myREML=myBurger$REMLs
  myVG=myBurger$vg
  myVE=myBurger$ve
}else{
  myREML=NA
  myVG=NA
  myVE=NA
}

#Recoding the optimum GK
if(count==1){
  GK.save=GK
  LL.save=myREML
  	SUPER_optimum_GD=SUPER_GD     ########### get SUPER GD

}else{
  if(myREML<LL.save){
    GK.save=GK
    LL.save=myREML
	SUPER_optimum_GD=SUPER_GD     ########### get SUPER GD
  }
}
  
#Put to storage
Compression[count,1]=1
Compression[count,2]=bin
Compression[count,3]=inc
Compression[count,4]=myREML
Compression[count,5]=myVG
Compression[count,6]=myVG
print(Compression[count,]) 

#print("---------------SUPER 2nd stage: calculating LL ------------------------")


}   # end of if(byPass)

}#end of for (ca in kinship.cluster)

#Skip the rest group in case group 1 is finished
if(group==1) break #To skip the rest group interations

}#end of for (group in GROUP)
}#end of for (kt in kinship.group)

  
}#end of for (inc in inclosure)
}#end of for (bin in bin.level)


if(Model.selection == TRUE){ 

  print("------------------------Model selection for optimal number of PCs and Covariates-------------------------------------------------")
  #update KI with the best likelihood
  KI=KI.save
  if(numSetting>1){
  Compression=Compression[order(as.numeric(Compression[,4]),decreasing = FALSE),]  #sort on REML
  kt=Compression[1,1]
  ca=Compression[1,2]
  group=Compression[1,3]
  }

  cp <- GAPIT.Compress(KI=KI,kinship.cluster=ca,kinship.group=kt,GN=group,Timmer=Timmer,Memory=Memory)
  Timmer=cp$Timmer
  Memory=cp$Memory

  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="PreP3D 2_cp")
  Memory=GAPIT.Memory(Memory=Memory,Infor="PreP3D 2_cp")
  
  bk <- GAPIT.Block(Z=Z,GA=cp$GA,KG=cp$KG)
  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="PreP3D 2_bk")
  Memory=GAPIT.Memory(Memory=Memory,Infor="PreP3D 2 bk")

  zc <- GAPIT.ZmatrixCompress(Z=Z,GAU =bk$GA)

  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="PreP3D 2_zc")
  Memory=GAPIT.Memory(Memory=Memory,Infor="PreP3D 2 zc")

  z0=as.matrix(zc$Z[,-1])
  Z1=matrix(as.numeric(z0),nrow=nrow(z0),ncol=ncol(z0))


  
  BIC <- rep(NA,ncol(X0))
  LogLike <- rep(NA, ncol(X0))
  for(i in 1:ncol(X0)){#1 because the first column of X0 is the intercept

    X0.test <- as.matrix(X0[,1:i]) 
    
    #print("The dim of bk$KW is ")
    #print(dim(bk$KW))
    #print(dim(X0.test))
    #print(dim(CVI))

    p3d <- GAPIT.EMMAxP3D(ys=ys,xs=as.matrix(as.data.frame(GD[,1])),K = as.matrix(bk$KW) ,Z=Z1,X0=X0.test,CVI=CVI,CV.Inheritance=CV.Inheritance,GI=GI,SNP.P3D=SNP.P3D,Timmer=Timmer,Memory=Memory,fullGD=fullGD,
            SNP.permutation=SNP.permutation, GP=GP,
			      file.path=file.path,file.from=file.from,file.to=file.to,file.total=file.total, file.fragment = file.fragment, byFile=byFile, file.G=file.G,file.Ext.G=file.Ext.G,file.GD=file.GD, file.GM=file.GM, file.Ext.GD=file.Ext.GD,file.Ext.GM=file.Ext.GM,
            GTindex=GTindex,genoFormat=genoFormat,optOnly=TRUE,SNP.effect=SNP.effect,SNP.impute=SNP.impute,name.of.trait=name.of.trait, Create.indicator = Create.indicator, Major.allele.zero = Major.allele.zero)

    
    
    k.num.param <- 2+i
    #k is (i-1) because we have the following parameters in the likelihood function:
    #  intercept
    #  (i-1) covariates
    #  sigma_g
    #  delta
    
    #print(paste("The value of round(p3d$REMLs,5) is ", round(p3d$REMLs,5), sep = ""))
    #print(paste("The value of log(GTindex) is ", log(GTindex), sep = ""))
    #print(paste("The value of 0.5*k.num.param*log(GTindex) is ", 0.5*k.num.param*log(nrow(Z1)), sep = ""))
    
    LogLike[i] <- p3d$logLM
    BIC[i] <- p3d$logLM -(0.5*k.num.param*log(nrow(Z1)))
    
    #print("The value of k.num.param  is: ")
    #print(k.num.param)
    
    #print(paste("The value of nrow(Z1) is ", nrow(Z1), sep = ""))  
    
    }   
    Optimum.from.BIC <- which(BIC == max(BIC))
    
    print(paste("-----------------------The optimal number of PCs/covariates is ", (Optimum.from.BIC-1)," -------------------------", sep = ""))
    
    BIC.Vector <- cbind(as.matrix(rep(0:(ncol(X0)-1))), as.matrix(BIC), as.matrix(LogLike))

           
    #print(seq(0:ncol(X0)))
    
       #print(BIC.Vector)
 
    colnames(BIC.Vector) <- c("Number of PCs/Covariates", "BIC (larger is better) - Schwarz 1978", "log Likelihood Function Value")
    
    utils::write.table(BIC.Vector, paste("GAPIT.", name.of.trait, ".BIC.Model.Selection.Results.csv", sep = ""), quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
    
    #print(BIC.Vector)
    
    X0 <- X0[,1:(Optimum.from.BIC)]
    
    if(Optimum.from.BIC == 1){
    X0 <- as.matrix(X0)
    }
    print("The dimension of X0 after model selection is:")
    print(dim(X0))
    
    print("The head of X0 after model selection is")
    print(utils::head(X0))
    

} # where does it start: 522

print("---------------------Sandwich bottom bun-------------------------------")
print("Compression") 
print(Compression)

#Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Compression")
#Memory=GAPIT.Memory(Memory=Memory,Infor="Copmression")

if(numSetting==1)
{
  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="GWAS")
  Memory=GAPIT.Memory(Memory=Memory,Infor="GWAS")
}
  
#Perform GWAS with the optimum setting
#This section is omited if there is only one setting
if((numSetting>1)| (!is.null(sangwich.bottom)&!byPass) | Model.selection) {
  print("Genomic screening..." )
  
optOnly=FALSE  #set default to false and change it to TRUE in these situations:
if(!hasGenotype) optOnly=TRUE
if(!SNP.test) optOnly=TRUE

if(optOnly){
 colInclude=1
}else{
 colInclude=c(1:ncol(GD))
}

if(numSetting>1){
#Find the best ca,kt and group
#print(paste(as.numeric(Compression[1,4]))) ###added by Jiabo Wang 2015.7.20
#print(paste(min(as.numeric(Compression[,4]),rm.na=TRUE)))
adjust_value=as.numeric(Compression[1,4])-min(as.numeric(Compression[,4]),rm.na=TRUE)
# nocompress_value=as.numeric(Compression[1,4])
# REML_storage=as.numeric(Compression[,4])

adjust_sq=sqrt(stats::var(as.numeric(Compression[,4])))
# threshold=adjust_mean*0.1       
if(which.min(as.numeric(Compression[,4]))!=1)     ###added by Jiabo Wang 2015.7.20
{
if(which.min(as.numeric(Compression[,4]))==which.max(as.numeric(Compression[,5])))
{
  kt=Compression[which.min(as.numeric(Compression[,4])),1]
  ca=Compression[which.min(as.numeric(Compression[,4])),2]
  group=Compression[which.min(as.numeric(Compression[,4])),3]
  va=Compression[which.min(as.numeric(Compression[,4])),5]
  ve=Compression[which.min(as.numeric(Compression[,4])),6]
}else{
  # Compression0=Compression
  cnn=which.min(as.numeric(Compression[,4]))
  if(cnn-which.min(as.numeric(Compression[-cnn,4]))<2)
    {
      kt=Compression[which.min(as.numeric(Compression[,4])),1]
      ca=Compression[which.min(as.numeric(Compression[,4])),2]
      group=Compression[which.min(as.numeric(Compression[,4])),3]
      va=Compression[which.min(as.numeric(Compression[,4])),5]
      ve=Compression[which.min(as.numeric(Compression[,4])),6]
    }else{
      kt=Compression[1,1]
      ca=Compression[1,2]
      group=Compression[1,3]
      va=Compression[1,5]
      ve=Compression[1,6]
      print("The difference of compression is not enough!!")
    }
  
}



# Compression=Compression0

print(paste("Compress Optimum: ",ca,kt,group,va,va,ve,sep = " "))
}else{
Compression=Compression[order(as.numeric(Compression[,4]),decreasing = FALSE),]  #sort on REML

kt=Compression[1,1]
ca=Compression[1,2]
group=Compression[1,3]
print(paste("Optimum: ",Compression[1,2],Compression[1,1],Compression[1,3],Compression[1,5],Compression[1,6],Compression[1,4],sep = " "))
}
}#end  if(numSetting>1)
Compression=Compression[order(as.numeric(Compression[,4]),decreasing = FALSE),]
print(Compression)


print("--------------  Sandwich bottom ------------------------") 

if(!byPass) 
{ 
print("--------------  Sandwich bottom with raw burger------------------------") 

 if(Model.selection == FALSE){
  #update KI with the best likelihood
  if(is.null(sangwich.bottom)) KI=KI.save

  cp <- GAPIT.Compress(KI=KI,kinship.cluster=ca,kinship.group=kt,GN=group,Timmer=Timmer,Memory=Memory)
  Timmer=cp$Timmer
  Memory=cp$Memory
  
  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="PreP3D 2_cp")
  Memory=GAPIT.Memory(Memory=Memory,Infor="PreP3D 2_cp")
  
  bk <- GAPIT.Block(Z=Z,GA=cp$GA,KG=cp$KG)
  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="PreP3D 2_bk")
  Memory=GAPIT.Memory(Memory=Memory,Infor="PreP3D 2 bk")
  
  zc <- GAPIT.ZmatrixCompress(Z=Z,GAU =bk$GA)
  
  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="PreP3D 2_zc")
  Memory=GAPIT.Memory(Memory=Memory,Infor="PreP3D 2 zc")
  
  #Reform KW and Z into EMMA format
  
  z0=as.matrix(zc$Z[,-1])   
  Z1=matrix(as.numeric(z0),nrow=nrow(z0),ncol=ncol(z0))
 }
 
 print("--------------EMMAxP3D with the optimum setting-----------------------") 
 #print(dim(ys))
 #print(dim(as.matrix(as.data.frame(GD[GTindex,colInclude]))))
  p3d <- GAPIT.EMMAxP3D(ys=ys,xs=as.matrix(as.data.frame(GD[GTindex,colInclude]))   ,K = as.matrix(bk$KW) ,Z=Z1,X0=as.matrix(X0),CVI=CVI, CV.Inheritance=CV.Inheritance,GI=GI,SNP.P3D=SNP.P3D,Timmer=Timmer,Memory=Memory,fullGD=fullGD,
          SNP.permutation=SNP.permutation, GP=GP,
    			 file.path=file.path,file.from=file.from,file.to=file.to,file.total=file.total, file.fragment = file.fragment, byFile=byFile, file.G=file.G,file.Ext.G=file.Ext.G,file.GD=file.GD, file.GM=file.GM, file.Ext.GD=file.Ext.GD,file.Ext.GM=file.Ext.GM,
           GTindex=GTindex,genoFormat=genoFormat,optOnly=optOnly,SNP.effect=SNP.effect,SNP.impute=SNP.impute,name.of.trait=name.of.trait, Create.indicator = Create.indicator, Major.allele.zero = Major.allele.zero)  
    
  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="GWAS")
  Memory=GAPIT.Memory(Memory=Memory,Infor="GWAS")  
 print("--------------EMMAxP3D with the optimum setting done------------------") 
  
}#end of if(!byPass) 
}#end of if(numSetting>1 & hasGenotype & !SNP.test)  

#print("Screening wiht the optimum setting done") 

if(byPass)
{
print("---------------Sandwich bottom with grilled burger---------------------") 
print("---------------Sandwich bottom: reload bins ---------------------------")

#SUPER: Final screening
  GK=GK.save
  # print(GK)
  myBread=GAPIT.Bread(Y=Y,CV=CV,Z=Z,GK=GK,GD=cbind(as.data.frame(GT[GTindex]),as.data.frame(GD)),GM=GI,method=sangwich.bottom,GTindex=GTindex,LD=LD,file.output=file.output)
  
  print("SUPER saving results...")

  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="GWAS")
  Memory=GAPIT.Memory(Memory=Memory,Infor="GWAS")  

   
}   #end of if(byPass)

print("--------------------Final results presentations------------------------")



#Plotting optimum group kinship
if(!byPass) {
  if(length(bk$KW)>1 &length(bk$KW)<length(KI) & length(bk$KW)<1000 & GAPIT3.output){
#    if( file.output == TRUE ){
        grDevices::pdf(paste("GAPIT.",name.of.trait,".Kin.Optimum.pdf",sep=""), width = 12, height = 12)
        graphics::par(mar = c(25,25,25,25))
        gplots::heatmap.2(as.matrix(bk$KW),  cexRow =.2, cexCol = 0.2, col=rev(grDevices::heat.colors(256)), scale="none", symkey=FALSE, trace="none")
        grDevices::dev.off()
#    }
  }
}


#Merge GWAS resultss from files to update ps,maf and nobs in p3d
if(byFile&!fullGD)
{
print("Loading GWAS results from file...")
for (file in file.from:file.to)
{

#Initicalization
frag=1
numSNP=file.fragment

while(numSNP==file.fragment) {     #this is problematic if the read end at the last line  

#Initicalization GI to detect reading empty line
#theGI=NULL
#theP=NULL
#theMAF=NULL
#thenobs=NULL

 
#reload results from files
print(paste("Current file ",file,"Fragment: ",frag))

theGI <- try(utils::read.table(paste("GAPIT.TMP.GI.",name.of.trait,file,".",frag,".txt",sep=""), head = TRUE)   ,silent=TRUE)
theP <- try(utils::read.table(paste("GAPIT.TMP.ps.",name.of.trait,file,".",frag,".txt",sep=""), head = FALSE)   ,silent=TRUE)
theMAF <- try(utils::read.table(paste("GAPIT.TMP.maf.",name.of.trait,file,".",frag,".txt",sep=""), head = FALSE),silent=TRUE)
thenobs <- try(utils::read.table(paste("GAPIT.TMP.nobs.",name.of.trait,file,".",frag,".txt",sep=""),head= FALSE),silent=TRUE)
thersquare_base <- try(utils::read.table(paste("GAPIT.TMP.rsquare.base.",name.of.trait,file,".",frag,".txt",sep=""),head= FALSE),silent=TRUE)
thersquare <- try(utils::read.table(paste("GAPIT.TMP.rsquare.",name.of.trait,file,".",frag,".txt",sep=""),head= FALSE),silent=TRUE)
          thedf  <- try(utils::read.table(paste("GAPIT.TMP.df.",name.of.trait,file,".",frag,".txt",sep=""),head= FALSE),silent=TRUE)
          thetvalue  <- try(utils::read.table(paste("GAPIT.TMP.tvalue.",name.of.trait,file,".",frag,".txt",sep=""),head= FALSE),silent=TRUE)
          thestderr  <- try(utils::read.table(paste("GAPIT.TMP.stderr.",name.of.trait,file,".",frag,".txt",sep=""),head= FALSE),silent=TRUE)
theeffect.est <- try(utils::read.table(paste("GAPIT.TMP.effect.est.",name.of.trait,file,".",frag,".txt",sep=""),head= FALSE),silent=TRUE)

if(inherits(theGI, "try-error"))  {
#if(nrow(theGI)<1){
  numSNP=0
  #print("This fragment is empty.")
}else{



#print("Records loaded for this fragment.")
  numSNP=nrow(theGI)  
  colnames(theP)="P"
  colnames(theMAF )="MAF"
  colnames(thenobs )="nobs"
  colnames(thersquare_base) = "Base.Model.R.square"  
  colnames(thersquare) = "Model.R.square"
            colnames(thedf) = "Model.DF"
            colnames(thetvalue) = "Model.tvalue"
            colnames(thestderr) = "Model.stderr"
  colnames(theeffect.est) = "Effect.Est"    
  colnames(theGI) = colnames(GI)
 



#Merge results  
  if(file==file.from & frag==1){

    GI=theGI  
    #print(dim(GI))
    allP=theP
    #print(head(theP))
    allMAF=theMAF
    allnobs=thenobs
    allrsquare_base=thersquare_base
    allrsquare=thersquare
              alldf=thedf
              alltvalue=thetvalue
              allstderr=thestderr
    alleffect.est=theeffect.est

  }else{
    allP=as.data.frame(rbind(as.matrix(allP),as.matrix(theP))  )
    allMAF=as.data.frame(rbind(as.matrix(allMAF),as.matrix(theMAF)) )
    allnobs=as.data.frame(rbind(as.matrix(allnobs),as.matrix(thenobs)))
    allrsquare_base=as.data.frame(rbind(as.matrix(allrsquare_base),as.matrix(thersquare_base)))
    allrsquare=as.data.frame(rbind(as.matrix(allrsquare),as.matrix(thersquare)))
              alldf=as.data.frame(rbind(as.matrix(alldf),as.matrix(thedf)))
              alltvalue=as.data.frame(rbind(as.matrix(alltvalue),as.matrix(thetvalue)))
              allstderr=as.data.frame(rbind(as.matrix(allstderr),as.matrix(thestderr)))
    alleffect.est=as.data.frame(rbind(as.matrix(alleffect.est),as.matrix(theeffect.est)))
    #print("!!!!!!!!!!!!!!!")
    #print(dim(GI))
    #print(dim(theGI))
    GI=as.data.frame(rbind(as.matrix(GI),as.matrix(theGI)))
  }

}#end of  if(inherits(theGI, "try-error")) (else section)

#setup for next fragment
frag=frag+1   #Progress to next fragment 

}#end of loop on fragment: while(numSNP==file.fragment)
}#end of loop on file

#update p3d with components from files

  p3d$ps=allP
  p3d$maf=allMAF
  p3d$nobs=allnobs
  p3d$rsquare_base=allrsquare_base
  p3d$rsquare=allrsquare
      p3d$df=alldf
      p3d$tvalue=alltvalue
      p3d$stderr=allstderr
  p3d$effect.est=alleffect.est
  
#Delete all the GAPIT.TMP files
theFile=paste("GAPIT.TMP.",name.of.trait,".*")
  system('cmd /c del "GAPIT.TMP*.*"') 
  system('cmd /c del "GAPIT.TMP*.*"') 
  print("GWAS results loaded from all files succesfully!")
} #end of if(byFile)

#--------------------------------------------------------------------------------------------------------------------#
#Final report   
print("Generating summary" )
GWAS=NULL
GPS=NULL
rm(zc)
gc()

Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Final")
Memory=GAPIT.Memory(Memory=Memory,Infor="Final")

#genomic prediction
print("Genomic Breeding Values (GBV) ..." )
#print(p3d$BLUP)
gs=NULL
if(!byPass) 
{

if(length(bk$KW)>ncol(X0)) {
    gs <- GAPIT.GS(KW=bk$KW,KO=bk$KO,KWO=bk$KWO,GAU=bk$GAU,UW=cbind(p3d$BLUP,p3d$PEV))
}

print("Writing GBV and Acc..." )

GPS=NULL
if(length(bk$KW)>ncol(X0)) GPS=gs$BLUP
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="GPS")
Memory=GAPIT.Memory(Memory=Memory,Infor="GPS")

#Make heatmap for distribution of BLUP and PEV
print("GBV and accuracy distribution..." )
if(length(bk$KW)>ncol(X0) &file.output) {
  GAPIT.GS.Visualization(gsBLUP = gs$BLUP, BINS=BINS,name.of.trait = name.of.trait)
}

#Make a plot Summarzing the Compression Results, if more than one "compression level" has been assessed
print("Compression portfolios..." )
#print(Compression)
if(file.output){
  GAPIT.Compression.Visualization(Compression = Compression, 
                                  name.of.trait = name.of.trait)
}
print("Compression Visualization done")

if(length(Compression)<1){
  h2.opt= NULL
}else{
  print(Compression)
if(length(Compression)<6) Compression=t(as.matrix(Compression[which(Compression[,4]!="NULL" | Compression[,4]!="NaN"),]))
if(length(Compression)==6) Compression=matrix(Compression,1,6) 
if(length(Compression)>6) Compression=Compression[which(Compression[,4]!="NULL" | Compression[,4]!="NaN"),]
Compression.best=Compression[1,] 
variance=as.numeric(Compression.best[5:6])
varp=variance/sum(variance)
h2.opt= varp[1]
}

Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Compression.Visualization")
Memory=GAPIT.Memory(Memory=Memory,Infor="Compression.Visualization")
# print("$$$$$")
# print(str(p3d))

ps=p3d$ps
nobs=p3d$nobs
maf=p3d$maf
rsquare_base=p3d$rsquare_base
rsquare=p3d$rsquare
      df=p3d$df
      tvalue=p3d$tvalue
      stderr=p3d$stderr
effect.est=p3d$effect.est
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Extract p3d results")
Memory=GAPIT.Memory(Memory=Memory,Infor="Extract p3d results")
print("p3d objects transfered")  

#where does it start: 936
}else{  #byPass
    #print("The head of myBread$GWAS is")
  #print(head(myBread$GWAS))
  GPS=myBread$BLUP
  ps=myBread$GWAS[,4]
  nobs=myBread$GWAS[,6]
  #print(dim(GI))
  #print(head())
  Bread_index=match(as.character(myBread$GWAS[,1]),as.character(GI[,1]))
  #print(GD[1:5,1:5])
  Bread_X=GD[,Bread_index]
  #print(dim(Bread_X))
  maf=apply(Bread_X,2,function(one) abs(1-sum(one)/(2*nrow(Bread_X))))
  maf[maf>0.5]=1-maf[maf>0.5]
  rsquare_base=rep(NA,length(ps))
  rsquare=rep(NA,length(ps))
  df=rep(NA,length(nobs))
  tvalue=rep(NA,length(nobs))
  stderr=rep(NA,length(nobs))
  effect.est=rep(NA,length(nobs))
  
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Extract bread results")
Memory=GAPIT.Memory(Memory=Memory,Infor="Extract bread results")
 
}
print("Merge BLUP and BLUE")
#print(head(ps))
#Merge BLUP and BLUE
Pred=NULL
if((!byPass)&(!Model.selection)){
 print("GAPIT before BLUP and BLUE")
 #print(dim(p3d$BLUE))
 BLUE=data.frame(cbind(data.frame(CV.taxa),data.frame(p3d$BLUE)))
 colnames(BLUE)=c("Taxa","BLUE")
 
 #Initial BLUP as BLUe and add additional columns
 gs.blup=cbind(BLUE,NA,NA,0,NA)
 
 if(!is.null(gs))gs.blup=gs$BLUP
 BB= merge(gs.blup, BLUE, by.x = "Taxa", by.y = "Taxa")
 if (is.null(my_allCV)){my_allX=matrix(1,length(my_taxa),1)
 }else{
     # my_allX=as.matrix(my_allCV[,-1])
     my_allX=cbind(1,as.matrix(my_allCV[,-1]))
	}
	
    #print(dim(my_allX))
    #print(head(my_allX))
    #print(dim(BB))
    #print(CV.Inheritance)
 if(is.null(CV.Inheritance))
 
   {
   Prediction=BB[,5]+BB[,7]
   Pred_Heritable=Prediction
   }
 if(!is.null(CV.Inheritance))
   {
       #inher_CV=my_allX[,1:(1+CV.Inheritance)]
       #beta.Inheritance=p3d$effect.cv[1:(1+CV.Inheritance)]
    #print(beta.Inheritance)
    #if(length(beta)==1)CV=X
    all_BLUE=try(my_allX%*%p3d$effect.cv,silent=T)
    if(inherits(BLUE, "try-error")) all_BLUE = NA
    

    Pred_Heritable=BB[,5]+BB[,7]
    Prediction=BB[,5]+all_BLUE
   }
   #print("@@@@@@@@@@")
 #print(dim(CVI))
 #print(BB)
 #CV.Inheritance
 #Pred_Heritable=p3d$effect.cv[CV.Inheritance]%*%CVI[CV.Inheritance]+BB[,7]
 Pred=data.frame(cbind(BB,data.frame(Prediction)),data.frame(Pred_Heritable))
 if(noCV)
    {
    if(NOBLUP)
    {Pred=NA
    }else{
    BLUE=Pred$BLUE[1]
    prediction=as.matrix(GPS$BLUP)+(BLUE)
    Pred=cbind(GPS,BLUE,prediction)
 colnames(Pred)=c("Taxa","Group","RefInf","ID","BLUP","PEV","BLUE","Prediction")
    }#end NOBLUP
    }#end noCV
 print("GAPIT after BLUP and BLUE")
}

#Export BLUP and PEV
if(!byPass &GAPIT3.output) 
{
print("Exporting BLUP and Pred")
  #try(write.table(gs$BLUP, paste("GAPIT.", name.of.trait,".BLUP.csv" ,sep = ""), quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE))
  try(utils::write.table(Pred, paste("GAPIT.", name.of.trait,".PRED.csv" ,sep = ""), quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE))
}

if(byPass) 
{
  theK.back=NULL
}else{
  theK.back=cp$KG
}
if(byPass)Compression[1,4]=0 #create a fake value to aloow output of SUPER 

#Export GWAS results
PWI.Filtered=NULL
if(hasGenotype &SNP.test &!is.na(Compression[1,4]))     #require not NA REML 
{
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Extract GWAS start")
Memory=GAPIT.Memory(Memory=Memory,Infor="Extract GWAS start")


  #print("Filtering SNPs with MAF..." )
	#index=maf>=SNP.MAF	   
  
	PWI.Filtered=cbind(GI,ps,maf,nobs,rsquare_base,rsquare,effect.est)#[index,]
	#print(dim(PWI.Filtered))
	colnames(PWI.Filtered)=c("SNP","Chromosome","Position ","P.value", "maf", "nobs", "Rsquare.of.Model.without.SNP","Rsquare.of.Model.with.SNP","effect")

if(!byPass){  
   if(Create.indicator){
    #Add a counter column for GI
    GI.counter <- cbind(GI, seq(1:nrow(GI))) 
    
    #Turn GI and effect.est into data frames
    GI.counter.data.frame <- data.frame(GI.counter)
    colnames(GI.counter.data.frame) <- c("X1", "X2", "X3", "X4")
    
    effect.est.data.frame <- data.frame(effect.est)
    colnames(effect.est.data.frame) <- c("X1", "X2", "X3")
    print(utils::head(GI.counter.data.frame))
    print(utils::head(effect.est.data.frame))
    #Do a merge statement
    GWAS.2 <- merge(GI.counter.data.frame, effect.est.data.frame, by.x = "X4", by.y = "X1")
    
    #Remove the counter column
    GWAS.2 <- GWAS.2[,-1]
    
    #Add column names
    colnames(GWAS.2) <- c("SNP","Chromosome","Position ", "Genotype", "Allelic Effect Estimate")
    
    
   }
   if(!Create.indicator){ 
    GWAS.2 <- PWI.Filtered[,c(1:3,9)]
    colnames(GWAS.2) <- c("SNP","Chromosome","Position ", "Allelic Effect Estimate")
   } 
}
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="MAF filtered")
Memory=GAPIT.Memory(Memory=Memory,Infor="MAF filtered")
		     
  #print("SNPs filtered with MAF")
   
  
  if(!is.null(PWI.Filtered))
  {

  #Run the BH multiple correction procedure of the results
  #Create PWIP, which is a table of SNP Names, Chromosome, bp Position, Raw P-values, FDR Adjusted P-values
  #print("Calculating FDR..." )

  PWIP <- GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure(PWI = PWI.Filtered, FDR.Rate = FDR.Rate, FDR.Procedure = "BH")
  
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Multiple Correction")
Memory=GAPIT.Memory(Memory=Memory,Infor="Multiple Correction")


  #QQ plots
  #print("QQ plot..." )
  if(file.output) GAPIT.QQ(P.values = PWIP$PWIP[,4], name.of.trait = name.of.trait,DPP=DPP)


Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="QQ plot")
Memory=GAPIT.Memory(Memory=Memory,Infor="QQ plot")


  #Manhattan Plots
  
  
   #print("Manhattan plot (Genomewise)..." )
#  if(file.output) GAPIT.Manhattan(GI.MP = PWIP$PWIP[,2:4], name.of.trait = name.of.trait, DPP=DPP, plot.type = "Genomewise",cutOff=cutOff)
#  if(file.output) GAPIT.Manhattan(GI.MP = PWIP$PWIP[,2:4], name.of.trait = name.of.trait, DPP=DPP, plot.type = "Genomewise",cutOff=cutOff,seqQTN=QTN.position)  #QTN does not work with sorted P
 if(file.output) GAPIT.Manhattan(GI.MP = PWI.Filtered[,2:4], name.of.trait = name.of.trait, DPP=DPP, plot.type = "Genomewise",cutOff=cutOff,seqQTN=QTN.position,plot.style=plot.style,plot.bin=plot.bin,chor_taxa=chor_taxa)

 #print("Manhattan plot (Chromosomewise)..." )
 
  #if(file.output) GAPIT.Manhattan(GI.MP = PWIP$PWIP[,2:4], name.of.trait = name.of.trait, DPP=DPP, plot.type = "Chromosomewise",cutOff=cutOff)
 if(file.output&SNP.fraction==1) GAPIT.Manhattan(GI.MP = PWI.Filtered[,2:4],GD=GD,CG=CG, name.of.trait = name.of.trait, DPP=DPP, plot.type = "Chromosomewise",cutOff=cutOff,plot.bin=plot.bin,chor_taxa=chor_taxa)

Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Manhattan plot")
Memory=GAPIT.Memory(Memory=Memory,Infor="Manhattan plot")


  #Association Table
  #print("Association table..." )
  #print(dim(PWIP$PWIP))
  #GAPIT.Table(final.table = PWIP$PWIP, name.of.trait = name.of.trait,SNP.FDR=SNP.FDR)
  #print(head(PWIP$PWIP))
  GWAS=PWIP$PWIP[PWIP$PWIP[,9]<=SNP.FDR,]
  #print("Joining tvalue and stderr" )
  
        DTS=cbind(GI,df,tvalue,stderr,effect.est)
        colnames(DTS)=c("SNP","Chromosome","Position","DF","t Value","std Error","effect")	

  #print("Creating ROC table and plot" )
	if(file.output) myROC=GAPIT.ROC(t=tvalue,se=stderr,Vp=stats::var(ys),trait=name.of.trait)
  #print("ROC table and plot created" )

  #MAF plots
  #print("MAF plot..." )
   if(file.output) myMAF1=GAPIT.MAF(MAF=GWAS[,5],P=GWAS[,4],E=NULL,trait=name.of.trait)


  #print(dim(GWAS))

  if(file.output){
   utils::write.table(GWAS, paste("GAPIT.", name.of.trait, ".GWAS.Results.csv", sep = ""), quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
   utils::write.table(DTS, paste("GAPIT.", name.of.trait, ".Df.tValue.StdErr.csv", sep = ""), quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
   if(!byPass) utils::write.table(GWAS.2, paste("GAPIT.", name.of.trait, ".Allelic_Effect_Estimates.csv", sep = ""), quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
  }


  
  } #end of if(!is.null(PWI.Filtered))
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Extract GWAS end")
Memory=GAPIT.Memory(Memory=Memory,Infor="Extract GWAS end")

  
} #end of if(hasGenotype )

#Log
if(GAPIT3.output) log=GAPIT.Log(Y=Y,KI=KI,Z=Z,CV=CV,SNP.P3D=SNP.P3D,
				group.from = group.from ,group.to =group.to ,group.by = group.by ,kinship.cluster = kinship.cluster, kinship.group= kinship.group,
                      	ngrid = ngrid , llin = llin , ulim = ulim , esp = esp ,name.of.trait = name.of.trait)
#Memory usage
#GAPIT.Memory.Object(name.of.trait=name.of.trait)

#Timming
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Report")
Memory=GAPIT.Memory(Memory=Memory,Infor="Report")
if(file.output){
file=paste("GAPIT.", name.of.trait,".Timming.csv" ,sep = "")
utils::write.table(Timmer, file, quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)

file=paste("GAPIT.", name.of.trait,".Memory.Stage.csv" ,sep = "")
utils::write.table(Memory, file, quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
}
print(paste(name.of.trait, "has been analyzed successfully!") )
print(paste("The results are saved in the directory of ", getwd()) )



#print("==========================================================================================")
TV<-list()
TV$ps=ps
TV$nobs=nobs
TV$maf=maf
TV$rsquare_base=rsquare_base
TV$rsquare=rsquare
TV$df=df
TV$tvalue=tvalue
TV$stderr=stderr
TV$effect.est=effect.est
#print("!!!!!!!!!!!!!")
#print(head(effect.est))
#print(head(DTS[,7]))
#print(ys)
if(byPass | Model.selection) Pred <- NA
print("before ending GAPIT.Main")
#print(dim(Compression))
return (list(Timmer=Timmer,Compression=Compression,kinship.optimum=theK.back, kinship=KI,PC=PC,GWAS=PWI.Filtered, GPS=GPS,Pred=Pred,REMLs=Compression[count,4],Timmer=Timmer,Memory=Memory,SUPER_GD=SUPER_GD,P=ps,effect.snp=DTS[,7],effect.cv=p3d$effect.cv,h2= h2.opt,TV=TV))
} #end if non-SUPER.GS situation, this is a long if statement, structure needs improvement
}#The function GAPIT.Main ends here
#=============================================================================================


