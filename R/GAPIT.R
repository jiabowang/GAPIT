
#' GAPIT Genome Association and Prediction Integrated Tools
#' 
#' @description 
#' GAPIT analyzes phenotypic and genotypics data to infer association.
#' 
#' 
#' @param Y = NULL,
#' @param G = NULL,
#' @param GD = NULL,
#' @param GM = NULL,
#' @param KI = NULL,
#' @param Z = NULL,
#' @param CV = NULL,
#' @param CV.Inheritance = NULL,
#' @param GP = NULL,
#' @param GK = NULL,
#' @param testY = NULL,
#' @param group.from = 1e+06,
#' @param group.to = 1e+06,
#' @param group.by = 20,
#' @param DPP = 1e+05,
#' @param kinship.cluster = "average",
#' @param kinship.group = "Mean",
#' @param kinship.algorithm = "VanRaden",
#' @param buspred = FALSE,
#' @param lmpred = FALSE,
#' @param FDRcut = FALSE,
#' @param bin.from = 10000,
#' @param bin.to = 10000,
#' @param bin.by = 10000,
#' @param inclosure.from = 10,
#' @param inclosure.to = 10,
#' @param inclosure.by = 10,
#' @param SNP.P3D = TRUE,
#' @param SNP.effect = "Add",
#' @param SNP.impute = "Middle",
#' @param PCA.total = 0,
#' @param SNP.fraction = 1,
#' @param seed = NULL,
#' @param BINS = 20,
#' @param SNP.test = TRUE,
#' @param SNP.MAF = 0,
#' @param FDR.Rate = 1,
#' @param SNP.FDR = 1,
#' @param SNP.permutation = FALSE,
#' @param SNP.CV = NULL,
#' @param SNP.robust = "GLM",
#' @param file.from = 1,
#' @param file.to = 1,
#' @param file.total = NULL,
#' @param file.fragment = 99999,
#' @param file.path = NULL,
#' @param file.G = NULL,
#' @param file.Ext.G = NULL,
#' @param file.GD = NULL,
#' @param file.GM = NULL,
#' @param file.Ext.GD = NULL,
#' @param file.Ext.GM = NULL,
#' @param ngrid = 100,
#' @param llim = -10,
#' @param ulim = 10,
#' @param esp = 1e-10,
#' @param LD.chromosome = NULL,
#' @param LD.location = NULL,
#' @param LD.range = NULL,
#' @param PCA.col = NULL,
#' @param PCA.3d = FALSE,
#' @param NJtree.group = NULL,
#' @param NJtree.type = c("fan", "unrooted"),
#' @param sangwich.top = NULL,
#' @param sangwich.bottom = NULL,
#' @param QC = TRUE,
#' @param GTindex = NULL,
#' @param LD = 0.1,
#' @param plot.bin = 10^5,
#' @param file.output = TRUE,
#' @param cutOff = 0.05,
#' @param Model.selection = FALSE,
#' @param output.numerical = FALSE,
#' @param output.hapmap = FALSE,
#' @param Create.indicator = FALSE,
#' @param Multi_iter = FALSE,
#' @param num_regwas = 10,
#' @param opt = "extBIC",
#' @param QTN = NULL,
#' @param QTN.round = 1,
#' @param QTN.limit = 0,
#' @param QTN.update = TRUE,
#' @param QTN.method = "Penalty",
#' @param Major.allele.zero = FALSE,
#' @param Random.model = FALSE,
#' @param method.GLM = "FarmCPU.LM",
#' @param method.sub = "reward",
#' @param method.sub.final = "reward",
#' @param method.bin = "static",
#' @param bin.size = c(1e+06),
#' @param bin.selection = c(10, 20, 50, 100, 200, 500, 1000),
#' @param memo = NULL,
#' @param Prior = NULL,
#' @param ncpus = 1,
#' @param maxLoop = 3,
#' @param threshold.output = 0.01,
#' @param Inter.Plot = FALSE,
#' @param Inter.type = c("m", "q"),
#' @param WS = c(1, 1000, 10000, 1e+05, 1e+06, 1e+07),
#' @param alpha = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
#' @param maxOut = 100,
#' @param QTN.position = NULL,
#' @param CG = NULL,
#' @param converge = 1,
#' @param iteration.output = FALSE,
#' @param acceleration = 0,
#' @param iteration.method = "accum",
#' @param PCA.View.output = TRUE,
#' @param Geno.View.output = TRUE,
#' @param plot.style = "Oceanic",
#' @param SUPER_GD = NULL,
#' @param SUPER_GS = FALSE,
#' @param h2 = NULL,
#' @param NQTN = NULL,
#' @param QTNDist = "normal",
#' @param effectunit = 1,
#' @param category = 1,
#' @param r = 0.25,
#' @param cveff = NULL,
#' @param a2 = 0,
#' @param adim = 2,
#' @param Multiple_analysis = FALSE,
#' @param model = "MLM", options: MLM, GLM, CMLM, MMLM, SUPER, FarmCPU, gBLUP, or cBLUP
#' @param Para = NULL
#' 
#' @seealso 
#' GAPIT.DP(), GAPIT.Phenotype.View(), GAPIT.judge(), GAPIT.IC(), GAPIT.SS(), GAPIT.ID().
#' 
#' @author Zhiwu Zhang and Jiabo Wang
#' 
#' @examples 
#' \dontrun{
#' 
#' myPhenoFile <- system.file("extdata", "mdp_traits.txt.gz", package = "GAPIT3")
#' myGenoFile <- system.file("extdata", "mdp_genotype_test.hmp.txt.gz", package = "GAPIT3")
#' myPhenotypes <- read.table(myPhenoFile, header = TRUE)
#' myGenotypes  <- read.table(myGenoFile, header = FALSE)
#' 
#' myGAPIT <- GAPIT(
#'   Y = myPhenotypes,
#'   G = myGenotypes,
#'   PCA.total = 3,
#'   file.output = FALSE,
#'   model = "MLM"
#' )
#' }
#'
#'
#' @export
`GAPIT` <-
function(Y=NULL,G=NULL,GD=NULL,GM=NULL,KI=NULL,Z=NULL,CV=NULL,CV.Inheritance=NULL,GP=NULL,GK=NULL,testY=NULL,
 group.from=1000000 ,group.to=1000000,group.by=20,DPP=100000, 
 kinship.cluster="average", kinship.group='Mean',kinship.algorithm="VanRaden", buspred=FALSE,lmpred=FALSE,FDRcut=FALSE,
 bin.from=10000,bin.to=10000,bin.by=10000,inclosure.from=10,inclosure.to=10,inclosure.by=10,
 SNP.P3D=TRUE,SNP.effect="Add",SNP.impute="Middle",PCA.total=0, SNP.fraction = 1, seed = NULL, BINS = 20,SNP.test=TRUE,
 SNP.MAF=0,FDR.Rate = 1, SNP.FDR=1,SNP.permutation=FALSE,SNP.CV=NULL,SNP.robust="GLM",
 file.from=1, file.to=1, file.total=NULL, file.fragment = 99999,file.path=NULL, 
 file.G=NULL, file.Ext.G=NULL,file.GD=NULL, file.GM=NULL, file.Ext.GD=NULL,file.Ext.GM=NULL, 
 ngrid = 100, llim = -10, ulim = 10, esp = 1e-10,LD.chromosome=NULL,LD.location=NULL,LD.range=NULL,PCA.col=NULL,PCA.3d=FALSE,NJtree.group=NULL,NJtree.type=c("fan","unrooted"),
 sangwich.top=NULL,sangwich.bottom=NULL,QC=TRUE,GTindex=NULL,LD=0.1,plot.bin=10^5,
 file.output=TRUE,cutOff=0.05, Model.selection = FALSE,output.numerical = FALSE,
 output.hapmap = FALSE, Create.indicator = FALSE,Multi_iter=FALSE,num_regwas=10,opt="extBIC",
  QTN=NULL, QTN.round=1,QTN.limit=0, QTN.update=TRUE, QTN.method="Penalty", Major.allele.zero = FALSE,Random.model=FALSE,
  method.GLM="FarmCPU.LM",method.sub="reward",method.sub.final="reward",method.bin="static",bin.size=c(1000000),bin.selection=c(10,20,50,100,200,500,1000),
  memo=NULL,Prior=NULL,ncpus=1,maxLoop=3,threshold.output=.01,Inter.Plot=FALSE,Inter.type=c("m","q"),
  WS=c(1e0,1e3,1e4,1e5,1e6,1e7),alpha=c(.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),maxOut=100,QTN.position=NULL,CG=NULL,
  converge=1,iteration.output=FALSE,acceleration=0,iteration.method="accum",PCA.View.output=TRUE,Geno.View.output=TRUE,plot.style="Oceanic",SUPER_GD=NULL,SUPER_GS=FALSE,
		    h2=NULL,NQTN=NULL,QTNDist="normal",effectunit=1,category=1,r=0.25,cveff=NULL,a2=0,adim=2,Multiple_analysis=FALSE,
  model="MLM",Para=NULL
		){
#Object: To perform GWAS and GPS (Genomic Prediction/Selection)
#Designed by Zhiwu Zhang
#Writen by Jiabo Wang
#Last update: Novenber 3, 2016
##############################################################################################
print("--------------------- Welcome to GAPIT ----------------------------")
echo=TRUE
all.memo=NULL

GAPIT.Version=GAPIT.0000()
#if(!is.null(model))if(!match(model,c("MLM","CMLM","SUPER","GLM","FarmCPU","Blink","BlinkC","MLMM","gBLUP","cBLUP","sBLUP"))) stop(paste("PLease choose one model from ","MLM","CMLM","SUPER","GLM","FarmCPU","Blink","gBLUP","cBLUP","sBLUP",sep=""))
#Allow either KI or K, but not both
if(model%in%c("gBLUP","cBLUP","sBLUP"))
  {
    SNP.test=FALSE
    SUPER_GS=TRUE
  }
if(!is.null(KI)&is.null(GD)&is.null(G)&is.null(file.G)&is.null(file.GD)) SNP.test=FALSE
model_store=model
KI0=KI


print(model_store)
if(!is.null(Y))
  {
     for(m in 1:length(model_store))
        {
        model=model_store[m]
        if(toupper(model)=="BLINK") model="Blink"
        if(toupper(model)=="FARMCPU") model="FarmCPU"
        if(toupper(model)=="BLINKC") model="BlinkC"
        if(toupper(model)=="GBLUP") model="gBLUP"
        if(toupper(model)=="CBLUP") model="cBLUP"
        if(toupper(model)=="SBLUP") model="sBLUP"
        if(toupper(model)=="FARMCPU2") 
        {model="FarmCPU2"
         Multi_iter=TRUE
         # memo=paste(memo,"_Back",sep="")
        }
        if(toupper(model)=="BLINK2") 
        {model="Blink2"
         Multi_iter=TRUE
         # memo=paste(memo,"_Back",sep="")
        }
        if(toupper(model)=="MLMM2") 
        {model="MLMM2"
         Multi_iter=TRUE
         # memo=paste(memo,"_Back",sep="")
        }

        if(group.from<nrow(Y)) model="CMLM"
  # }  
        if(group.to!=group.from)model="CMLM"
        if(group.to==1&group.from==1)model="GLM"
        if(!is.null(sangwich.bottom)&!is.null(sangwich.bottom))model="SUPER"
        if(model=="gBLUP") model="MLM"
        if(model=="cBLUP") model="CMLM"
        if(model=="sBLUP") 
          { model="SUPER"
            Para$group.from=1000000
            Para$group.to=1000000
            Para$group.by=nrow(Y)/10
          }
#CMLM
        if(model=="GLM")
          {
            Para$group.from=1
            Para$group.to=1
            Para$group.by=group.by
            Para$kinship.algorithm="VanRaden"
          }
        if(model=="MLM")
          {
            Para$group.from=1000000
            Para$group.to=1000000
            Para$group.by=group.by
            Para$kinship.algorithm="VanRaden"
          }
        if(model=="CMLM")
          {
            if(group.from>=group.to)Para$group.from=1
            Para$group.to=group.to
            Para$group.by=group.by
            Para$kinship.algorithm="VanRaden"
#if(Para$group.from==Para$group.to)Para$group.from=10
            print(group.from)
            print(group.to)
          }
        if(model=="SUPER")
          {
            if(!is.null(inclosure.from)&is.null(Para$inclosure.from))Para$inclosure.from=inclosure.from
            if(is.null(Para$inclosure.from))Para$inclosure.from=10
            if(!is.null(inclosure.to)&is.null(Para$inclosure.to))Para$inclosure.to=inclosure.to
            if(is.null(Para$inclosure.to))Para$inclosure.to=100
            if(!is.null(inclosure.by)&is.null(Para$inclosure.by))Para$inclosure.by=inclosure.by
            if(is.null(Para$inclosure.by))Para$inclosure.by=10
            if(!is.null(bin.from)&is.null(Para$bin.from))Para$bin.from=bin.from  
            if(is.null(Para$bin.from))Para$bin.from=10000
            if(!is.null(bin.to)&is.null(Para$bin.to))Para$bin.to=bin.to  
            if(is.null(Para$bin.to))Para$bin.to=10000
            if(!is.null(bin.by)&is.null(Para$bin.by))Para$bin.by=bin.by  
            if(is.null(Para$bin.by))Para$bin.by=10000
            if(!is.null(sangwich.top)&is.null(Para$sangwich.top))Para$sangwich.top=sangwich.top  
            if(is.null(Para$sangwich.top))Para$sangwich.top="MLM"
            if(!is.null(sangwich.bottom)&is.null(Para$sangwich.bottom))Para$sangwich.bottom=sangwich.bottom  
            if(is.null(Para$sangwich.bottom))Para$sangwich.bottom="SUPER"
            Para$kinship.algorithm="VanRaden"
          }
        if(model=="FarmCPU")Para$kinship.algorithm="FarmCPU"
        if(model=="MLMM")Para$kinship.algorithm="MLMM"
        if(model=="Blink")Para$kinship.algorithm="Blink"
        if(model=="FarmCPU2")
        {Para$kinship.algorithm="FarmCPU"
         Para$Multi_iter=TRUE}
        if(model=="MLMM2")
        {Para$kinship.algorithm="MLMM"
        Para$Multi_iter=TRUE}
        if(model=="Blink2")
        {Para$kinship.algorithm="Blink"
        Para$Multi_iter=TRUE}
        if(model=="BlinkC")Para$kinship.algorithm="BlinkC"
        if(is.null(memo))
            {
                Para$memo=model
            }else{
                # print(memo)
                # print(model)
                Para$memo=paste(memo,".",model,sep="")
            }
        all.memo=c(all.memo,Para$memo)
# print(Para$memo)
GAPIT_list=list(group.from=group.from ,group.to=group.to,group.by=group.by,DPP=DPP,kinship.cluster=kinship.cluster, kinship.group=kinship.group,kinship.algorithm=kinship.algorithm, FDRcut=FDRcut,
         bin.from=bin.from,bin.to=bin.to,bin.by=bin.by,inclosure.from=inclosure.from,inclosure.to=inclosure.to,inclosure.by=inclosure.by,SNP.P3D=SNP.P3D,SNP.effect=SNP.effect,SNP.impute=SNP.impute,PCA.total=PCA.total, SNP.fraction = SNP.fraction, seed = seed, BINS = 20,SNP.test=SNP.test,
         SNP.MAF=SNP.MAF,FDR.Rate = FDR.Rate, SNP.FDR=SNP.FDR,SNP.permutation=SNP.permutation,SNP.CV=NULL,SNP.robust="GLM",file.from=file.from, file.to=file.to, file.total=file.total, file.fragment = file.fragment,file.path=file.path, 
         file.G=file.G, file.Ext.G=file.Ext.G,file.GD=file.GD, file.GM=file.GM, file.Ext.GD=file.Ext.GD,file.Ext.GM=file.Ext.GM,ngrid = 100, llim = -10, ulim = 10, esp = 1e-10,Inter.Plot=Inter.Plot,Inter.type=Inter.type,
         LD.chromosome=LD.chromosome,LD.location=LD.location,LD.range=LD.range,PCA.col=PCA.col,PCA.3d=PCA.3d,NJtree.group=NJtree.group,NJtree.type=NJtree.type,opt=opt,
         sangwich.top=sangwich.top,sangwich.bottom=sangwich.bottom,QC=QC,GTindex=GTindex,LD=LD,plot.bin=plot.bin,file.output=file.output,cutOff=cutOff, Model.selection = Model.selection,output.numerical = output.numerical,
         output.hapmap = output.hapmap, Create.indicator = Create.indicator,QTN=QTN, QTN.round=1,QTN.limit=0, QTN.update=TRUE, QTN.method="Penalty", Major.allele.zero = Major.allele.zero,
         method.GLM=method.GLM,method.sub=method.sub,method.sub.final="reward",method.bin="static",bin.size=bin.size,bin.selection=bin.selection,model=model,Random.model=Random.model,
         h2=h2,NQTN=NQTN,QTNDist="normal",effectunit=effectunit,category=category,r=r,cveff=NULL,a2=0,adim=2,Multi_iter=Multi_iter,num_regwas=num_regwas,
         memo="",Prior=NULL,ncpus=1,maxLoop=maxLoop,threshold.output=threshold.output,WS=c(1e0,1e3,1e4,1e5,1e6,1e7),alpha=alpha,maxOut=100,QTN.position=QTN.position,CG=CG,
         converge=converge,iteration.output=iteration.output,acceleration=0,iteration.method="accum",PCA.View.output=PCA.View.output,Geno.View.output=Geno.View.output,plot.style="Oceanic",SUPER_GD=NULL,SUPER_GS=SUPER_GS,Multiple_analysis=Multiple_analysis)
        
        G_list_M=rownames(as.matrix(GAPIT_list))
        P_list_M=rownames(as.matrix(Para))

        Para=c(GAPIT_list[!G_list_M%in%P_list_M],Para)
#print(Para$kinship.algorithm)
        if(SUPER_GS==TRUE)Para$SNP.test=FALSE
        IC=NULL
#GAPIT.Version=GAPIT.0000()
        print("--------------------Processing traits----------------------------------")
        # if(!is.null(Y)){
        print("Phenotype provided!")
        if(ncol(Y)<2)  stop ("Phenotype should have taxa name and one trait at least. Please correct phenotype file!")
        print(paste("The ",m," model in all.",sep=""))
        print(model)
        if(m==1)
          {
            DP=GAPIT.DP(G=G,GD=GD,GM=GM,KI=KI0,Z=Z,CV=CV,CV.Inheritance=Para$CV.Inheritance,GP=GP,GK=GK,
            group.from=Para$group.from ,group.to= Para$group.to,group.by=Para$group.by,DPP= Para$DPP, FDRcut=Para$FDRcut,
            kinship.cluster=Para$kinship.cluster, kinship.group=Para$kinship.group,kinship.algorithm=Para$ kinship.algorithm, NJtree.group=Para$NJtree.group,NJtree.type=Para$NJtree.type,plot.bin=Para$plot.bin,PCA.col=Para$PCA.col,PCA.3d=Para$PCA.3d,
             sangwich.top=Para$sangwich.top,sangwich.bottom=Para$sangwich.bottom,LD=Para$LD,bin.from= Para$bin.from,bin.to= Para$bin.to,bin.by= Para$bin.by,inclosure.from= Para$inclosure.from,inclosure.to= Para$inclosure.to,inclosure.by= Para$inclosure.by,
             SNP.P3D= Para$SNP.P3D,SNP.effect= Para$SNP.effect,SNP.impute= Para$SNP.impute,PCA.total= Para$PCA.total, SNP.fraction = Para$SNP.fraction, seed = Para$seed, 
             BINS = Para$BINS,SNP.test=Para$SNP.test, SNP.MAF= Para$SNP.MAF,FDR.Rate = Para$FDR.Rate, SNP.FDR= Para$SNP.FDR,SNP.permutation= Para$SNP.permutation,opt=Para$opt,
             SNP.CV= Para$SNP.CV,SNP.robust= Para$SNP.robust,   Inter.Plot=Para$Inter.Plot,  Inter.type=Para$Inter.type,   
             file.from= Para$file.from, file.to=Para$file.to, file.total= Para$file.total, file.fragment = Para$file.fragment,file.path= Para$file.path, 
             file.G= Para$file.G, file.Ext.G= Para$file.Ext.G,file.GD= Para$file.GD, file.GM= Para$file.GM, file.Ext.GD= Para$file.Ext.GD,file.Ext.GM= Para$file.Ext.GM, 
             ngrid = Para$ngrid, llim = Para$llim, ulim = Para$ulim, esp = Para$esp,Multi_iter=Para$Multi_iter,num_regwas=Para$num_regwas,
             LD.chromosome= Para$LD.chromosome,LD.location= Para$LD.location,LD.range= Para$LD.range,
             QC= Para$QC,GTindex= Para$GTindex,cutOff=Para$cutOff, Model.selection = Para$Model.selection,output.numerical = Para$output.numerical,Random.model=Para$Random.model,
             Create.indicator = Para$Create.indicator,QTN= Para$QTN, QTN.round= Para$QTN.round,QTN.limit= Para$QTN.limit, QTN.update= Para$QTN.update, QTN.method= Para$QTN.method, Major.allele.zero = Para$Major.allele.zero,
             method.GLM=Para$method.GLM,method.sub= Para$method.sub,method.sub.final= Para$method.sub.final,
             method.bin= Para$method.bin,bin.size= Para$bin.size,bin.selection= Para$bin.selection,
             memo= Para$memo,Prior= Para$Prior,ncpus=Para$ncpus,maxLoop= Para$maxLoop,threshold.output= Para$threshold.output,
             WS= Para$WS,alpha= Para$alpha,maxOut= Para$maxOut,QTN.position= Para$QTN.position, converge=Para$converge,iteration.output= Para$iteration.output,acceleration=Para$acceleration,
             iteration.method= Para$iteration.method,PCA.View.output=Para$PCA.View.output, 
             output.hapmap = Para$output.hapmap, file.output= Para$file.output,Geno.View.output=Para$Geno.View.output,plot.style=Para$plot.style,SUPER_GD= Para$SUPER_GD,SUPER_GS= Para$SUPER_GS,CG=Para$CG,model=model)
          }else{ 
             DP$kinship.algorithm=Para$kinship.algorithm
             DP$group.from=Para$group.from
             DP$group.to=Para$group.to
             DP$group.by=Para$group.by
             DP$sangwich.top=Para$sangwich.top
             DP$sangwich.bottom=Para$sangwich.bottom
             DP$bin.from= Para$bin.from
             DP$bin.to= Para$bin.to
             DP$bin.by= Para$bin.by
             DP$inclosure.from = Para$inclosure.from
             DP$inclosure.to= Para$inclosure.toDP$inclosure.by= Para$inclosure.by
             DP$Multi_iter=Para$Multi_iter
          }

        for (trait in 2: ncol(Y))  
          {
             traitname=colnames(Y)[trait]
###Statistical distributions of phenotype
###Correlation between phenotype and principal components
             print(paste("Processing trait: ",traitname,sep=""))
             if(!is.null(Para$memo)) traitname=paste(Para$memo,".",traitname,sep="")
             if(!is.null(Y) & Para$file.output)ViewPhenotype<-GAPIT.Phenotype.View(myY=Y[,c(1,trait)],traitname=traitname,memo=Para$memo)
             if(!Para$kinship.algorithm%in%c("FarmCPU","MLMM","Blink","BlinkC")&is.null(DP$KI))
             {
                myKI_test=GAPIT.kinship.VanRaden(snps=as.matrix(DP$GD[,-1]))     #  build kinship
                colnames(myKI_test)=as.character(DP$GD[,1])
                KI0=cbind(as.character(DP$GD[,1]),as.data.frame(myKI_test))
             }
             if(!is.null(KI0))DP$KI=KI0 
             Judge=GAPIT.Judge(Y=Y[,c(1,trait)],G=DP$G,GD=DP$GD,KI=DP$KI,GM=DP$GM,group.to=DP$group.to,group.from=DP$group.from,sangwich.top=DP$sangwich.top,sangwich.bottom=DP$sangwich.bottom,kinship.algorithm=DP$kinship.algorithm,PCA.total=DP$PCA.total,model=DP$model,SNP.test=DP$SNP.test)
             DP$group.from=Judge$group.from
             DP$group.to=Judge$group.to
             DP$name.of.trait=traitname
             DP$Y=Y[!is.na(Y[,trait]),c(1,trait)]
             DP$model=model
# print(Para$SNP.test)
             IC=GAPIT.IC(DP=DP)
             SS=GAPIT.SS(DP=DP, IC=IC, buspred=buspred, lmpred=lmpred)
             if(Para$SNP.test&Para$file.output)ID=GAPIT.ID(DP=DP,IC=IC,SS=SS)
          }#for loop trait
#print(SNP.test)
        print("GAPIT accomplished successfully for multiple traits. Result are saved")
        print("It is OK to see this: 'There were 50 or more warnings (use warnings() to see the first 50)'")
        out <- list()
        out$QTN<-QTN.position
        out$GWAS<-SS$GWAS
        out$Pred<-SS$Pred
        out$QTN<-IC$QTN
        out$Power<-SS$Power
        out$FDR<-SS$FDR
        out$Power.Alpha<-SS$Power.Alpha
        out$alpha<-SS$alpha
        out$mc=SS$mc
        out$bc=SS$bc
        out$mp=SS$mp
        out$h2=SS$h2
        out$PCA=IC$myallCV
        out$GD=DP$GD
        out$GM=DP$GM
        out$KI=IC$K
        out$GM=DP$GM
        out$Compression=SS$Compression
        if(Para$SNP.test==TRUE)names(out$GWAS$P.value)="mp"
        if(kinship.algorithm=="FarmCPU")names(out$Pred)=c("Taxa",traitname,"Prediction")
#return (out)
        }#end of model loop
  }else{# is.null(Y)
  #print(Para$SNP.MAF)
        out <- list()
        GAPIT_list=list(group.from=group.from ,group.to=group.to,group.by=group.by,DPP=DPP,kinship.cluster=kinship.cluster, kinship.group=kinship.group,kinship.algorithm=kinship.algorithm, 
         bin.from=bin.from,bin.to=bin.to,bin.by=bin.by,inclosure.from=inclosure.from,inclosure.to=inclosure.to,inclosure.by=inclosure.by,SNP.P3D=SNP.P3D,SNP.effect=SNP.effect,SNP.impute=SNP.impute,PCA.total=PCA.total, SNP.fraction = SNP.fraction, seed = seed, BINS = 20,SNP.test=SNP.test,
         SNP.MAF=SNP.MAF,FDR.Rate = FDR.Rate, SNP.FDR=SNP.FDR,SNP.permutation=SNP.permutation,SNP.CV=NULL,SNP.robust="GLM",file.from=file.from, file.to=file.to, file.total=file.total, file.fragment = file.fragment,file.path=file.path, 
         file.G=file.G, file.Ext.G=file.Ext.G,file.GD=file.GD, file.GM=file.GM, file.Ext.GD=file.Ext.GD,file.Ext.GM=file.Ext.GM,ngrid = 100, llim = -10, ulim = 10, esp = 1e-10,Inter.Plot=Inter.Plot,Inter.type=Inter.type,
         LD.chromosome=LD.chromosome,LD.location=LD.location,LD.range=LD.range,PCA.col=PCA.col,PCA.3d=PCA.3d,NJtree.group=NJtree.group,NJtree.type=NJtree.type,opt=opt,
         sangwich.top=sangwich.top,sangwich.bottom=sangwich.bottom,QC=QC,GTindex=GTindex,LD=LD,plot.bin=plot.bin,file.output=file.output,cutOff=cutOff, Model.selection = Model.selection,output.numerical = output.numerical,
         output.hapmap = output.hapmap, Create.indicator = Create.indicator,QTN=QTN, QTN.round=1,QTN.limit=0, QTN.update=TRUE, QTN.method="Penalty", Major.allele.zero = Major.allele.zero,
         method.GLM=method.GLM,method.sub=method.sub,method.sub.final="reward",method.bin="static",bin.size=bin.size,bin.selection=bin.selection,model=model,Random.model=Random.model,
         h2=h2,NQTN=NQTN,QTNDist="normal",effectunit=effectunit,category=category,r=r,cveff=NULL,a2=0,adim=2,Multi_iter=Multi_iter,num_regwas=num_regwas,
         memo="",Prior=NULL,ncpus=1,maxLoop=maxLoop,threshold.output=threshold.output,WS=c(1e0,1e3,1e4,1e5,1e6,1e7),alpha=alpha,maxOut=100,QTN.position=QTN.position,CG=CG,
         converge=converge,iteration.output=iteration.output,acceleration=0,iteration.method="accum",PCA.View.output=PCA.View.output,Geno.View.output=Geno.View.output,plot.style="Oceanic",SUPER_GD=NULL,SUPER_GS=SUPER_GS,Multiple_analysis=Multiple_analysis)
        if(model=="MLM")
          {
            Para$group.from=1000000
            Para$group.to=1000000
            Para$group.by=group.by
          }
        G_list_M=rownames(as.matrix(GAPIT_list))
        P_list_M=rownames(as.matrix(Para))
        if(is.null(memo))
            {
                Para$memo=model
            }else{
                Para$memo=paste(memo,".",mode,sep="")
            }
        all.memo=c(all.memo,Para$memo)
        Para=c(GAPIT_list[!G_list_M%in%P_list_M],Para)
        myGenotype<-GAPIT.Genotype(G=G,GD=GD,GM=GM,KI=KI,kinship.algorithm=kinship.algorithm,PCA.total=PCA.total,SNP.fraction=SNP.fraction,SNP.test=SNP.test,
 file.path=file.path,file.from=file.from, file.to=file.to, file.total=file.total, file.fragment = file.fragment, file.G=file.G, 
 file.Ext.G=file.Ext.G,file.GD=file.GD, file.GM=file.GM, file.Ext.GD=file.Ext.GD,file.Ext.GM=file.Ext.GM,
 SNP.MAF=SNP.MAF,FDR.Rate = FDR.Rate,SNP.FDR=SNP.FDR,SNP.effect=SNP.effect,SNP.impute=SNP.impute,NJtree.group=NJtree.group,NJtree.type=NJtree.type,
 LD.chromosome=LD.chromosome,LD.location=LD.location,LD.range=LD.range,GP=GP,GK=GK,bin.size=NULL,inclosure.size=NULL, 
 sangwich.top=NULL,sangwich.bottom=sangwich.bottom,GTindex=NULL,file.output=file.output, Create.indicator = Create.indicator, Major.allele.zero = Major.allele.zero,Geno.View.output=Geno.View.output,PCA.col=PCA.col,PCA.3d=PCA.3d)
        GD=myGenotype$GD
        GI=myGenotype$GI
        GT=myGenotype$GT
#G=myGenotype$G
        chor_taxa=myGenotype$chor_taxa
        rownames(GD)=GT
        colnames(GD)=GI[,1]
        taxa=GT
   if(!is.null(chor_taxa))
   {
     chro=as.numeric(as.matrix(GI[,2]))
     for(i in 1:length(chro))
     {
      chro[chro==i]=chor_taxa[i]
     }
     GI[,2]=chro
   }
#print(GD[1:5,1:5])
        if(output.numerical) 
          {
            utils::write.table(cbind(taxa,GD),  "GAPIT.Genotype.Numerical.txt", quote = FALSE, sep = "\t", row.names = F,col.names = T)
            utils::write.table(GI,  "GAPIT.Genotype.map.txt", quote = FALSE, sep = "\t", row.names = F,col.names = T)
          }
        if(output.hapmap) utils::write.table(myGenotype$G,  "GAPIT.Genotype.hmp.txt", quote = FALSE, sep = "\t", row.names = FALSE,col.names = FALSE)
#GD=cbind(as.data.frame(GT),GD)
        if(!is.null(seed))set.seed(seed)
#print(Para$NQTN)
        if(!is.null(Para$NQTN)&!is.null(Para$h2))
          {
            myG_simulation<-GAPIT.Phenotype.Simulation(GD=cbind(as.data.frame(myGenotype$GT),myGenotype$GD),GM=myGenotype$GI,h2=Para$h2,NQTN=Para$NQTN,QTNDist=Para$QTNDist,effectunit=Para$effectunit,category=Para$category,r=Para$r,cveff=Para$cveff,a2=Para$a2,adim=Para$adim)
            out=c(out,myG_simulation)
          }
        out$GD=data.frame(cbind(as.data.frame(GT),as.data.frame(GD)))
        out$GM=GI
        out$G=myGenotype$G
        out$kinship=myGenotype$KI
        out$PCA=myGenotype$PC
        out$chor_taxa=chor_taxa
  }# is.null(Y)

#print(tail(IC$GM))
model_store=all.memo
if(!is.null(Y)&SNP.test)if(Multiple_analysis&Para$file.output&length(model_store)*(ncol(Y)-1)>1&length(model_store)*(ncol(Y)-1)<9)
  { 
  #print(DP$QTN.position)
   GMM=GAPIT.Multiple.Manhattan(model_store=model_store,Y=Y,GM=IC$GM,seqQTN=DP$QTN.position,cutOff=DP$cutOff)
#print(str(GMM$multip_mapP))
   GAPIT.Circle.Manhatton.Plot(band=1,r=3,GMM$multip_mapP,plot.type=c("c","q"),signal.line=1,xz=GMM$xz,threshold=DP$cutOff)
  }# end of mutiple manhantton plot

if(file.output&!SNP.test&model_store%in%c("gBLUP","cBLUP","sBLUP")&Inter.Plot)
  { 
  print("here will start interactive for GS !!!")
  GAPIT.Interactive.GS(model_store=model_store,Y=Y)
  if(!is.null(testY))GAPIT.Interactive.GS(model_store=model_store,Y=Y,testY=testY)

#print(str(GMM$multip_mapP))
#GAPIT.Circle.Manhatton.Plot(band=1,r=3,GMM$multip_mapP,plot.type=c("c","q"),signal.line=1,xz=GMM$xz,threshold=DP$cutOff)
  }# end of mutiple manhantton plot





return (out)
}  #end of GAPIT function
