
#' GAPIT Genome Association and Prediction Integrated Tools
#' 
#' @description 
#' GWAS and GS procedure using the Multiple models (General Linear Model, Mixed Linear Model, Compression Mixed Linear Model, SUPER, Multiple Loci Mixed linear Model, FarmCPU, and BLINK)
#' 
#' 
#' @param Y  data.frame of phenotype data where each row is a sample and each column is a trait, the first column is the sample names
#' @param G  data.frame of genotypic data in HAPMAP format
#' @param GD data.frame of genetic data in numerical format, where each row is a sample and each column is a variant.
#' @param GM a data.frame of genomic coordinates for the genetic map
#' @param KI an $NxN$ matrix of kinship coefficients
#' @param Z  an $NxN$ (for MLM) or an $NxN`$ (CMLM) matrix of index, which is made with 0 and 1 value to indicate indivdual belong to each group.
#' @param CV Covariance matrix
#' @param testY data.frame of phenotype data in testing population, where each row is a sample and each column is a trait, the first column is the sample names.
#' @param group.from integer, minimum number of group(s) to consider in CMLM
#' @param group.to integer, maximum number of group(s) to consider in CMLM
#' @param group.by integer, increment for evaluating group size in CMLM
#' @param kinship.cluster algorithm for calculating kinship centroid (options: "average", "complete", "ward", "single", "mcquitty", "median", and "centroid") 
#' @param kinship.group method for calculating group membership (options: "Mean", "Max", "Min", and "Median")
#' @param kinship.algorithm algorithm to calculate the kinship matrix (options: "VanRaden", "EMMA", "Loiselle", and "Zhang")
#' @param buspred logical, option for prediction after GWAS.
#' @param lmpred logical (vector), option for seletion of linear model prediction or (and) ABLUP.
#' @param FDRcut logical, filter pseudo QTN based on FDR cut-off in BLINK
#' @param bin.from integer, minimum number of bin(s) to consider in SUPER
#' @param bin.to integer, maximum number of bin(s) to consider in SUPER
#' @param bin.by integer, increment for evaluating bin size in SUPER
#' @param inclosure.from integer, minimum number of pesudo QTNs to consider in SUPER
#' @param inclosure.to integer, maximum number of pesudo QTNs to consider in SUPER
#' @param inclosure.by integer, increment for evaluating number of pesudo QTNs in SUPER
#' @param SNP.P3D logical, to use P3D or Not for Testing SNPs
#' @param SNP.effect genetic model for coding the SNP effect (options: "Add" (additive), "Dom", "Left", and "Right")
#' @param SNP.impute SNP imputation method (options: "Middle", "Major", and "Minor")
#' @param PCA.total integer, number of principal components to include in Q matrix (can be zero)
#' @param SNP.fraction numerical input between 0 and 1, fraction of SNPs Sampled to Estimate Kinship and PCs
#' @param SNP.MAF numerical input between 0 and 1, minor allele frequency to filter SNPs in GWAS reports
#' @param SNP.FDR numerical input between 0 and 1, false discovery rate for filtering SNPs
#' @param PCA.col list for points color in PCA plot. The total length of PCA.col should be equal to the number of individuals in the GD or G file.
#' @param PCA.3d logical, whether output 3D PCA plot.
#' @param NJtree.group numeric, set the number of clustering groups in the NJtree plot.
#' @param NJtree.type type of neighbor joining tree (options: "fan" and "unrooted")
#' @param sangwich.top Model type to run in the first iteration of SUPER, (options: "MLM", "GLM", "CMLM","Fast-LMM")
#' @param sangwich.bottom Model type to run in the last iteration of SUPER, (options: "MLM", "GLM", "CMLM","Fast-LMM")
#' @param file.output logical, whether output all result files.
#' @param cutOff numeric value, the threshold for filtering significant markers from all. It would be transfor as Bornferrni cutoff in Manhattan plots.
#' @param Model.selection logical, whether evaluate optimum number of CV file. If TRUE, all likelyhood values should be evaluated for each CV combination.
#' @param output.numerical logical,whether output numerical genotype file from HapMap file.
#' @param output.hapmap logical,whether output numerical HapMap file from numerical genotype file.
#' @param Multi_iter logical, whether add more iterations for FarmCPU, BLINK.
#' @param num_regwas numeric, the maximum number of selective significant markers into re-GWAS model.
#' @param Major.allele.zero logical, whether set major allele as 0, and minor allele as 2, if FALSE, they will be set as reverse.
#' @param Random.model logical, whether ran random model to estimate PVE values for significant markers after GWAS.
#' @param memo text, from users to remark for output files.
#' @param Inter.Plot logical, whether to output the interactive Manhattan and QQ plots.
#' @param Inter.type Interactive plot type for Manhattan and QQ plots."m" indicate manhattan plot and "q" indicate QQ plot.
#' @param WS numeric or numeric vector, the distance between detected markers and real QTN should be recognized as a real power.
#' @param WS0 numeric, the cutoff threshold for distance between markers to display in GAPIT.Genotype.Distance_R_Chro.pdf file.
#' @param Aver.Dis numeric, average display windowsize in LD decay plot,
#' @param maxOut numeric, set the number of markers in the power calculation, the top maxOut number of P values markers should be selected.
#' @param QTN.position numeric vector, set where are the QTNs' position. Its maximun values should be equal to total marker number, and its length should be equal to the NQTN.
#' @param PCA.View.output logical, whether to output the PCA view
#' @param Geno.View.output logical whether to output the Genotype analysis including MAF, heterzygosity, LD decay, and other genotype distribution output.
#' @param h2 numeric value, to set simulation phenotype heritability. It ranged from 0 to 1 means 0\% to 100\%.
#' @param NQTN numeric value, to set simulation number of QTN. It ranged from 1 to the total markers number.
#' @param QTNDist option for distribution of simulated QTN genetic effect in the simulation,(options: "normal" and "geometry")
#' @param effectunit numeric value, the effect unit of the first choosed marker in the simulation pheotype. default as 1
#' @param Multiple_analysis logical, whether to output the mulitple mahattan and QQ plots. default as TRUE 
#' @param model model type to run, (options: "MLM", "GLM", "CMLM", "MMLM", "SUPER", "FarmCPU", "gBLUP",  "cBLUP", and "sBLUP"
#' @param Predict.type option to display which type predicted factor again real phenotype in the GAPIT.Association.Prediction pdf file.(options: "GEBV","BLUP" and "BLUE")
#' @param SNP.test logical, whether to do GWAS or GS.
#' @param seq.cutoff numeric value, the threshold for filtering significant markers from all. It would be transfor as Bornferrni cutoff in GGS.
#'
#' @param CV.Extragenetic param
#' @param bin.regwas param
#' @param N4 param
#' @param N.sig param
#' @param QC.Y param
#' @param QTN.gs param
#' @param r param
#' @param seq.num param
#' @param plot.bin param
#' @param PCA.legend param
#' @param Phenotype.View param
#'


#' @details 
#' Genome Association and Prediction Integrated Tools
#' Available models: MLM, GLM, CMLM, MMLM, SUPER, FarmCPU, gBLUP, cBLUP
#' 
#' 
#' @return 
#' A list
#' including some of the following elements:MLM, GLM, CMLM, MMLM, SUPER, FarmCPU, gBLUP, cBLUP
#'
#'
#' @seealso 
#' GAPIT.DP(), GAPIT.Phenotype.View(), GAPIT.judge(), GAPIT.IC(), GAPIT.SS(), GAPIT.ID().
#' 
#' 
#' library(help = "GAPIT")
#' 
#' @author Zhiwu Zhang and Jiabo Wang
#' 
#' @examples 
#' \dontrun{
#' 
#' myPhenoFile <- system.file("extdata", "mdp_traits.txt.gz", package = "GAPIT")
#' myGenoFile <- system.file("extdata", "mdp_genotype_test.hmp.txt.gz", package = "GAPIT")
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
`GAPIT` <- function(
  Y = NULL, #phenotype
  G = NULL, #hapmap genotype
  GD = NULL, #numeric genotype
  GM = NULL, #genotype map information
  KI = NULL, #kinship
  Z = NULL, #Z matrix for MLM, cMLM, encMLM
  CV = NULL, #corvariance matrix
  Aver.Dis=1000,
  # a2 = 0,
  # adim = 2,
  # acceleration = 0,
  # alpha = c(.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1), # confidence coefficient
  buspred = FALSE, #Bus prediction
  bin.from = 10000, #SUPER 
  bin.to = 10000, #SUPER
  bin.by = 10000, #SUPER
  # bin.size = c(1000000), 
  # bin.selection = c(10,20,50,100,200,500,1000),
  # BINS = 20,
  # converge = 1,
  cutOff = 0.05, #threshold for significant
  # category = 1, #Simulation phenotype
  # cveff = NULL, #Simulation phenotype
  # Create.indicator = FALSE, #
  # CG = NULL, #candidate gene matrix for relationship
  CV.Extragenetic = 0, # the top number of no-inheritance columns in CV
  # Cross.Vali=TRUE,
  # color0=NULL,
  # DPP = 100000, #content points in Manhattan Plot
  # DP=NULL,
  # esp = 1e-10,
  effectunit = 1, #Simulation phenotype
  # file.from = 1,  #read seqarated data files
  # file.to = 1, #read seqarated data files
  # file.total = NULL, #read seqarated data files
  # file.fragment = 99999,#read seqarated data files
  # file.path = NULL, #read seqarated data files
  # file.G = NULL, #read seqarated data files
  # file.Ext.G = NULL,#read seqarated data files
  # file.GD = NULL, #read seqarated data files
  # file.GM = NULL, #read seqarated data files
  # file.Ext.GD = NULL,#read seqarated data files
  # file.Ext.GM = NULL, #read seqarated data files
  file.output = TRUE, #output option
  # FDR.Rate = 1, # filter FDR
  FDRcut = FALSE, # filter pseudo QTN based on cutOff in blink
  group.from = 1000000,#MLM
  group.to = 1000000,#MLM
  group.by = 50,#MLM
  # GTindex = NULL,
  Geno.View.output = TRUE,#genotype analysis option
  # GP = NULL,
  # GK = NULL, #group kinship
  h2 = NULL, #simulation phenotype heritability
  inclosure.from = 10, #SUPER
  inclosure.to = 10, #SUPER
  inclosure.by = 10, #SUPER
  # iteration.output = FALSE,
  # iteration.method = "accum",
  # inpch=NULL, # in pch of S manhattans
  Inter.Plot = FALSE, #Interactive plot option
  Inter.type = c("m","q"), #Interactive plot type for Manhattan and QQ plots
  kinship.cluster = "average", #cMLM
  kinship.group = 'Mean',#cMLM
  kinship.algorithm = "Zhang",#cMLM
  # llim = -10, 
  lmpred = FALSE, #option for linear model prediction or ABLUP prediction, that could be set as multiple parameters
  # LD.chromosome = NULL, #LD plot of markers in significant marker region
  # LD.location = NULL, #LD plot of markers in significant marker region
  # LD.range = NULL, #LD plot of markers in significant marker region
  # LD = 0.1, #SUPER
  model = "MLM",# model or method in GWAS or GS
  # method.GLM = "FarmCPU.LM", 
  # method.sub = "reward",
  # method.sub.final = "reward",
  # method.bin = "static",
  maxOut = 100, # power for top number of markers in the GWAS
  memo = NULL, #label for marking
  # maxLoop = 3,
  Model.selection = FALSE,# optimum number of CV and PCAs
  Multi_iter = TRUE, #Multiple step for FarmCPU and BLink
  Major.allele.zero = FALSE, #convert hapmap file to numeric file, set major marker as 0
  Multiple_analysis = TRUE, #option for multiple Manhattan and QQ plots
  num_regwas = 10,# the max number of Multiple markers 
  bin.regwas = 100000,
  # ncpus = 1,
  # ngrid = 100, 
  N4=FALSE,
  NQTN = NULL, #Simulation phenotype, number of QTN
  N.sig=NULL, #Random.model, Number of significant markers
  NJtree.group = NULL, #NJtree set number of cluster group
  NJtree.type = c("fan","unrooted"),#NJtree type
  # opt = "extBIC",
  output.numerical = FALSE,# option for output numeric files
  output.hapmap = FALSE, # option for output hapmap files
  # outpch=NULL, # out pch of S manhattans
  # QTN = NULL, 
  # QTN.round = 1,
  # QTN.limit = 0, 
  # QTN.update = TRUE, 
  # QTN.method = "Penalty", 
  # QC = TRUE,
  QC.Y=FALSE,
  QTN.position = NULL, #Simulation phenotype, QTN position in the order of map file
  QTN.gs = 0, # The number of QTNs in the CV file
  QTNDist = "normal",
  r = 0.25,
  Random.model = TRUE, #Random.model to calculate PVE
  sangwich.top = NULL, #SUPER
  sangwich.bottom = NULL,#SUPER
  seq.cutoff=NULL,
  seq.num=50, # number of selected sig markers into GS
  # seed = NULL, 
  SNP.P3D = TRUE,
  SNP.effect = "Add",
  SNP.impute = "Middle",
  SNP.fraction = 1, 
  SNP.test = TRUE,
  SNP.MAF = 0,
  SNP.FDR = 1,
  # SNP.permutation = FALSE,
  # SNP.CV = NULL,
  # SNP.robust = "GLM",
  # SUPER_GD = NULL,
  # SUPER_GS = FALSE,
  testY = NULL,
  # plot.style = "Oceanic",
  plot.bin = 10^5,
  PCA.total = 0, # PCA number
  PCA.col = NULL, #indicater colors for individuals in PCA plot
  PCA.3d = FALSE, #3D PCA plot option
  PCA.legend=NULL, # PCA legend list
  PCA.View.output = TRUE, #option for PCA plot
  Phenotype.View= TRUE, # option for phenotype view plot
  # Prior = NULL,
  # Para = NULL,
  Predict.type="GEBV",
  # ulim = 10, 
  WS = c(1e0,1e3,1e4,1e5,1e6,1e7),
  WS0 = 10000
	){
#Object: To perform GWAS and GPS (Genomic Prediction/Selection)
#Designed by Zhiwu Zhang
#Writen by Jiabo Wang
#Last update: Mar 8, 2023
##############################################################################################
print("--------------------- Welcome to GAPIT ----------------------------")
all.memo=NULL
GAPIT.Version=GAPIT.0000()
#Allow either KI or K, but not both

if(!is.null(KI)&is.null(GD)&is.null(G)) SNP.test=FALSE
# model_store=model
KI0=KI
model_store=append(model[!model%in%c("gBLUP","cBLUP","sBLUP")],model[model%in%c("gBLUP","cBLUP","sBLUP")])

print(model_store)

if(!is.null(Y))
  {
     for(m in 1:length(model_store))
        {
        # print(model_store)
        model=model_store[m]
        # print(model)
        
        if(toupper(model)=="BLINK") model="BLINK"
        if(toupper(model)=="FARMCPU") model="FarmCPU"
        if(toupper(model)=="BLINKC") model="BLINKC"
        if(toupper(model)=="GBLUP") model="gBLUP"
        if(toupper(model)=="CBLUP") model="cBLUP"
        if(toupper(model)=="SBLUP") model="sBLUP"
        if(toupper(model)=="FARMCPU2") 
        {model="FarmCPU2"
         Multi_iter=TRUE
        }
        if(toupper(model)=="BLINK2") 
        {model="BLINK2"
         Multi_iter=TRUE
        }
        if(toupper(model)=="MLMM2") 
        {model="MLMM2"
         Multi_iter=TRUE
        }
        if(model%in%c("gBLUP","cBLUP","sBLUP"))
        {
          SNP.test=FALSE
          SUPER_GS=TRUE
        }else{
          # SNP.test=TRUE
          SUPER_GS=FALSE
        }
        # if(group.to!=group.from)model="CMLM"
        # if(group.to==1&group.from==1)model="GLM"
        # if(!is.null(sangwich.bottom)&!is.null(sangwich.bottom))model="SUPER"
        if(model=="GLM")
          {
            group.from=1
            group.to=1
            if(is.null(kinship.algorithm))kinship.algorithm="Zhang"
            Random.model=FALSE
          }
        if(model=="MLM"|model=="gBLUP")
          {
            group.from=1000000
            group.to=1000000
            if(is.null(kinship.algorithm))kinship.algorithm="Zhang"
          }
        if(model=="CMLM"|model=="cBLUP")
          {
            if(group.from>=group.to)group.from=1
            print(group.from)
            print(group.to)
            if(is.null(kinship.algorithm))kinship.algorithm="Zhang"
          }
        if(model=="SUPER"|model=="sBLUP")
          {
            if(is.null(inclosure.from))inclosure.from=10
            if(is.null(inclosure.to))inclosure.to=100
            if(is.null(inclosure.by))inclosure.by=10
            if(is.null(bin.from))bin.from=10000
            if(is.null(bin.to))bin.to=10000
            if(is.null(bin.by))bin.by=10000
            if(is.null(sangwich.top))sangwich.top="MLM"
            if(is.null(sangwich.bottom))sangwich.bottom="SUPER"
            if(is.null(kinship.algorithm))kinship.algorithm="Zhang"
            group.from=1000000
            group.to=1000000
            group.by=nrow(Y)/10
          }
        if(model=="FarmCPU")kinship.algorithm="FarmCPU"
        if(model=="MLMM")kinship.algorithm="MLMM"
        if(model=="BLINK")kinship.algorithm="BLINK"
        if(model=="FarmCPU2")
          {
            kinship.algorithm="FarmCPU"
            Multi_iter=TRUE
          }
        if(model=="MLMM2")
          {
            kinship.algorithm="MLMM"
            Multi_iter=TRUE
          }
        if(model=="BLINK2")
          {
            kinship.algorithm="BLINK"
            Multi_iter=TRUE
          }
        if(model=="BLINKC")kinship.algorithm="BLINKC"
        if(is.null(memo))
          {
            memo0=model
          }else{
            memo0=paste(memo,".",model,sep="")
          }
        all.memo=c(all.memo,memo0)
        if(SUPER_GS==TRUE)SNP.test=FALSE
        IC=NULL
#GAPIT.Version=GAPIT.0000()
        print("--------------------Processing traits----------------------------------")
        # if(!is.null(Y)){
        print("Phenotype provided!")
        if(ncol(Y)<2)  stop ("Phenotype should have taxa name and one trait at least. Please correct phenotype file!")
        print(paste("The ",m," model in all.",sep=""))
        print(model)
        # print(dim(KI0))
        if(m==1)
          {
            DP=GAPIT.DP(G=G,GD=GD,GM=GM,KI=KI0,Z=Z,CV=CV,CV.Extragenetic=CV.Extragenetic,
            group.from=group.from ,group.to= group.to,group.by=group.by,FDRcut=FDRcut,Major.allele.zero=Major.allele.zero,
            kinship.cluster=kinship.cluster, kinship.group=kinship.group,kinship.algorithm=kinship.algorithm, NJtree.group=NJtree.group,NJtree.type=NJtree.type,PCA.col=PCA.col,PCA.3d=PCA.3d,
             sangwich.top=sangwich.top,sangwich.bottom=sangwich.bottom,bin.from=bin.from,bin.to=bin.to,bin.by=bin.by,inclosure.from=inclosure.from,inclosure.to=inclosure.to,inclosure.by=inclosure.by,
             SNP.P3D=SNP.P3D,SNP.effect=SNP.effect,SNP.impute=SNP.impute,PCA.total=PCA.total, SNP.fraction =SNP.fraction, seed =NULL, 
             SNP.test=SNP.test, SNP.MAF=SNP.MAF,FDR.Rate =1, SNP.FDR=SNP.FDR,
             Inter.Plot=Inter.Plot,  Inter.type=Inter.type,N.sig=N.sig,seq.num=seq.num,
             Multi_iter=Multi_iter,num_regwas=num_regwas,QTN.gs=QTN.gs,bin.regwas=bin.regwas,
             cutOff=cutOff, Model.selection =Model.selection,output.numerical =output.numerical,Random.model=Random.model,
             PCA.legend=PCA.legend,PCA.View.output=PCA.View.output, 
             WS0=WS0,Aver.Dis=Aver.Dis,memo=memo0,WS=WS,maxOut=maxOut,QTN.position=QTN.position, 
             output.hapmap =output.hapmap, file.output= file.output,Geno.View.output=Geno.View.output,SUPER_GS=SUPER_GS,model=model)
          }else{ 
             DP$kinship.algorithm=kinship.algorithm
             DP$group.from=group.from
             DP$group.to=group.to
             DP$group.by=group.by
             DP$sangwich.top=sangwich.top
             DP$sangwich.bottom=sangwich.bottom
             DP$bin.from=bin.from
             DP$bin.to=bin.to
             DP$bin.by=bin.by
             DP$inclosure.from =inclosure.from
             DP$inclosure.to=inclosure.to
             DP$inclosure.by=inclosure.by
             DP$Multi_iter=Multi_iter
             DP$file.output=file.output
             DP$SNP.test=SNP.test
             DP$model=model
          }

        for (trait in 2: ncol(Y))  
          {
             traitname=colnames(Y)[trait]
             traitname0=colnames(Y)[trait]
###Statistical distributions of phenotype
###Correlation between phenotype and principal components
             print(paste("Processing trait: ",traitname,sep=""))
             if(!is.null(Y) & file.output&Phenotype.View&m==1)ViewPhenotype<-GAPIT.Phenotype.View(myY=Y[,c(1,trait)],traitname=traitname)
             if(!is.null(memo0)) traitname=paste(memo0,".",traitname,sep="")
             # print(DP$kinship.algorithm)
             if(!DP$kinship.algorithm%in%c("FarmCPU","MLMM","BLINK","BLINKC")&is.null(DP$KI)&!is.null(DP$GD))
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
             DP$Y=Judge$Y
             if(QC.Y) DP$Y[,2]=GAPIT.Remove.outliers(DP$Y[,2])
             DP$model=model
# print(DP$name.of.trait)
             IC=GAPIT.IC(DP=DP)
             SS=GAPIT.SS(DP=DP, IC=IC, buspred=buspred, lmpred=lmpred)
             if(SNP.test&DP$file.output)ID=GAPIT.ID(DP=DP,IC=IC,SS=SS,testY=testY)
          }#for loop trait
#print(SNP.test)
        print("GAPIT accomplished successfully for multiple traits. Result are saved")
#        print("It is OK to see this: 'There were 50 or more warnings (use warnings() to see the first 50)'")
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
        if(SNP.test)names(out$GWAS$P.value)="mp"
        # if(kinship.algorithm=="FarmCPU")names(out$Pred)=c("Taxa",traitname,"Prediction")
        kinship.algorithm=NULL
        }#end of model loop
  }else{# is.null(Y)
  #print(Para$SNP.MAF)
        SNP.test=FALSE
        out <- list()
        if(model=="MLM")
          {
            group.from=1000000
            group.to=1000000
          }
        if(is.null(memo))
          {
            memo=model
          }else{
            memo=paste(memo,".",model,sep="")
          }
        all.memo=c(all.memo,memo)
        myGenotype<-GAPIT.Genotype(G=G,GD=GD,GM=GM,KI=KI,kinship.algorithm=kinship.algorithm,PCA.total=PCA.total,SNP.fraction=SNP.fraction,SNP.test=SNP.test,
                          file.from=1, file.to=1, file.total=NULL, file.fragment =9999,file.path=NULL, 
             file.G=NULL, file.Ext.G=NULL,file.GD=NULL, file.GM=NULL, file.Ext.GD=NULL,file.Ext.GM= NULL,WS0=WS0,Aver.Dis=Aver.Dis,
                          SNP.MAF=SNP.MAF,FDR.Rate = 1,SNP.FDR=SNP.FDR,SNP.effect=SNP.effect,SNP.impute=SNP.impute,NJtree.group=NJtree.group,NJtree.type=NJtree.type,
                          GP=NULL,GK=NULL,bin.size=NULL,inclosure.size=NULL, PCA.legend=PCA.legend,
                          sangwich.top=NULL,sangwich.bottom=sangwich.bottom,GTindex=NULL,file.output=file.output, Create.indicator = FALSE,
                          Major.allele.zero = Major.allele.zero,Geno.View.output=Geno.View.output,PCA.col=PCA.col,PCA.3d=PCA.3d)
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
        if(output.numerical) 
          {
            utils::write.table(cbind(taxa,GD),  "GAPIT.Genotype.Numerical.txt", quote = FALSE, sep = "\t", row.names = F,col.names = T)
            utils::write.table(GI,  "GAPIT.Genotype.map.txt", quote = FALSE, sep = "\t", row.names = F,col.names = T)
          }
        if(output.hapmap) utils::write.table(myGenotype$G,  "GAPIT.Genotype.hmp.txt", quote = FALSE, sep = "\t", row.names = FALSE,col.names = FALSE)
        # if(!is.null(seed))set.seed(seed)
        if(!is.null(NQTN)&!is.null(h2))
          {
            myG_simulation<-GAPIT.Phenotype.Simulation(GD=cbind(as.data.frame(myGenotype$GT),myGenotype$GD),GM=myGenotype$GI,h2=h2,NQTN=NQTN,QTNDist=QTNDist,effectunit=effectunit)
            out=c(out,myG_simulation)
            if(file.output)ViewPhenotype<-GAPIT.Phenotype.View(myY=myG_simulation$Y,traitname="Simulated.Phenotype",memo=memo0)
          }
        print("Now the GAPIT is cbind taxa and numeric genotype...")
        out$GD=data.frame(cbind(as.data.frame(GT),as.data.frame(GD)))
        out$GM=GI
        out$G=myGenotype$G
        out$kinship=myGenotype$KI
        out$PCA=myGenotype$PC
        out$chor_taxa=chor_taxa
  }# is.null(Y)
# model_store=all.memo
if(!is.null(Y)) 
  {
    if(SNP.test&Multiple_analysis&DP$file.output)
      {
        all.memo=all.memo[!model_store%in%c("gBLUP","cBLUP","sBLUP")]
#BJK        if(length(all.memo)==0) break
        GMM=GAPIT.Multiple.Manhattan(model_store=all.memo,
                Y.names=colnames(Y)[-1],GM=IC$GM,seqQTN=DP$QTN.position,
                cutOff=DP$cutOff,plot.type=c("s"))
        print("GAPIT has output Multiple Manhattan figure with Symphysic type!!!")
        if(length(all.memo)*(ncol(Y)-1)>1&length(all.memo)*(ncol(Y)-1)<9)
          {
            print(all.memo)
            GMM=GAPIT.Multiple.Manhattan(model_store=all.memo,Y.names=colnames(Y)[-1],GM=IC$GM,seqQTN=QTN.position,cutOff=cutOff,plot.type=c("w","h"))
            print("GAPIT has output Multiple Manhattan figures with Wide and High types!!!")
            GAPIT.Circle.Manhattan.Plot(band=1,r=3,GMM$multip_mapP,plot.type=c("c","q"),signal.line=1,xz=GMM$xz,threshold=cutOff)
            print("GAPIT has output Multiple Manhattan and QQ figures with Circle types!!!")
          }
      } 
    if(!SNP.test|buspred)
      {     
        if(!is.null(testY)) GAPIT.PagainstP(Y=Y,testY=testY,model_store=model_store,traitname0=traitname0,lmpred=lmpred,type=Predict.type)
      }
    if(file.output&!SNP.test)
      { 
        model_store.gs=model_store[model_store%in%c("gBLUP","cBLUP","sBLUP")]  
        print("Here will start interactive for GS!!!")
        if(Inter.Plot)
          {
            GAPIT.Interactive.GS(model_store=model_store.gs,Y=Y)
            if(!is.null(testY))GAPIT.Interactive.GS(model_store=model_store.gs,Y=Y,testY=testY)
          }
      }# file.output&!SNP.test
  } # !is.null(Y)
options(warn = 0)
print("GAPIT has done all analysis!!!")
if(file.output) 
{
  print("Please find your all results in :")
  print(paste(getwd()))
}
return (out)
}  #end of GAPIT function
