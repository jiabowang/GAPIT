#'
#'
#'
#'
#'
#'
#'
`GAPIT.EMMAxP3D` <-function(
  ys,
  xs,
  K=NULL,
  Z=NULL,
  X0=NULL,
  CVI=NULL,
  CV.Inheritance=NULL,
  GI=NULL,
  GP=NULL,
	file.path=NULL,
  file.from=NULL,
  file.to=NULL,
  file.total=1, 
  genoFormat="Hapmap", 
  file.fragment=NULL,
  byFile=FALSE,
  fullGD=TRUE,
  SNP.fraction=1,
  file.G=NULL,
  file.Ext.G=NULL,
  GTindex=NULL,
  file.GD=NULL, 
  file.GM=NULL, 
  file.Ext.GD=NULL,
  file.Ext.GM=NULL,
  SNP.P3D=TRUE,
  Timmer,
  Memory,
  optOnly=TRUE,
  SNP.effect="Add",
  SNP.impute="Middle", 
  SNP.permutation=FALSE,
  ngrids=100,
  llim=-10,
  ulim=10,
  esp=1e-10,
  name.of.trait=NULL, 
  Create.indicator = FALSE, 
  Major.allele.zero = FALSE){
#Object: To esimate variance component by using EMMA algorithm and perform GWAS with P3D/EMMAx
#Output: ps, REMLs, stats, dfs, vgs, ves, BLUP,  BLUP_Plus_Mean, PEV
#Authors: Feng Tian, Alex Lipka and Zhiwu Zhang
# Last update: April 6, 2016
# Library used: EMMA (Kang et al, Genetics, Vol. 178, 1709-1723, March 2008)
# Note: This function was modified from the function of emma.REML.t from the library
##############################################################################################
#print("EMMAxP3D started...")
  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="P3D Start")
  Memory=GAPIT.Memory(Memory=Memory,Infor="P3D Start")


  #When numeric genotypes are selected, impute the missing SNPs with the allele indicated by the "SNP.impute" value
  if(!optOnly){
   if(SNP.impute == "Major") xs[which(is.na(xs))] = 2
   if(SNP.impute == "Minor") xs[which(is.na(xs))] = 0
   if(SNP.impute == "Middle") xs[which(is.na(xs))] = 1
  }


#--------------------------------------------------------------------------------------------------------------------<
  #Change data to matrix format if they are not
  if(is.null(dim(ys)) || ncol(ys) == 1)  ys <- matrix(ys, 1, length(ys))
  if(is.null(X0)) X0 <- matrix(1, ncol(ys), 1)

  #handler of special Z and K
  if(!is.null(Z)){ if(ncol(Z) == nrow(Z)) Z = NULL }
  if(!is.null(K)) {if(length(K)<= 1) K = NULL}

  #Extract dimension information
  g <- nrow(ys) #number of traits
  n <- ncol(ys) #number of observation

  q0 <- ncol(X0)#number of fixed effects
  q1 <- q0 + 1  #Nuber of fixed effect including SNP

  nr=n
  if(!is.null(K)) tv=ncol(K)

  #decomposation without fixed effect
  #print("Caling emma.eigen.L...")
  if(!is.null(K)) eig.L <- emma.eigen.L(Z, K) #this function handle both NULL Z and non-NULL Z matrix

  #eig.L$values[eig.L$values<0]=0
  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="eig.L")
  Memory=GAPIT.Memory(Memory=Memory,Infor="eig.L")

  #decomposation with fixed effect (SNP not included)
  #print("Calling emma.eigen.R.w.Z...")
  X <-  X0 #covariate variables such as population structure
  if(!is.null(Z) & !is.null(K)) eig.R <- try(emma.eigen.R.w.Z(Z, K, X),silent=TRUE) #This will be used to get REstricted ML (REML)
  if(is.null(Z)  & !is.null(K)) eig.R <- try(emma.eigen.R.wo.Z(   K, X),silent=TRUE) #This will be used to get REstricted ML (REML)

  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="eig.R")
  Memory=GAPIT.Memory(Memory=Memory,Infor="eig.R")

#eig.R$values[eig.R$values<0]=0
#print(labels(eig.R))
#print(length(eig.R$values))
#print(dim(eig.R$vectors))
#print("emma.eigen.R.w.Z called!!!")
#Handler of error in emma
#print("!!!!!!")
  if(!is.null(K)){
  if(inherits(eig.R, "try-error"))
       return(list(ps = NULL, REMLs = NA, stats = NULL,
                   effect.est = NULL, dfs = NULL, maf = NULL,
                   nobs = NULL, Timmer = Timmer, Memory=Memory,
                   vgs = NA, ves = NA, BLUP = NULL, BLUP_Plus_Mean = NULL,
                   PEV = NULL, BLUE=NULL
                   )
              )

#print("@@@@@")
   }
#-------------------------------------------------------------------------------------------------------------------->
#print("Looping through traits...")
#Loop on Traits
  for (j in 1:g){

  if(optOnly){

  #REMLE <- GAPIT.emma.REMLE(ys[j,], X, K, Z, ngrids, llim, ulim, esp, eig.R)
  #vgs <- REMLE$vg
  #ves <- REMLE$ve
  #REMLs <- REMLE$REML
  #REMLE_delta=REMLE$delta

    if(!is.null(K)){
      REMLE <- GAPIT.emma.REMLE(ys[j,], X, K, Z, ngrids, llim, ulim, esp, eig.R)

      Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="REML")
      Memory=GAPIT.Memory(Memory=Memory,Infor="REML")

      rm(eig.R)
      gc()
      Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="eig.R removed")
      Memory=GAPIT.Memory(Memory=Memory,Infor="eig.R removed")

      vgs <- REMLE$vg
      ves <- REMLE$ve
      REMLs <- REMLE$REML
      REMLE_delta=REMLE$delta

      rm(REMLE)
      gc()
    }


    vids <- !is.na(ys[j,])
    yv <- ys[j, vids]

    if(!is.null(Z) & !is.null(K))  U <- eig.L$vectors * matrix(c(sqrt(1/(eig.L$values + REMLE_delta)),rep(sqrt(1/REMLE_delta),nr - tv)),nr,((nr-tv)+length(eig.L$values)),byrow=TRUE)
    if( is.null(Z) & !is.null(K))  U <- eig.L$vectors * matrix(  sqrt(1/(eig.L$values + REMLE_delta)),nr,length(eig.L$values),byrow=TRUE)

    if( !is.null(Z) & !is.null(K)) eig.full.plus.delta <- as.matrix(c((eig.L$values + REMLE_delta), rep(REMLE_delta,(nr - tv))))
    if( is.null(Z) & !is.null(K))  eig.full.plus.delta <- as.matrix((eig.L$values + REMLE_delta))

    if(!is.null(K)){
      if(length(which(eig.L$values < 0)) > 0 ){
      #print("---------------------------------------------------The group kinship matrix at this compression level is not positive semidefinite. Please select another compression level.---------------------------------------------------")
             #return(list(ps = NULL, REMLs = 999999, stats = NULL, effect.est = NULL, dfs = NULL,maf=NULL,nobs = NULL,Timmer=Timmer,Memory=Memory,
             #vgs = 1.000, ves = 1.000, BLUP = NULL, BLUP_Plus_Mean = NULL,
             #PEV = NULL, BLUE=NULL))
      }
    }

  #Calculate the log likelihood function for the intercept only model

    X.int <- matrix(1,nrow(as.matrix(yv)),ncol(as.matrix(yv)))
    iX.intX.int <- solve(crossprod(X.int, X.int))
    iX.intY <- crossprod(X.int, as.matrix(as.matrix(yv)))
    beta.int <- crossprod(iX.intX.int, iX.intY)  #Note: we can use crossprod here becase iXX is symmetric
    X.int.beta.int <- X.int%*%beta.int


    logL0 <- 0.5*((-length(yv))*log(((2*pi)/length(yv))
                 *crossprod((yv-X.int.beta.int),(yv-X.int.beta.int)))
                  -length(yv))

    #print(paste("The value of logL0 inside of the optonly template is is",logL0, sep = ""))
  #print(paste("The value of nrow(as.matrix(ys[!is.na(ys)])) is ",nrow(as.matrix(ys[!is.na(ys)])), sep = ""))

    if(!is.null(K)){
      yt <- yt <- crossprod(U, yv)
      X0t <- crossprod(U, X0)

      X0X0 <- crossprod(X0t, X0t)
      X0Y <- crossprod(X0t,yt)
      XY <- X0Y

      iX0X0 <- try(solve(X0X0),silent=TRUE)
      if(inherits(iX0X0, "try-error")){
        iX0X0 <- MASS::ginv(X0X0)
        print("At least two of your covariates are linearly dependent. Please reconsider the covariates you are using for GWAS and GPS")
      }
      iXX <- iX0X0
    }

    if(is.null(K)){
      iXX <- try(solve(crossprod(X,X)),silent=TRUE)
      if(inherits(iXX, "try-error"))iXX <- MASS::ginv(crossprod(X,X))
      XY = crossprod(X,yv)
    }
    beta <- crossprod(iXX,XY) #Note: we can use crossprod here because iXX is symmetric
    X.beta <- X%*%beta
      
    beta.cv=beta
    BLUE=X.beta
      
    if(!is.null(K)){
      U.times.yv.minus.X.beta <- crossprod(U,(yv-X.beta))
      logLM <- 0.5*(-length(yv)*log(((2*pi)/length(yv))*crossprod(U.times.yv.minus.X.beta,U.times.yv.minus.X.beta))
                    - sum(log(eig.full.plus.delta)) - length(yv))
    }


    if(is.null(K)){
      U.times.yv.minus.X.beta <- yv-X.beta
      logLM <- 0.5*(-length(yv)*log(((2*pi)/length(yv))*crossprod(U.times.yv.minus.X.beta,U.times.yv.minus.X.beta)) - length(yv))
    }

  }#End if(optOnly)



#--------------------------------------------------------------------------------------------------------------------<
  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Trait")
  Memory=GAPIT.Memory(Memory=Memory,Infor="Trait")

  if(!is.null(K)){
    REMLE <- GAPIT.emma.REMLE(ys[j,], X, K, Z, ngrids, llim, ulim, esp, eig.R)

    Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="REML")
    Memory=GAPIT.Memory(Memory=Memory,Infor="REML")

    rm(eig.R)
    gc()
    Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="eig.R removed")
    Memory=GAPIT.Memory(Memory=Memory,Infor="eig.R removed")

    vgs <- REMLE$vg
    ves <- REMLE$ve
    REMLs <- REMLE$REML
    REMLE_delta=REMLE$delta

    rm(REMLE)
    gc()
  }
  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="REMLE removed")
  Memory=GAPIT.Memory(Memory=Memory,Infor="REMLE removed")

  if(!is.null(Z) & !is.null(K))  U <- eig.L$vectors * matrix(c(sqrt(1/(eig.L$values + REMLE_delta)),rep(sqrt(1/REMLE_delta),nr - tv)),nr,((nr-tv)+length(eig.L$values)),byrow=TRUE)
  if( is.null(Z) & !is.null(K))  U <- eig.L$vectors * matrix(  sqrt(1/(eig.L$values + REMLE_delta)),nr,length(eig.L$values),byrow=TRUE)

  if( !is.null(Z) & !is.null(K)) eig.full.plus.delta <- as.matrix(c((eig.L$values + REMLE_delta), rep(REMLE_delta,(nr - tv))))
  if( is.null(Z) & !is.null(K))  eig.full.plus.delta <- as.matrix((eig.L$values + REMLE_delta))



  if(!is.null(K)){
    if(length(which(eig.L$values < 0)) > 0 ){
 #print("---------------------------------------------------The group kinship matrix at this compression level is not positive semidefinite. Please select another compression level.---------------------------------------------------")
       #return(list(ps = NULL, REMLs = 999999, stats = NULL, effect.est = NULL, dfs = NULL,maf=NULL,nobs = NULL,Timmer=Timmer,Memory=Memory,
        #vgs = 1.000, ves = 1.000, BLUP = NULL, BLUP_Plus_Mean = NULL,
        #PEV = NULL, BLUE=NULL))
    }
  }


  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="U Matrix")
  Memory=GAPIT.Memory(Memory=Memory,Infor="U Matrix")

  if(SNP.P3D == TRUE)rm(eig.L)
  gc()

  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="eig.L removed")
  Memory=GAPIT.Memory(Memory=Memory,Infor="eig.L removed")

  #-------------------------------------------------------------------------------------------------------------------->

  #The cases that go though multiple file once
  file.stop=file.to
  if(optOnly) file.stop=file.from
  if(fullGD)  file.stop=file.from
  if(!fullGD & !optOnly) {print("Screening SNPs from file...")}

  #Add loop for genotype data files
  for (file in file.from:file.stop){
    Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="New Genotype file")
    Memory=GAPIT.Memory(Memory=Memory,Infor="New Genotype file")

    frag=1
    numSNP=file.fragment
    myFRG=NULL
    while(numSNP==file.fragment) {     #this is problematic if the read end at the last line

      #initial previous SNP storage
      x.prev <- vector(length = 0)

      #force to skip the while loop if optOnly
      if(optOnly) numSNP=0

      #Determine the case of first file and first fragment: skip read file
      if(file==file.from & frag==1& SNP.fraction<1){
        firstFileFirstFrag=TRUE
      }else{
        firstFileFirstFrag=FALSE
      }

      #In case of xs is not full GD, replace xs from file
      if(!fullGD & !optOnly & !firstFileFirstFrag ){

        Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Clean myFRG")
        Memory=GAPIT.Memory(Memory=Memory,Infor="Clean myFRG")

        #update xs for each file
        rm(xs)
        rm(myFRG)
        gc()
        print(paste("Current file: ",file," , Fragment: ",frag,sep=""))

        Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Read file fragment")
        Memory=GAPIT.Memory(Memory=Memory,Infor="Read file fragment")

    myFRG=GAPIT.Fragment( file.path = file.path,  
                          file.total = file.total,
                          file.G = file.G,
                          file.Ext.G = file.Ext.G,
#                          seed = seed,
                          SNP.fraction = SNP.fraction,
                          SNP.effect = SNP.effect,
                          SNP.impute = SNP.impute,
                          genoFormat = genoFormat,
                          file.GD = file.GD,
                          file.Ext.GD = file.Ext.GD,
                          file.GM = file.GM,
                          file.Ext.GM = file.Ext.GM,
                          file.fragment = file.fragment,
                          file = file,
                          frag = frag, 
                          Create.indicator = Create.indicator, 
                          Major.allele.zero = Major.allele.zero)

        Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Genotype file converted")
        Memory=GAPIT.Memory(Memory=Memory,Infor="Genotype file converted")

      #print("-----------------------------------------------------------------")

        if(is.null(myFRG$GD)){
          xs=NULL
        }else{
          xs=myFRG$GD[GTindex,]
        }


        if(!is.null(myFRG$GI))    {
          colnames(myFRG$GI)=c("SNP","Chromosome","Position")
          GI=as.matrix(myFRG$GI)
        }
        

        if(!is.null(myFRG$GI))    {
          numSNP=ncol(myFRG$GD)
        } else {
          numSNP=0
        }
        if(is.null(myFRG))numSNP=0  #force to end the while loop
  } # end of if(!fullGD)

  if(fullGD)numSNP=0  #force to end the while loop

#Skip REML if xs is from a empty fragment file
if(!is.null(xs))  {

   
  if(is.null(dim(xs)) || nrow(xs) == 1)  xs <- matrix(xs, length(xs),1)
  
  xs <- as.matrix(xs)
  
    if(length(which(is.na(xs)))>0){    #for the case where fragments are read in
     if(SNP.impute == "Major") xs[which(is.na(xs))] = 2
     if(SNP.impute == "Minor") xs[which(is.na(xs))] = 0
     if(SNP.impute == "Middle") xs[which(is.na(xs))] = 1
    }

  
  m <- ncol(xs) #number of SNPs
  t <- nrow(xs) #number of individuals
 
Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Before cleaning")
Memory=GAPIT.Memory(Memory=Memory,Infor="Before cleaning")
  #allocate spaces for SNPs
  #rm(dfs)
  #rm(stats)
  #rm(effect.est)
  #rm(ps)
  #rm(nobs)
  #rm(maf)
  #rm(rsquare_base)
  #rm(rsquare)
  #          rm(df)
  #          rm(tvalue)
  #          rm(stderr)
  gc()

Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="After cleaning")
Memory=GAPIT.Memory(Memory=Memory,Infor="After cleaning")

  dfs <- matrix(nrow = m, ncol = g)
  stats <- matrix(nrow = m, ncol = g)
  if(!Create.indicator) effect.est <- matrix(nrow = m, ncol = g)
  if(Create.indicator) effect.est <- NULL
  ps <- matrix(nrow = m, ncol = g)
  nobs <- matrix(nrow = m, ncol = g)
  maf <- matrix(nrow = m, ncol = g)
  rsquare_base <- matrix(nrow = m, ncol = g)
  rsquare <- matrix(nrow = m, ncol = g)
            df <- matrix(nrow = m, ncol = g)
            tvalue <- matrix(nrow = m, ncol = g)
            stderr <- matrix(nrow = m, ncol = g)
  #print(paste("Memory allocated.",sep=""))
  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Memory allocation")
  Memory=GAPIT.Memory(Memory=Memory,Infor="Memory allocation")

  if(optOnly)mloop=0
  if(!optOnly)mloop=m

  #Loop on SNPs
  #print(paste("Number of SNPs is ",mloop," in genotype file ",file, sep=""))

#set starting point of loop
if(file==file.from&frag==1){loopStart=0}else{loopStart=1}

for (i in loopStart:mloop){
#print(i)
#--------------------------------------------------------------------------------------------------------------------<
    normalCase=TRUE

    if((i >0)&(floor(i/1000)==i/1000)) {print(paste("Genotype file: ", file,", SNP: ",i," ",sep=""))}
    # To extract current snp. It save computation for next one in case they are identical
    if(i ==0&file==file.from&frag==1){
      #For the model without fitting SNP
      vids <- !is.na(ys[j,]) #### Feng changed
      xv <- ys[j, vids]*0+1 #### Feng changed
    }

    if(i >0 | file>file.from | frag>1){
      if(Create.indicator){ #I need create indicators and then calculate the minor allele frequency
       condition.temp <- unique(xs[vids,i])
       #Define what a bit is
       
       bit=nchar(as.character(xs[vids[1],i]))
       
       #Expand on the "which" statement below to include all instances of missing data
       
       if(bit==1)  condition <-  condition.temp[-which(condition.temp == "N")]
       if(bit==2)  condition <-  condition.temp[-which(condition.temp == "NN")]
       
       #print("condition.temp is ")
       #print(condition.temp)
                                                                                                                        
       #print("condition is")
       #print(condition)
       
       #print(paste("The value of i is ", i, sep = "")) 
        
       
       if(length(condition) <= 1){
        dfs[i, ] <- rep(NA, g)
        stats[i, ] <- rep(NA, g)
        effect.est <- rbind(effect.est, c(i,rep(NA, g), rep(NA, g)))
        ps[i, ] = rep(1, g)
        rsquare[i, ] <- rep(NA,g)
        rsquare_base[i, ]<-rep(NA,g)
        maf[i, ] <- rep(0, g)
                    df[i, ] <- rep(NA,g)
                    tvalue[i, ] <- rep(NA,g)
                    stderr[i, ] <- rep(NA,g)
        normalCase=FALSE
        x.prev= vector(length = 0)
       }
       
       }
       if(normalCase){
       #print("The head of xs[vids,i] is")
       #print(head(xs[vids,i]))
      
      if(Create.indicator){     #I need create indicators and then calculate the minor allele frequency
       
       indicator <-  GAPIT.Create.Indicator(xs[vids,i], SNP.impute = SNP.impute)
       xv <- indicator$x.ind
       vids <- !is.na(xv[,1]) #### Feng changed
      
       vids.TRUE=which(vids==TRUE)
       vids.FALSE=which(vids==FALSE)
       ns=nrow(xv)
       ss=sum(xv[,ncol(xv)])

       maf[i]=min(ss/ns,1-ss/ns)
       nobs[i]=ns
       
        q1 <- q0 + ncol(xv)    # This is done so that parameter estimates for all indicator variables are included

       
        #These two matrices need to be reinitiated for each SNP.
        Xt <- matrix(NA,nr, q1)
        iXX=matrix(NA,q1,q1)
       }       
      }
     
      if(!Create.indicator){ #### Feng changed
	   #print(xs[1:10,1:10])

       xv <- xs[vids,i]
       vids <- !is.na(xs[,i]) #### Feng changed
      
       vids.TRUE=which(vids==TRUE)
       vids.FALSE=which(vids==FALSE)
       ns=length(xv)
	   #print(xv))
       ss=sum(xv)

       maf[i]=min(.5*ss/ns,1-.5*ss/ns)
       nobs[i]=ns
      }

     nr <- sum(vids)
     if(i ==1 & file==file.from&frag==1 & !Create.indicator) {
       Xt <- matrix(NA,nr, q1)
       iXX=matrix(NA,q1,q1)
     }

    }

    #Situation of no variation for SNP except the fisrt one(synthetic for EMMAx/P3D)
    if((min(xv) == max(xv)) & (i >0 | file>file.from |frag>1))
    {
      dfs[i, ] <- rep(NA, g)
      stats[i, ] <- rep(NA, g)
      if(!Create.indicator) effect.est[i,] <- rep(NA, g)
      if(Create.indicator) effect.est <- rbind(effect.est, c(i,rep(NA, g),rep(NA, g)))
      ps[i, ] = rep(1, g)
      rsquare[i, ] <- rep(NA,g)
      rsquare_base[i, ]<-rep(NA,g)
                df[i, ] <- rep(NA,g)
                tvalue[i, ] <- rep(NA,g)
                stderr[i, ] <- rep(NA,g)
      normalCase=FALSE
    }else if(identical(x.prev, xv))     #Situation of the SNP is identical to previous
    {
      if(i >1 | file>file.from | frag>1){
        dfs[i, ] <- dfs[i - 1, ]
        stats[i, ] <- stats[i - 1, ]
        if(!Create.indicator) effect.est[i, ] <- effect.est[i - 1, ]
        if(Create.indicator) effect.est <- rbind(effect.est, c(i, rep(NA, g), rep(NA, g))) #If the previous SNP is idnetical, indicate this by "NA"
        ps[i, ] <- ps[i - 1, ]
        rsquare[i, ] <- rsquare[i - 1, ]
        rsquare_base[i, ] <-rsquare_base[i - 1, ]
                  df[i, ] <- df[i - 1, ]
                  tvalue[i, ] <- tvalue[i - 1, ]
                  stderr[i, ] <- stderr[i - 1, ]
        normalCase=FALSE
      }
    }
#-------------------------------------------------------------------------------------------------------------------->
  if(i == 0 &file==file.from &frag==1){

   #Calculate the log likelihood function for the intercept only model

   #vids <- !is.na(ys[j,])
   yv <- ys[j, vids]

   X.int <- matrix(1,nrow(as.matrix(yv)),ncol(as.matrix(yv)))
   iX.intX.int <- solve(crossprod(X.int, X.int))
   iX.intY <- crossprod(X.int, as.matrix(as.matrix(yv)))
   beta.int <- crossprod(iX.intX.int, iX.intY)  #Note: we can use crossprod here becase iXX is symmetric
   X.int.beta.int <- X.int%*%beta.int




   #X.int <- matrix(1,nrow(as.matrix(ys[!is.na(ys)])),ncol(as.matrix(ys[!is.na(ys)])))
   #iX.intX.int <- solve(crossprod(X.int, X.int))
   #iX.intY <- crossprod(X.int, as.matrix(ys[!is.na(ys)]))
   #beta.int <- crossprod(iX.intX.int, iX.intY)  #Note: we can use crossprod here becase iXX is symmetric
   #X.int.beta.int <- X.int%*%beta.int

   logL0 <- 0.5*((-length(yv))*log(((2*pi)/length(yv))
                 *crossprod((yv-X.int.beta.int),(yv-X.int.beta.int)))
                  -length(yv))


   #logL0 <- 0.5*((-nrow(as.matrix(ys[!is.na(ys)])))*log(((2*pi)/nrow(ys))
   # *crossprod(((as.matrix(ys[!is.na(ys)]))-X.int.beta.int),((as.matrix(ys[!is.na(ys)]))-X.int.beta.int)))
   # -nrow(as.matrix(ys[!is.na(ys)])))

    #print(paste("The value of logL0 inside of the calculating SNPs loop is", logL0, sep = ""))
   }

    #Normal case
    if(normalCase)
    {

#--------------------------------------------------------------------------------------------------------------------<
      #nv <- sum(vids)
      yv <- ys[j, vids] #### Feng changed
      nr <- sum(vids) #### Feng changed
      if(!is.null(Z) & !is.null(K))
      {
        r<- ncol(Z) ####Feng, add a variable to indicate the number of random effect
        vran <- vids[1:r] ###Feng, add a variable to indicate random effects with nonmissing genotype
        tv <- sum(vran)  #### Feng changed
      }



#-------------------------------------------------------------------------------------------------------------------->

#--------------------------------------------------------------------------------------------------------------------<



      if(i >0 | file>file.from|frag>1) dfs[i, j] <- nr - q1
    	if(i >0 | file>file.from|frag>1){ 
        if(!Create.indicator) X <- cbind(X0[vids, , drop = FALSE], xs[vids,i])
        if(Create.indicator){
          X <- cbind(X0[vids, , drop = FALSE], xv)
          #if(i == 1) {print("the head of X for running GWAS is")}
          #if(i == 1) {print(head(X))}
        }       
        
      } 
       #Recalculate eig and REML if not using P3D  NOTE THIS USED TO BE BEFORE the two solid lines
      if(SNP.P3D==FALSE & !is.null(K))
      {
        if(!is.null(Z)) eig.R <- emma.eigen.R.w.Z(Z, K, X) #This will be used to get REstricted ML (REML)
        if(is.null(Z)) eig.R <- emma.eigen.R.wo.Z( K, X) #This will be used to get REstricted ML (REML)
        if(!is.null(Z)) REMLE <- GAPIT.emma.REMLE(ys[j,], X, K, Z, ngrids, llim, ulim, esp, eig.R)
        if(is.null(Z)) REMLE <- GAPIT.emma.REMLE(ys[j,], X, K, Z = NULL, ngrids, llim, ulim, esp, eig.R)
        if(!is.null(Z) & !is.null(K)) U <- eig.L$vectors * matrix(c(sqrt(1/(eig.L$values + REMLE$delta)),rep(sqrt(1/REMLE$delta),nr - tv)),nr,((nr-tv)+length(eig.L$values)),byrow=TRUE)
        if(is.null(Z) & !is.null(K)) U <- eig.L$vectors * matrix( sqrt(1/(eig.L$values + REMLE$delta)),nr,length(eig.L$values),byrow=TRUE)

        vgs <- REMLE$vg
        ves <- REMLE$ve
        REMLs <- REMLE$REML
        REMLE_delta=REMLE$delta

      }

      if(n==nr)
      {
        if(!is.null(K))
        {
            yt <- crossprod(U, yv)
            if(i == 0 &file==file.from &frag==1){
             X0t <- crossprod(U, X0)
             Xt <- X0t
            }
            if(i > 0 | file>file.from |frag>1){
              #if(i ==1 & file==file.from&frag==1) Xt <- matrix(NA,nr, q1)
             
             if(Create.indicator){
                xst <- crossprod(U, X[,(q0+1):q1])
                Xt[1:nr,1:q0] <- X0t
                Xt[1:nr,(q0+1):q1] <- xst
               
             }
             
              #print(paste("i:",i,"q0:",q0,"q1:",q1,"nt:",nr,"XT row",nrow(Xt),"XT col",ncol(Xt),sep=" "))
             if(!Create.indicator){
                xst <- crossprod(U, X[,ncol(X)])
                Xt[1:nr,1:q0] <- X0t
                Xt[1:nr,q1] <- xst
             }
            }
        }else{
        yt=yv
        if(i == 0 &file==file.from &frag==1) X0t <- X0
        if(i > 0 | file>file.from |frag>1) xst <- X[,ncol(X)]
        }

        if(i == 0 &file==file.from &frag==1){
         X0X0 <- crossprod(X0t, X0t)
         #XX <- X0X0
        }
        if(i > 0 | file>file.from |frag>1){
         #if(i == 1)XX=matrix(NA,q1,q1)


         X0Xst <- crossprod(X0t,xst)
         XstX0 <- t(X0Xst)
         xstxst <- crossprod(xst, xst)
         # if(i == 1){
         # Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Calculate_X0Xst_XstX0_xstxst")
         # Memory=GAPIT.Memory(Memory=Memory,Infor="Calculate_X0Xst_XstX0_xstxst")
         # }
         #XX <- rbind(cbind(X0X0, X0Xst), cbind(XstX0, xstxst))

         #XX[1:q0,1:q0] <- X0X0
         #XX[q1,1:q0] <- X0Xst
         #XX[1:q0,q1] <- X0Xst
         #XX[q1,q1] <- xstxst


        }


        if(X0X0[1,1] == "NaN")
        {
          Xt[which(Xt=="NaN")]=0
          yt[which(yt=="NaN")]=0
          XX=crossprod(Xt, Xt)
        }
        if(i == 0 &file==file.from & frag==1){
         X0Y <- crossprod(X0t,yt)
         XY <- X0Y
        }
        if(i > 0 | file>file.from |frag>1){
         xsY <- crossprod(xst,yt)
         XY <- c(X0Y,xsY)
# if(i == 1){
# Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Calculate_xsY_X0Y")
# Memory=GAPIT.Memory(Memory=Memory,Infor="Calculate_xsY_X0Y")
# }
        }
        #XY = crossprod(Xt,yt)
      }

      #Missing SNP
      if(n>nr)
      {
       UU=crossprod(U,U)
       A11=UU[vids.TRUE,vids.TRUE]
       A12=UU[vids.TRUE,vids.FALSE]
       A21=UU[vids.FALSE,vids.TRUE]
       A22=UU[vids.FALSE,vids.FALSE]
       A22i =try(solve(A22),silent=TRUE )
       if(inherits(A22i, "try-error")) A22i <- MASS::ginv(A22)

       F11=A11-A12%*%A22i%*%A21
       XX=crossprod(X,F11)%*%X
       XY=crossprod(X,F11)%*%yv
      }
      if(i == 0 &file==file.from &frag==1){
       iX0X0 <- try(solve(X0X0),silent=TRUE)
       if(inherits(iX0X0, "try-error")){
         iX0X0 <- MASS::ginv(X0X0)
         print("At least two of your covariates are linearly dependent. Please reconsider the covariates you are using for GWAS and GPS")
       }
       iXX <- iX0X0
      }
      if(i > 0 | file>file.from |frag>1){

      #if(i ==1 &file==file.from &frag==1) iXX=matrix(NA,q1,q1)
        if(Create.indicator){
          B22 <- xstxst - XstX0%*%iX0X0%*%X0Xst
          invB22 <- solve(B22)
          B21 <- tcrossprod(XstX0, iX0X0)
          NeginvB22B21 <- crossprod(-invB22,B21)
          B11 <- iX0X0 + as.numeric(invB22)*crossprod(B21,B21)



          iXX[1:q0,1:q0]=B11
          iXX[(q0+1):q1,(q0+1):q1]=solve(B22)  
          iXX[(q0+1):q1,1:q0]=NeginvB22B21
          iXX[1:q0,(q0+1):q1]=t(NeginvB22B21)

        }

      
        if(!Create.indicator){
          B22 <- xstxst - XstX0%*%iX0X0%*%X0Xst
          invB22 <- 1/B22
          #B12 <- crossprod(iX0X0,X0Xst)
          B21 <- tcrossprod(XstX0, iX0X0)
          NeginvB22B21 <- crossprod(-invB22,B21)
          #B11 <- iX0X0 + B12%*%invB22%*%B21
          B11 <- iX0X0 + as.numeric(invB22)*crossprod(B21,B21)
          #iXX <- rbind(cbind(B11,t(NeginvB22B21)), cbind(NeginvB22B21,invB22))

          iXX[1:q0,1:q0]=B11
          iXX[q1,q1]=1/B22
          iXX[q1,1:q0]=NeginvB22B21
          iXX[1:q0,q1]=NeginvB22B21
          

        }
        #if(i == 1){
        # Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Calculate_iXX")
        # Memory=GAPIT.Memory(Memory=Memory,Infor="Calculate_iXX")
        #}

      }

      if(is.null(K)){
        iXX <- try(solve(crossprod(X,X)),silent=TRUE)
        if(inherits(iXX, "try-error"))iXX <- MASS::ginv(crossprod(X,X))
        XY = crossprod(X,yv)
      }

      #iXX <- try(solve(XX))
      #if(inherits(iXX, "try-error")) iXX <- ginv(crossprod(Xt, Xt))
      #print("The dimension if iXX is")
      #print(dim(iXX))
      #print("The length of XY is")
      #print(length(XY))
      
      beta <- crossprod(iXX,XY) #Note: we can use crossprod here becase iXX is symmetric
      #print("beta was estimated")

#-------------------------------------------------------------------------------------------------------------------->

       
#--------------------------------------------------------------------------------------------------------------------<
      if(i ==0 &file==file.from &frag==1 & !is.null(K))
      {
        Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="ReducedModel")
Memory=GAPIT.Memory(Memory=Memory,Infor="ReducdModel")

        #beta.cv=beta

        
        
        XtimesBetaHat <- X%*%beta

        YminusXtimesBetaHat <- ys[j,]- XtimesBetaHat
        vgK <- vgs*K
        Dt <- crossprod(U, YminusXtimesBetaHat)

        if(!is.null(Z)) Zt <- crossprod(U, Z)
        if(is.null(Z)) Zt <- t(U)

        if(X0X0[1,1] == "NaN")
        {
        Dt[which(Dt=="NaN")]=0
        Zt[which(Zt=="NaN")]=0
        }

        BLUP <- K %*% crossprod(Zt, Dt) #Using K instead of vgK because using H=V/Vg

#print("!!!!")
      #Clean up the BLUP starf to save memory
      Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="before Dt clean")
      Memory=GAPIT.Memory(Memory=Memory,Infor="before Dt clean")
#      rm(Dt)
      gc()
      Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Dt clean")
      Memory=GAPIT.Memory(Memory=Memory,Infor="Dt clean")





        grand.mean.vector <- rep(beta[1], length(BLUP))
        BLUP_Plus_Mean <- grand.mean.vector + BLUP
    	  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="BLUP")
        Memory=GAPIT.Memory(Memory=Memory,Infor="BLUP")

        #PEV
        C11=try(vgs*solve(crossprod(Xt,Xt)),silent=TRUE)
        if(inherits(C11, "try-error")) C11=vgs*MASS::ginv(crossprod(Xt,Xt))

        C21=-K%*%crossprod(Zt,Xt)%*%C11
		Kinv=try(solve(K)  ,silent=TRUE  ) 
        if(inherits(Kinv, "try-error")) Kinv=MASS::ginv(K)
        
        if(!is.null(Z)) term.0=crossprod(Z,Z)/ves
        if(is.null(Z)) term.0=diag(1/ves,nrow(K))

        term.1=try(solve(term.0+Kinv/vgs ) ,silent=TRUE )
        if(inherits(term.1, "try-error")) term.1=MASS::ginv(term.0+Kinv/vgs )

        term.2=C21%*%crossprod(Xt,Zt)%*%K
        C22=(term.1-term.2 )
        PEV=as.matrix(diag(C22))
 #print(paste("The value of is.na(CVI) is", is.na(CVI),  sep = ""))
        # CVI may be > 1 element long
        #if(!is.na(CVI)){
        if(any(!is.na(CVI))){
    		  XCV=as.matrix(cbind(1,data.frame(CVI[,-1])))
	
      		#CV.Inheritance specified
		      beta.Inheritance=beta
		      if(!is.null(CV.Inheritance)){
			      XCV=XCV[,1:(1+CV.Inheritance)]
			      beta.Inheritance=beta[1:(1+CV.Inheritance)]
		      }
		#Interception only
		if(length(beta)==1)XCV=X
		
        BLUE=try(XCV%*%beta.Inheritance,silent=TRUE)
        if(inherits(BLUE, "try-error")) BLUE = NA
     #print("GAPIT just after BLUE")
     Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="PEV")
        Memory=GAPIT.Memory(Memory=Memory,Infor="PEV")

      }#end of if(i ==0&file==file.from   & !is.null(K))

        # CVI may be > 1 element long.
        #if(is.na(CVI)) BLUE = NA
        if(any(is.na(CVI))) BLUE = NA
}#end if(!is.na(CVI))
#-------------------------------------------------------------------------------------------------------------------->

#--------------------------------------------------------------------------------------------------------------------<
      if(i ==0 &file==file.from &frag==1 & is.null(K))
      {
        YY=crossprod(yt, yt)
        ves=(YY-crossprod(beta,XY))/(n-q0)
        r=yt-X%*%iXX%*%XY
        REMLs=-.5*(n-q0)*log(det(ves)) -.5*n -.5*(n-q0)*log(2*pi)
# REMLs=-.5*n*log(det(ves)) -.5*log(det(iXX)/ves) -.5*crossprod(r,r)/ves -.5*(n-q0)*log(2*pi)
        vgs = 0
        BLUP = 0
        BLUP_Plus_Mean = NaN
        PEV = ves
        #print(paste("X row:",nrow(X)," col:",ncol(X)," beta:",length(beta),sep=""))
XCV=as.matrix(cbind(1,data.frame(CVI[,-1])))

#CV.Inheritance specified
beta.Inheritance=beta
if(!is.null(CV.Inheritance)){
XCV=XCV[,1:(1+CV.Inheritance)]
beta.Inheritance=beta[1:(1+CV.Inheritance)]
}
#Interception only
if(length(beta)==1)XCV=X


        #BLUE=XCV%*%beta.Inheritance   modified by jiabo wang 2016.11.21
        BLUE=try(XCV%*%beta.Inheritance,silent=TRUE)
        if(inherits(BLUE, "try-error")) BLUE = NA

      }


#Clean up the BLUP stuff to save memory
if(i ==0 &file==file.from &frag==1 & !is.null(K))
{
     Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="K normal")
        Memory=GAPIT.Memory(Memory=Memory,Infor="K normal")
if(SNP.P3D == TRUE) K=1  #NOTE: When SNP.P3D == FALSE, this line will mess up the spectral decomposition of the kinship matrix at each SNP.
#rm(Dt)
rm(Zt)            
rm(Kinv)
rm(C11)
rm(C21)
rm(C22)

gc()
     Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="K set to 1")
        Memory=GAPIT.Memory(Memory=Memory,Infor="K set to 1")
}

      if(i == 0 &file==file.from & frag==1){
      beta.cv=beta
      X.beta <- X%*%beta

      if(!is.null(K)){
              U.times.yv.minus.X.beta <- crossprod(U,(yv-X.beta))
              logLM_Base <- 0.5*(-length(yv)*log(((2*pi)/length(yv))*crossprod(U.times.yv.minus.X.beta,U.times.yv.minus.X.beta))
                    - sum(log(eig.full.plus.delta)) - length(yv))

      }
      if(is.null(K)){
              U.times.yv.minus.X.beta <- yv-X.beta
              logLM_Base <- 0.5*(-length(yv)*log(((2*pi)/length(yv))*crossprod(U.times.yv.minus.X.beta,U.times.yv.minus.X.beta)) - length(yv))
      }
      rsquare_base_intitialized <- 1-exp(-(2/length(yv))*(logLM_Base-logL0))

      }


      #print(Create.indicator)
      #calculate t statistics and P-values
      if(i > 0 | file>file.from |frag>1)
      {
       if(!Create.indicator){
        #if(i<5)print(beta[q1])
        #if(i<5)print(iXX[q1, q1])
        if(!is.null(K)) stats[i, j] <- beta[q1]/sqrt(iXX[q1, q1] *vgs) 
        if(is.null(K)) stats[i, j] <- beta[q1]/sqrt(iXX[q1, q1] *ves)
        effect.est[i, ] <- beta[q1]
        ps[i, ] <- 2 * stats::pt(abs(stats[i, ]), dfs[i, ],lower.tail = FALSE)
        if(is.na(ps[i,]))ps[i,]=1
        #print(c(i,ps[i,],stats[i,],beta[q1],iXX[q1, q1]))
       } 
       if(Create.indicator){
       
        F.num.first.two <- crossprod(beta[(q0+1):q1], solve(iXX[(q0+1):q1,(q0+1):q1]))
        if(!is.null(K)) stats[i, j] <- (F.num.first.two %*% beta[(q0+1):q1])/(length((q0+1):q1)*vgs)
        if(is.null(K)) stats[i, j] <- (F.num.first.two %*% beta[(q0+1):q1])/(length((q0+1):q1)*ves)
        effect.est <- rbind(effect.est, cbind(rep(i,length((q0+1):q1)), indicator$unique.SNPs, beta[(q0+1):q1])) #Replace with rbind
        ps[i, ] <- stats::pf(stats[i, j], df1=length((q0+1):q1), df2=(nr-ncol(X)), lower.tail = FALSE) #Alex, are these denominator degrees of freedom correct?
        dfs[i,] <- nr-nrow(X)
        
       }
              #Calculate the maximum full likelihood function value and the r square

      X.beta <- X%*%beta
      if(!is.null(K)){
          U.times.yv.minus.X.beta <- crossprod(U,(yv-X.beta))
          logLM <- 0.5*(-length(yv)*log(((2*pi)/length(yv))*crossprod(U.times.yv.minus.X.beta,U.times.yv.minus.X.beta))
                       - sum(log(eig.full.plus.delta))- length(yv))
      }
      if(is.null(K)){
            U.times.yv.minus.X.beta <- yv-X.beta
            logLM <- 0.5*(-length(yv)*log(((2*pi)/length(yv))*crossprod(U.times.yv.minus.X.beta,U.times.yv.minus.X.beta)) - length(yv))
      }

      rsquare_base[i, ] <- rsquare_base_intitialized
      rsquare[i, ] <- 1-exp(-(2/length(yv))*(logLM-logL0))

                  #Calculate df, t value and standard error _xiaolei changed
                  df[i,] <- dfs[i,]

                  tvalue[i,] <- stats[i, j]
                  stderr[i,] <- beta[ncol(CVI)+1]/stats[i, j]
                  #stderr[i,] <- sqrt(vgs)
                  # modified by Jiabo at 20191115
      }
      #print("!!!!!!!!!!!!!!!")
      #print(Create.indicator)
#-------------------------------------------------------------------------------------------------------------------->

    } # End of if(normalCase)
    x.prev=xv #update SNP

} # End of loop on SNPs

Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Screening SNPs")
Memory=GAPIT.Memory(Memory=Memory,Infor="Screening SNPs")
# print(head(tvalue))
# print(head(stderr))
# print(head(effect.est))
#output p value for the genotype file
if(!fullGD)
{ 
  #print("!!!!!!!!!!")
  #print(dim(GI))
  utils::write.table(GI, paste("GAPIT.TMP.GI.",name.of.trait,file,".",frag,".txt",sep=""), quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
  utils::write.table(ps, paste("GAPIT.TMP.ps.",name.of.trait,file,".",frag,".txt",sep=""), quote = FALSE, sep = "\t", row.names = FALSE,col.names = FALSE)
  utils::write.table(maf, paste("GAPIT.TMP.maf.",name.of.trait,file,".",frag,".txt",sep=""), quote = FALSE, sep = "\t", row.names = FALSE,col.names = FALSE)
  utils::write.table(nobs, paste("GAPIT.TMP.nobs.",name.of.trait,file,".",frag,".txt",sep=""), quote = FALSE, sep = "\t", row.names = FALSE,col.names = FALSE)
  utils::write.table(rsquare_base, paste("GAPIT.TMP.rsquare.base.",name.of.trait,file,".",frag,".txt",sep=""), quote = FALSE, sep = "\t", row.names = FALSE,col.names = FALSE)
  utils::write.table(rsquare, paste("GAPIT.TMP.rsquare.",name.of.trait,file,".",frag,".txt",sep=""), quote = FALSE, sep = "\t", row.names = FALSE,col.names = FALSE)
  utils::write.table(df, paste("GAPIT.TMP.df.",name.of.trait,file,".",frag,".txt",sep=""), quote = FALSE, sep = "\t", row.names = FALSE,col.names = FALSE)
  utils::write.table(tvalue, paste("GAPIT.TMP.tvalue.",name.of.trait,file,".",frag,".txt",sep=""), quote = FALSE, sep = "\t", row.names = FALSE,col.names = FALSE)
  utils::write.table(stderr, paste("GAPIT.TMP.stderr.",name.of.trait,file,".",frag,".txt",sep=""), quote = FALSE, sep = "\t", row.names = FALSE,col.names = FALSE)
  utils::write.table(effect.est, paste("GAPIT.TMP.effect.est.",name.of.trait,file,".",frag,".txt",sep=""), quote = FALSE, sep = "\t", row.names = FALSE,col.names = FALSE)
 
  #rm(dfs,stats,ps,nobs,maf,GI)   #This cause problem on return
  #gc()
}
 
    frag=frag+1   #Progress to next fragment

} #end of if(!is.null(X))

} #end of repeat on fragment



} # Ebd of loop on file
  } # End of loop on traits

  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="GWAS done for this Trait")
  Memory=GAPIT.Memory(Memory=Memory,Infor="GWAS done for this Trait")
  #print("GAPIT.EMMAxP3D accomplished successfully!")

    return(list(ps = ps, REMLs = -2*REMLs, stats = stats, effect.est = effect.est, rsquare_base = rsquare_base, rsquare = rsquare, dfs = dfs, df = df, tvalue = tvalue, stderr = stderr,maf=maf,nobs = nobs,Timmer=Timmer,Memory=Memory,
        vgs = vgs, ves = ves, BLUP = BLUP, BLUP_Plus_Mean = BLUP_Plus_Mean,
        PEV = PEV, BLUE=BLUE, logLM = logLM,effect.snp=effect.est,effect.cv=beta.cv))

}#end of GAPIT.EMMAxP3D function
#=============================================================================================
