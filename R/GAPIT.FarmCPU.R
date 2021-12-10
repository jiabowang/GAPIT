#'
#'
#'

`FarmCPU.0000` <-
function(){
    #################################################################
    #FarmCPU: Fixed and random model Circuitous Probability Unification
    #This is an R package to perform GWAS and genome prediction
    #Designed by Zhiwu Zhang
    #Writen by Xiaolei Liu and Zhiwu Zhang
    #Thanks for Aaron Kusmec pointing out the bug in 'FarmCPU.Burger' function
    FarmCPU.Version="FarmCPU v1.02, Dec 21, 2016"
    return(FarmCPU.Version)
}

`FarmCPU.BIN` <-function(
    Y = NULL,
    GDP = NULL,
    GM = NULL,
    CV = NULL,
    P = NULL,
    orientation = "col",
    method = "random",
    b = c(5e5,5e6,5e7),
    s = seq(10,100,10), theLoop = NULL, bound = NULL){
    #Input: Y - n by 2 matrix with fist column as taxa name and second as trait
    #Input: GDP - n by m+1 matrix. The first colum is taxa name. The rest are m genotype
    #Input: GM - m by 3  matrix for SNP name, chromosome and BP
    #Input: CV - n by t matrix for t covariate variables.
    #Input: P - m by 1 matrix containing probability
    #Input: method - options are "static", "optimum", and "integral"
    #Input: b - vecter of length>=1 for bin size
    #Input: s - vecter of length>=1 for size of complexity (number of QTNs)
    #Requirement: Y, GDP and CV have same taxa order. GDP and GM have the same order on SNP
    #Requirement: P and GM are in the same order
    #Requirement: No missing data
    #Output: bin - n by s matrix of genotype
    #Output: binmap - s by 3 matrix for map of bin
    #Output: seqQTN - s by 1 vecter for index of QTN on GM (+1 for GDP column wise)
    #Relationship: bin=GDP[,c(seqQTN)], binmap=GM[seqQTN,]
    #Authors: Zhiwu Zhang
    # Last update: Febuary 28, 2013
    ##############################################################################
    #print("FarmCPU.BIN Started")
    
    #print("bin size")
    #print(b)
    #print("bin selection")
    #print(s)
    
    #print("method specified:")
    #print(method)
    if(is.null(P)) return(list(bin=NULL,binmap=NULL,seqQTN=NULL))
    
    #Set upper bound for bin selection to squareroot of sample size
    
    n=nrow(Y)
    #bound=round(sqrt(n)/log10(n))
    if(is.null(bound)){
        bound=round(sqrt(n)/sqrt(log10(n)))
    }
    #bound=round(sqrt(n))
    #bound=round(n/log10(n))
    #bound=n-1
    s[s>bound]=bound
    s=unique(s[s<=bound]) #keep the within bound
    
    #print("number of bins allowed")
    #print(s)
    
    optimumable=(length(b)*length(s)>1)
    if(!optimumable & method=="optimum"){
        #print("Warning: method was changed from optimum to static")
        method="static"
    }
    
    #print("method actually used:")
    #print(method)
    
    #Method of random
    #if(method=="random") seqQTN=sample(nrow(GM),s) #this is for test only
    #Method of static
    if(method=="static"){
        #print("Via static")
        if(theLoop==2){
            b=b[3]
        }else if(theLoop==3){
            b=b[2]
        }else{
            b=b[1]
        }
        s=bound
        #b=median(b)
        #s=median(s)
        s[s>bound]=bound
        #print("Bin : bin.size, bin.selection")
        #print(c(b,s))
        print("optimizing possible QTNs...")
        GP=cbind(GM,P,NA,NA,NA)
        mySpecify=GAPIT.Specify(GI=GM,GP=GP,bin.size=b,inclosure.size=s)
        seqQTN=which(mySpecify$index==TRUE)
        #print("Bin set through static")
    }
    #Method of optimum
    #============================optimum start============================================
    if(method=="optimum"&optimumable){
        #print("optimizing bins")
        #print("c(bin.size, bin.selection, -2LL, VG, VE)")
        print("optimizing possible QTNs...")
        count=0
        for (bin in b){
            for (inc in s){
                count=count+1
                GP=cbind(GM,P,NA,NA,NA)
                #print("debug in bin 000")
                
                #print("calling Specify")
                #print(date())
                
                mySpecify=GAPIT.Specify(GI=GM,GP=GP,bin.size=bin,inclosure.size=inc)
                
                #print("calling Specify done")
                #print(date())
                
                seqQTN=which(mySpecify$index==TRUE)
                #print("seqQTN")
                #print(seqQTN)
                if(orientation=="col"){
                    if(bigmemory::is.big.matrix(GDP)){
                        GK=bigmemory::deepcopy(GDP,cols=seqQTN)
                    }else{
                        GK=GDP[,seqQTN] #GK has the first as taxa in FarmCPU.Burger. But not get uesd.
                        #GK=GDP[,seqQTN]
                    }
                }else{
                    #if(is.big.matrix(GDP)){
                    #GK=bigmemory::deepcopy(GDP,rows=seqQTN)
                    #GK=t(GK)
                    #}else{
                    #GK=cbind(Y[,1],t(GDP[c(1,seqQTN),])) #GK has the first as taxa in FarmCPU.Burger. But not get uesd.
                    #some problem here
                    GK=t(GDP[seqQTN,])
                    #}
                }
                
                #print("GK")
                #print(GK)
                #print("calling Burger")
                #print(date())
                
                myBurger=FarmCPU.Burger(Y=Y[,1:2],CV=CV,GK=GK)
                
                #print("calling Burger done")
                #print(date())
                
                myREML=myBurger$REMLs
                myVG=myBurger$vg #it is unused
                myVE=myBurger$ve #it is unused
                
                #print("c(bin.size, bin.selection, -2LL, VG, VE)")
                print(c(bin,inc,myREML,myVG,myVE))
                #Recoding the optimum GK
                if(count==1){
                    seqQTN.save=seqQTN
                    LL.save=myREML
                    bin.save=bin
                    inc.save=inc
                    vg.save=myVG  # for genetic variance
                    ve.save=myVE  # for residual variance
                }else{
                    if(myREML<LL.save){
                        seqQTN.save=seqQTN
                        LL.save=myREML
                        bin.save=bin
                        inc.save=inc
                        vg.save=myVG  # for genetic variance
                        ve.save=myVE  # for residual variance
                    }
                } #end of if(count==1)
            }#loop on bin number
        }#loop on bin size
        seqQTN=seqQTN.save
        #ve.save=ve.save
        #vg.save=vg.save
        #print(seqQTN)
        #print("Bin optimized: -2LL, bin.size, bin.selection")
        #print(c(LL.save,bin.save,inc.save))
        #print(LL.save)
        #print("bin.save")
        #print(bin.save)
        #print("inc.save")
        #print(inc.save)
    }
    #============================end of optimum============================================
    
    bin=NULL
    binmap=NULL
    #The following are commented out as they will be finalized in Remove function
    #if(orientation=="col"){
    #  bin=GDP[,seqQTN]
    #}else{
    #  bin=t(GDP[seqQTN,] )
    #}
    #binmap=GM[seqQTN,]
    #print(length(seqQTN))
    
    #print("FarmCPU.Bin accomplished successfully!")
    return(list(bin=bin,binmap=binmap,seqQTN=seqQTN))
}#The function FarmCPU.BIN ends here
`FarmCPU.GLM` <-
function(Y=NULL,GDP=NULL,GM=NULL,CV=NULL,orientation="row",package="FarmCPU.LM",model="A",ncpus=1,seqQTN=NULL,npc=0){
    #Object: To perform GWAS with GLM model
    #Input: Y - n by 2 matrix with fist column as taxa name and second as trait
    #Input: GDP - n by m matrix. This is Genotype Data Pure (GDP). THERE IS NOT COLUMN FOR TAXA.
    #Input: GM - m by 3  matrix for SNP name, chromosome and BP
    #Input: CV - n by t matrix for t covariate variables.
    #Requirement: Y, GDP and CV have same taxa order. GDP and GM have the same order on SNP
    #Output: P - m by 4(t+1) matrix containing estimate, tvalue, stderr and pvalue for covariates and SNP
    #Authors: Xiaolei Liu and Zhiwu Zhang
    # Last update: may 9, 2012
    ##############################################################################
    if(is.null(Y)) return(NULL)  #Y is required
    if(is.null(GDP) & is.null(CV)) return(NULL)  #Need to have either genotype of CV
    #print("FarmCPU.GLM Started")
    #print("Dimention of Y, GDP, GM, and CV")
    #print(dim(Y))
    #print(dim(GDP))
    #print(dim(GM))
    #print(dim(CV))
    #print(head(CV))
    #print("Solving equation (This may take a while)...")
    if(!is.null(CV)){
        CV=as.matrix(CV)
        nf=ncol(CV)
    }else{
        nf=0
    }
    print("number of covariates in current loop is:")
    print(nf)
    #Build model with SNP as the last variable
    #if(package!="FarmCPU.LM"){
    if(is.null(CV)) {
        myModel="Y [,2]~x"
        if(package!="fast.lm"){
            ccv=rep(1,nrow(Y))
        }
    }else{
        #CV=as.matrix(CV)
        seqCV=1:(ncol(CV))
        myModel=paste("Y[,2]~",paste("CV[,",(seqCV),"]",collapse= "+"),"+ x")
        #print(head(CV))
        #ccv=cbind(rep(1,nrow(Y)),as.matrix(CV[,2:ncol(CV)]))
        if(package!="fast.lm"){
            ccv=cbind(rep(1,nrow(Y)),as.matrix(CV))
        }
        #ccv=as.matrix(CV)
        #print("ccv")
        #print(head(ccv))
    }
    #}
    ##print("The model is: ")
    ##print(myModel)
    #===========================by lm=======================================
    if(package=="lm"){
        #Solve the  model with lm
        #P <- apply(GDP,2,function(x){
        #fmla <- formula(myModel)
        #myLM=lm(fmla)
        #lms=summary(myLM)
        #lmcoef=lms$coefficients
        #lmcoefOnly=lmcoef[-1,]  #remove intercept
        ##print(lmcoefOnly)
        #(as.numeric(lmcoefOnly))
        #In order of estimate, t, se and P
        #cbind(lmcoefOnly[,1],lmcoefOnly[,3],lmcoefOnly[,2],lmcoefOnly[,4])
        #})
        
        P<-matrix(NA,nrow=nrow(GDP),ncol=4*(nf+1))
        for(i in 1:nrow(GDP)){
            x <- GDP[i,]
            fmla <- stats::formula(myModel)
            myLM = stats::lm(fmla)
            lms=summary(myLM)
            lmcoef=lms$coefficients
            lmcoefOnly=lmcoef[-1,]  #remove intercept
            ##print(lmcoefOnly)
            #(as.numeric(lmcoefOnly))
            #In order of estimate, t, se and P
            #P[i,]=cbind(lmcoefOnly[,1],lmcoefOnly[,3],lmcoefOnly[,2],lmcoefOnly[,4])
            #P[i,1]=lmcoefOnly[1]
            #P[i,2]=lmcoefOnly[3]
            #P[i,3]=lmcoefOnly[2]
            #P[i,4]=lmcoefOnly[4]
            #print(lmcoefOnly)
            
            P[i,c(1:(nf+1))]=lmcoefOnly[1:(nf+1)]
            P[i,c((nf+2):(2*nf+2))]=lmcoefOnly[(2*nf+3):(3*nf+3)]
            P[i,c((2*nf+3):(3*nf+3))]=lmcoefOnly[(nf+2):(2*nf+2)]
            P[i,c((3*nf+4):(4*nf+4))]=lmcoefOnly[(3*nf+4):(4*nf+4)]
            
            #P[i,c(1:(nf+1))]=lmcoefOnly[1]
            #P[i,c((nf+2):(2*nf+2))]=lmcoefOnly[3]
            #P[i,c((2*nf+3):(3*nf+3))]=lmcoefOnly[2]
            #P[i,c((3*nf+4):(4*nf+4))]=lmcoefOnly[4]
            
            
        }
        #print(head(P))
        #convert list to numeric matrix
        #P=t(sapply(P, function(row, max_length) c(row, rep(NA, max_length - length(row))), max(sapply(P, length))))
        #the following two do not work
        #P=as.data.frame(do.call(rbind, P))
        #P=t(as.matrix(P))
        P0=NULL
        pred=NULL
        PF=P[,ncol(P)]
        myLM=list(P=P,P0=P0,PF=PF,Pred=pred)
    } #end of lm if statement
    #===========================by fast.lm=======================================
    if(package=="fast.lm"){
        #fast.lm does not all?ow missing values
        missing=is.na(Y[,2]) #index for missing phenotype
        Mtotal=ncol(GDP)
        Ym=Y[missing,]
        Y=Y[!missing,]
        ccv=ccv[!missing,]
        GDP=GDP[!missing,]
        
        #set index for markers with no variation
        varSNP=apply(GDP, 2, stats::var)
        indexSNP=which(varSNP!=0)
        
        P0 <- apply(GDP[,indexSNP],2,function(x){
            x = cbind(ccv,x)
            fast.lm = RcppArmadillo::fastLmPure(y=Y[,2],X = x)
            tvalue=fast.lm$coefficients[-1]/fast.lm$stderr[-1]
            pvalue = 2 * stats::pt(abs(tvalue), fast.lm$df.residual, lower.tail = FALSE)
            cbind(fast.lm$coefficients[-1],tvalue,fast.lm$stderr[-1],pvalue)
        })
        
        #convert list to numeric matrix, the last (t+1) columns are p values for SNPs
        P0=t(as.matrix(P0))
        #Restore in original oder
        mtotal=ncol(GDP)
        nfix=ncol(P0)
        P=matrix(NA,Mtotal,nfix)
        rownames(P)=colnames(GDP)  #This should be OK
        P[indexSNP,]=P0 #restore the order with markers without variation
        P0=NULL
        pred=NULL
        PF=P[,ncol(P)]
        myLM=list(P=P,P0=P0,PF=PF,Pred=pred)
    }# end of fast.lm if statement
    #===========================by FarmCPU.LM=======================================
    if(package=="FarmCPU.LM"){
        #print("Calling GLM")
        #print(date())
        #print("Memory used before calling LM")
        #print(memory.size())
        gc()
        theCV=NULL
        if(!is.null(CV)) {
            theCV=as.matrix(CV)#as.matrix(CV[,-1])#
            seqCV=1:(ncol(theCV))
            myModel=paste("y~",paste("w[,",(seqCV),"]",collapse= "+"),"+ x")
        }
        #print("theCV")
        #print(head(theCV))
        if(ncpus == 1){
          myLM = FarmCPU.LM(y = Y[,2],
                            w = theCV,
                            GDP = GDP, 
                            orientation = orientation,
                            model = model,
                            ncpus = ncpus,
                            myModel = myModel,
                            seqQTN = seqQTN,
                            npc = npc)
        }
        if(ncpus>1){myLM = FarmCPU.LM.Parallel(y=Y[,2],
                    w = theCV,
                    x = GDP,
                    orientation = orientation,
                    model = model,
                    ncpus = ncpus#,
                    #npc=npc
                    )
        }
        #print("Memory used after calling LM")
        #print(memory.size())
        gc()
        
    }# end of FarmCPU.lm if statement
    
    #print("FarmCPU.GLM accoplished")
    #print(date())
    gc()
    #return(list(P=myLM$P,P0=myLM$P0,PF=myLM$PF,Pred=myLM$pred))
    return(myLM)
}#The function FarmCPU.GLM ends here
`FarmCPU.Inv` <- function(A){
    #Object: To invert a 2 by 2 matrix quickly
    #intput: A -  2 by 2 matrix
    #Output: Inverse
    #Authors: Zhiwu Zhang
    # Last update: March 6, 2013
    ##############################################################################################
    detA=A[1,1]*A[2,2]-A[1,2]*A[2,1]
    temp=A[1,1]
    A=-A
    A[1,1]=A[2,2]
    A[2,2]=T
    return(A/detA)
}#The function FarmCPU.Inv ends here
`FarmCPU.LM.Parallel` <-
function(y,w=NULL,x,orientation="col",model="A",ncpus=2){
    #Object: 1. To quickly sovel LM with one variable substitute multiple times
    #Object: 2. To fit additive and additive+dominace model
    #intput: y - dependent variable
    #intput: w - independent variable
    #intput: x - independent variable of substitution (GDP)
    #intput: model - genetic effects. Options are "A" and "AD"
    #Output: estimate, tvalue, stderr and pvalue ( plus the P value of F test on both A and D)
    #Straitegy: 1. Separate constant covariates (w) and dynamic coveriates (x)
    #Straitegy: 2. Build non-x related only once
    #Straitegy: 3. Use apply to iterate x
    #Straitegy: 4. Derive dominance indicate d from additive indicate (x) mathmaticaly
    #Straitegy: 5. When d is not estimable, continue to test x
    #Authors: Xiaolei Liu and Zhiwu Zhang
    #Start  date: March 1, 2013
    #Last update: March 6, 2013
    ##############################################################################################
    print("FarmCPU.LM started")
    print(date())
    print(paste("No. Obs: ",length(y),sep=""))
    print("diminsion of covariates and markers")
    if(!is.null(w))print(dim(w))
    
    print("Memory used at begining of LM")
    if(.Platform$OS.type == "windows"){print(utils::memory.size())}
#    print(utils::memory.size())
    gc()
    #Constant section (non individual marker specific)
    #---------------------------------------------------------
    #Configration
    nd=20 #number of markes for checking A and D dependency
    threshold=.99 # not solving d if correlation between a and d is above this
    N=length(y) #Total number of taxa, including missing ones
    direction=2
    if(orientation=="row")direction=1
    print("direction")
    print(direction)
    #Handler of non numerical y a and w
    
    if(!is.null(w)){
        nf=length(w)/N
        w=matrix(as.numeric(as.matrix(w)),N,nf  )
        w=cbind(rep(1,N),w)#add overall mean indicator
        q0=ncol(w) #Number of fixed effect excluding gnetic effects
    }else{
        w=rep(1,N)
        nf=0
        q0=1
    }
    
    y=matrix(as.numeric(as.matrix(y)),N,1  )
    
    print("Adding overall mean")
    print(date())
    
    print("Build the static section")
    print(date())
    
    #n=nrow(w) #number of taxa without missing
    n=N
    if(nd>n)nd=n #handler of samples less than nd
    k=1 #number of genetic effect: 1 and 2 for A and AD respectively
    if(model=="AD")k=2
    
    q1=(q0+1) # vecter index for the posistion of genetic effect (a)
    q2=(q0+1):(q0+2) # vecter index for the posistion of genetic effect (a and d)
    df=n-q0-k #residual df (this should be varied based on validating d)
    
    iXX=matrix(0,q0+k,q0+k) #Reserve the maximum size of inverse of LHS
    #theNA=c(rep(NA,q0),rep(0,k)) # this should not be useful anymore
    
    ww=crossprod(w,w)
    wy=crossprod(w,y)
    yy=crossprod(y,y)
    # wwi=solve(ww) Revised by Jiabo on 2021.3.4
    wwi <- try(solve(ww),silent=TRUE)
     if(inherits(wwi, "try-error")){
      # print("!!!!!")
     wwi <- MASS::ginv(ww)
     }
    print("Prediction")
    print(date())
    
    #Statistics on the reduced model without marker
    rhs=wy
    beta <- crossprod(wwi,rhs)
    ve=(yy-crossprod(beta,rhs))/df
#    se=sqrt(diag(wwi)*ve)
    se=sqrt(diag(wwi) * as.vector(ve))
    tvalue=beta/se
    pvalue <- 2 * stats::pt(abs(tvalue), df,lower.tail = FALSE)
    P0=c(beta[-1],tvalue[-1],se[-1],pvalue[-1])
    yp=w%*%beta
    
    print("Detecting genotype coding system")
    print(date())
    
    #Finding the middle of genotype coding (1 for 0/1/2 and 0 for -1/0/1)
    s=5 # number of taxa sampled
    t0=which(x[1:s,]<0)
    t1=which(x[1:s,]>1)
    middle=0
    if(length(t0)<length(t1)) middle=1
    
    print("Memory used after setting LM")
    if(.Platform$OS.type == "windows"){print(utils::memory.size())}
#    print(utils::memory.size())
    gc()
    #Dynamic section on individual marker
    print("Iterating.................")
    print(date())
    print("dimension of GD")
    print(dim(x))
    print(methods::is(x))
    
    #sfInit(parallel=ncpus>1, cpus=ncpus)
    #print(sprintf('%s cpus are used', sfCpus()))
    
    #---------------------------------------------------------
    #P <- apply(x,direction,function(x){
    P <- snowfall::sfApply(x,direction,function(x){
        print("debug snowfall::sfApply")
        r=1 #initial creteria for correlation between a and d
        if(model=="AD"){
            d=1-abs(x-middle)
            r=abs(stats::cor(x[1:nd],d[1:nd]))
            if(is.na(r))r=1
            if(r<=threshold) x=cbind(x,d) # having both a and d as marker effects
        }
        print("make some noise here")
        #Process the edge (marker effects)
        xw=crossprod(w,x)
        xy=crossprod(x,y)
        xx=crossprod(x,x)
        
        B21 <- crossprod(xw, wwi)
        #t1=crossprod(xw,wwi)
        t2=B21%*%xw #I have problem of using crossprod and tcrossprod here
        B22 <- xx - t2
        
        #B22 can a scaler (A model) or 2 by2 matrix (AD model)
        if(model=="AD"&r<=threshold){
            invB22 <- FarmCPU.Inv(B22)
        }else{
            invB22=1/B22
        }
        
        NeginvB22B21 <- crossprod(-invB22,B21)
        
        if(model=="AD"&r<=threshold){
            B11 <- wwi + crossprod(B21,B21)
        }else{
            B11 <- wwi + as.numeric(invB22)*crossprod(B21,B21)
        }
        
        #Derive inverse of LHS with partationed matrix
        iXX[1:q0,1:q0]=B11
        
        if(r>threshold){
            iXX[q1,q1]=invB22
            iXX[q1,1:q0]=NeginvB22B21
            iXX[1:q0,q1]=NeginvB22B21
        }else{
            iXX[q2,q2]=invB22
            iXX[q2,1:q0]=NeginvB22B21
            iXX[1:q0,q2]=NeginvB22B21
        }
        
        #statistics
        rhs=c(wy,xy) #the size varied automaticly by A/AD model and validated d
        
        if(abs(r)>threshold & model=="AD"){
            beta <- crossprod(iXX[-(q0+k),-(q0+k)],rhs) #the last one (d) dose not count
            df=n-q0-1
        }else{
            beta <- crossprod(iXX,rhs)   #both a and d go in
            df=n-q0-2
        }
        if(model=="A") df=n-q0-1 #change it back for model A
        
        ve=(yy-crossprod(beta,rhs))/df #this is a scaler
        
        #using iXX in the same as above to derive se
        if(abs(r)>threshold & model=="AD"){
            se=sqrt(diag(iXX[-(q0+k),-(q0+k)])*ve)
            
        }else{
            #se=sqrt(diag(iXX)*ve)
            se = sqrt(diag(iXX) * c(ve))
        }
        
        tvalue=beta/se
        pvalue <- 2 * stats::pt(abs(tvalue), df,lower.tail = FALSE)
        
        #Handler of dependency between  marker are covariate
        #if(abs(B22[1,1])<10e-8)pvalue[]=NA
        
        #Calculate P value for A+D effect
        if(model=="AD"){
            #the last bit could be d or a, the second last may be marker effect not even not
            #In either case, calculate F and P value and correct them later
            markerbits=(length(beta)-1):length(beta)
            SSM=crossprod(beta[markerbits],rhs[markerbits])
            F=(SSM/2)/ve
            PF=df(F,2,df)
            
            #correcting PF with P from t value
            if(r>threshold) PF=pvalue[length(pvalue)]
        }
        
        #in case AD model and a/d dependent, add NA column at end
        if(r>threshold & model=="AD"){
            beta=c(beta,NA)
            tvalue=c(tvalue,NA)
            se=c(se,NA)
            pvalue=c(pvalue,NA)
        }
        
        if(model=="AD"){
            result=c(beta[-1],tvalue[-1],se[-1],pvalue[-1],PF)
        }else{
            result=c(beta[-1],tvalue[-1],se[-1],pvalue[-1])
        }
    }) #end of defyning apply function
    #sfStop()
    
    print("iteration accoplished")
    print(date())
    print("Memory used after iteration")
    if(.Platform$OS.type == "windows"){print(utils::memory.size())}
#    print(utils::memory.size())
    gc()
    
    #Final report
    #---------------------------------------------------------
    P=t(as.matrix(P))
    
    PF=P[,ncol(P)]
    if(model=="AD")P=P[,-ncol(P)]
    
    print("FarmCPU.LM accoplished")
    print(date())
    
    
    print("Memory used at end of LM")
    if(.Platform$OS.type == "windows"){print(utils::memory.size())}
#    print(utils::memory.size())
    gc()
    
    return(list(P=P,P0=P0,PF=PF,Pred=yp))
}
#)#end of cmpfun(
`FarmCPU.LM` <-
#cmpfun(
function(y,w=NULL,GDP,orientation="col",model="A",ncpus=2,myModel=NULL,seqQTN=NULL,npc=0){
    #Object: 1. To quickly sovel LM with one variable substitute multiple times
    #Object: 2. To fit additive and additive+dominace model
    #intput: y - dependent variable
    #intput: w - independent variable
    #intput: GDP - independent variable of substitution (GDP)
    #intput: model - genetic effects. Options are "A" and "AD"
    #Output: estimate, tvalue, stderr and pvalue ( plus the P value of F test on both A and D)
    #Straitegy: 1. Separate constant covariates (w) and dynamic coveriates (x)
    #Straitegy: 2. Build non-x related only once
    #Straitegy: 3. Use apply to iterate x
    #Straitegy: 4. Derive dominance indicate d from additive indicate (x) mathmaticaly
    #Straitegy: 5. When d is not estimable, continue to test x
    #Authors: Xiaolei Liu and Zhiwu Zhang
    #Start  date: March 1, 2013
    #Last update: March 6, 2013
    ##############################################################################################
    #print("FarmCPU.LM started")
    #print(date())
    #print(paste("No. Obs: ",length(y),sep=""))
    #print("diminsion of covariates and markers")
    if(!is.null(w))#print(dim(w))
    
    #print("Memory used at begining of LM")
    #print(memory.size())
    gc()
    #Constant section (non individual marker specific)
    #---------------------------------------------------------
    #Configration
    nd=20 #number of markes for checking A and D dependency
    threshold=.99 # not solving d if correlation between a and d is above this
    N=length(y) #Total number of taxa, including missing ones
    direction=2
    if(orientation=="row")direction=1
    #print("direction")
    #print(direction)
    #Handler of non numerical y a and w
    
    if(!is.null(w)){
        nf=length(w)/N
        w=matrix(as.numeric(as.matrix(w)),N,nf  )
        w=cbind(rep(1,N),w)#add overall mean indicator
        q0=ncol(w) #Number of fixed effect excluding gnetic effects
    }else{
        w=rep(1,N)
        nf=0
        q0=1
    }
    
    y=matrix(as.numeric(as.matrix(y)),N,1  )
    
    #print("Adding overall mean")
    #print(date())
    #print("Build the static section")
    #print(date())
    
    #n=nrow(w) #number of taxa without missing
    n=N
    if(nd>n)nd=n #handler of samples less than nd
    k=1 #number of genetic effect: 1 and 2 for A and AD respectively
    if(model=="AD")k=2
    
    q1=(q0+1) # vecter index for the posistion of genetic effect (a)
    q2=(q0+1):(q0+2) # vecter index for the posistion of genetic effect (a and d)
    df=n-q0-k #residual df (this should be varied based on validating d)
    
    iXX=matrix(0,q0+k,q0+k) #Reserve the maximum size of inverse of LHS
    #theNA=c(rep(NA,q0),rep(0,k)) # this should not be useful anymore
    
    ww=crossprod(w,w)
    wy=crossprod(w,y)
    yy=crossprod(y,y)
    # wwi=solve(ww) Revised by Jiabo on 2021.3.4
    wwi <- try(solve(ww),silent=TRUE)
     if(inherits(wwi, "try-error")){
      # print("!!!!!")
     wwi <- MASS::ginv(ww)
     }
    #print("Prediction")
    #print(date())
    
    #Statistics on the reduced model without marker
    rhs=wy
    beta <- crossprod(wwi,rhs)
    ve=(yy-crossprod(beta,rhs))/df
#    se=sqrt(diag(wwi)*ve)
    se=sqrt(diag(wwi) * as.vector(ve))
    tvalue=beta/se
    pvalue <- 2 * stats::pt(abs(tvalue), df,lower.tail = FALSE)
    P0=c(beta[-1],tvalue[-1],se[-1],pvalue[-1])
    yp=w%*%beta
    
    if(npc!=0){
        betapc = beta[2:(npc+1)]
        betapred = beta[-c(1:(npc+1))]
    }else{
        betapc = NULL
        betapred = beta[-1]
    }
    #print("Detecting genotype coding system")
    #print(date())
    
    #Finding the middle of genotype coding (1 for 0/1/2 and 0 for -1/0/1)
    s=5 # number of taxa sampled
    t0=which(GDP[1:s,]<0)
    t1=which(GDP[1:s,]>1)
    middle=0
    if(length(t0)<length(t1)) middle=1
    
    #print("Memory used after setting LM")
    #print(memory.size())
    gc()
    #Dynamic section on individual marker
    #print("Iterating.................")
    #print(date())
    #print("dimension of GD")
    #print(dim(x))
    #sfInit(parallel=ncpus>1, cpus=ncpus)
    ##print(sprintf('%s cpus are used', sfCpus()))
    
    #---------------------------------------------------------
    #P <- matrix(NA,nrow=nrow(GDP),ncol=4*(nf+1))
    if(orientation=="row"){
        P <- matrix(NA,nrow=nrow(GDP),ncol=nf+1)
        ntest=nrow(GDP)
    }else{
        P <- matrix(NA,nrow=ncol(GDP),ncol=nf+1)
        ntest=ncol(GDP)
    }
    
    if(orientation=="row"){
        B <- matrix(NA,nrow=nrow(GDP),ncol=1)
    }else{
        B <- matrix(NA,nrow=ncol(GDP),ncol=1)
    }
    
    for(i in 1:ntest){
        if(orientation=="row"){
            x=GDP[i,]
        }else{
            x=GDP[,i]
        }
        
        #P <- apply(x,direction,function(x){
        #P <- sfApply(x,direction,function(x){
        r=1 #initial creteria for correlation between a and d
        if(model=="AD"){
            d=1-abs(x-middle)
            r=abs(stats::cor(x[1:nd],d[1:nd]))
            if(is.na(r))r=1
            if(r<=threshold) x=cbind(x,d) # having both a and d as marker effects
        }
        
        #Process the edge (marker effects)
        xy=crossprod(x,y)
        xx=crossprod(x,x)
        
        if(model=="AD"&r<=threshold){
            xw=crossprod(x,w)
            wx=crossprod(w,x)
            iXX22 <- solve(xx-xw%*%wwi%*%wx)
            iXX12 <- (-wwi)%*%wx%*%iXX22
            iXX21 <- (-iXX22)%*%xw%*%wwi
            iXX11 <- wwi + wwi%*%wx%*%iXX22%*%xw%*%wwi
        }else{
            xw=crossprod(w,x)
            B21 <- crossprod(xw, wwi)
            t2=B21%*%xw #I have problem of using crossprod and tcrossprod here
            B22 <- xx - t2
            invB22=1/B22
            NeginvB22B21 <- crossprod(-invB22,B21)
            iXX11 <- wwi + as.numeric(invB22)*crossprod(B21,B21)
        }
        
        #Derive inverse of LHS with partationed matrix
        iXX[1:q0,1:q0]=iXX11
        
        if(r>threshold){
            iXX[q1,q1]=invB22
            iXX[q1,1:q0]=NeginvB22B21
            iXX[1:q0,q1]=NeginvB22B21
        }else{
            iXX[q2,q2]=iXX22
            iXX[q2,1:q0]=iXX21
            iXX[1:q0,q2]=iXX12
        }
        
        #statistics
        rhs=c(wy,xy) #the size varied automaticly by A/AD model and validated d
        
        if(abs(r)>threshold & model=="AD"){
            beta <- crossprod(iXX[-(q0+k),-(q0+k)],rhs) #the last one (d) dose not count
            df=n-q0-1
        }else{
            beta <- crossprod(iXX,rhs)   #both a and d go in
            df=n-q0-2
        }
        if(model=="A") df=n-q0-1 #change it back for model A
        
        ve=(yy-crossprod(beta,rhs))/df #this is a scaler
        
        #using iXX in the same as above to derive se
        if(abs(r)>threshold & model=="AD"){
            se=sqrt(diag(iXX[-(q0+k),-(q0+k)])*ve)
        }else{
            # browser()
            # se=sqrt(diag(iXX)*ve)
            # se = sqrt(diag(iXX) * c(ve))
            myDiag <- diag(iXX)
            myDiag[ myDiag < 0 ] <- 0
            ve[ ve < 0 ] <- 0
            se = sqrt(myDiag * c(ve))
        }
        
        tvalue=beta/se
        pvalue <- 2 * stats::pt(abs(tvalue), df,lower.tail = FALSE)
        
        #Handler of dependency between  marker are covariate
        if(!is.na(abs(B22[1,1]))){
            if(abs(B22[1,1])<10e-8)pvalue[]=NA}
        
        #Calculate P value for A+D effect
        if(model=="AD"){
            #the last bit could be d or a, the second last may be marker effect not even not
            #In either case, calculate F and P value and correct them later
            markerbits=(length(beta)-1):length(beta)
            SSM=crossprod(beta[markerbits],rhs[markerbits])
            F=(SSM/2)/ve
            PF=df(F,2,df)
            
            #correcting PF with P from t value
            if(r>threshold) PF=pvalue[length(pvalue)]
        }
        
        #in case AD model and a/d dependent, add NA column at end
        if(r>threshold & model=="AD"){
            beta=c(beta,NA)
            tvalue=c(tvalue,NA)
            se=c(se,NA)
            pvalue=c(pvalue,NA)
        }
        
        if(model=="AD"){
            result=c(beta[-1],tvalue[-1],se[-1],pvalue[-1],PF)
        }else{
            #result=c(beta[-1],tvalue[-1],se[-1],pvalue[-1])
            #P[i,]=c(beta[-1],tvalue[-1],se[-1],pvalue[-1])
            P[i,c(1:(nf+1))]=pvalue[-1]
            B[i,]=beta[length(beta)]
            #P[i,c(1:(nf+1))]=beta[-1]
            #P[i,c((nf+2):(2*nf+2))]=pvalue[-1]
            #P[i,c((nf+2):(2*nf+2))]=tvalue[-1]
            #P[i,c((2*nf+3):(3*nf+3))]=se[-1]
            #P[i,c((3*nf+4):(4*nf+4))]=pvalue[-1]
        }
    }
    #}
    #}) #end of defyning apply function
    #sfStop()
    
    #print("iteration accoplished")
    #print(date())
    #print("Memory used after iteration")
    #print(memory.size())
    gc()
    
    #Final report
    #---------------------------------------------------------
    #P=t(as.matrix(P))
    #P=as.matrix(P)
    
    PF=P[,ncol(P)]
    if(model=="AD")P=P[,-ncol(P)]
    
    #print("FarmCPU.LM accoplished")
    #print(date())
    
    #print(dim(P))
    #print(P[1:5,])
    #print("Memory used at end of LM")
    #print(memory.size())
    gc()
    #print(head(P))
    return(list(P=P,P0=P0,PF=PF,Pred=yp,betapc=betapc,betapred=betapred,B=B))
} #end of function(
#)#end of cmpfun(
`FarmCPU.Pred` <- function(pred=NULL,ypred=NULL,name.of.trait=""){
    #Object: To display the correlation between observed phenotype and predicted phenotype
    #Input 1: pred, the first column is taxa name, the second column is observed phenotype and the third column is predicted phenotype
    #Input 2: ypred, the first column is taxa name, the second column is observed phenotype and the third column is predicted phenotype, the different between pred and ypred is that pred is to predict phenotypes with observed values already, ypred is to predict phenotype that is NA
    #Output: cor:correlation between observed phenotype and real phenotype (comment: pred is to predict phenotypes with observed values already)
    #Output: ycor:correlation between observed phenotype and real phenotype (comment: ypred is to predict phenotype that is NA)
    #Output: A table and plot (pdf)
    #Requirment: NA
    #Authors: Xiaolei Liu
    #Start date: June 26, 2014
    #Last update: June 26, 2014
    ##############################################################################################
    #print("Create prediction table..." )
    cor=NA
    ycor=NA
    if(!is.null(pred)) {
        index=!is.na(pred[,2])
        utils::write.table(pred, paste("FarmCPU.", name.of.trait, ".Pred.csv", sep = ""), quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
        #pred=read.table("FarmCPU.Iteration_02.Farm-CPU.Sim1.Pred.csv",sep=",",header=T)
        
        grDevices::pdf(paste("FarmCPU.", name.of.trait,".Pred.pdf" ,sep = ""), width = 5,height=5)
        graphics::par(mar = c(5,6,5,3))
        pred.lm = stats::lm(pred[,3][index]~pred[,2][index])
        plot(pred[,3][index]~pred[,2][index],pch=20,col='black',ylab="Predicted phenotype",xlab="Observed phenotype",cex.axis=1,cex=1,cex.lab=1,las=1,bty='n',xlim=c(floor(min(pred[,2],na.rm=T)),ceiling(max(pred[,2],na.rm=T))*1.2),ylim=c(floor(min(pred[,3],na.rm=T)),ceiling(max(pred[,3],na.rm=T))*1.2),xaxs="i",yaxs="i")
        graphics::abline(pred.lm,lty=5,col='red',lwd=2)
        #legend(max(pred[,3])+1,max(pred[,2])+1, paste("R^2 = ", 0.5), col = 'black', text.col = "black", lty = 1, ncol=1, cex = 1, lwd=2, bty='o')
        cor=round(summary(pred.lm)$r.sq, 3)
        graphics::text(max(pred[,2],na.rm=T)*1, max(pred[,3],na.rm=T)*1, paste("R^2=", cor), col= "forestgreen", cex = 1, pos=3)
        #title(paste("R^2 = ", round(summary(pred.lm)$r.sq, 3)), col= "black", cex = 1)
        grDevices::dev.off()
    }
    #print("Create prediction table for unknown phenotype...")
    if(!is.null(ypred)){
        yindex=!is.na(ypred[,2])
        ypredrna=ypred[,2][yindex]
        utils::write.table(ypred, paste("FarmCPU.", name.of.trait, ".unknownY.Pred.csv", sep = ""), quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
        if(length(ypredrna)!=0){
            grDevices::pdf(paste("FarmCPU.", name.of.trait,".unknownY.Pred.pdf" ,sep = ""), width = 5,height=5)
            graphics::par(mar = c(5,6,5,3))
            ypred.lm = stats::lm(ypred[,3][yindex]~ypredrna)
            plot(ypred[,3][yindex]~ypredrna,pch=20,col='black',ylab="Predicted phenotype",xlab="Observed phenotype",cex.axis=1,cex=1,cex.lab=1,las=1,bty='n',xlim=c(floor(min(pred[,2],na.rm=T)),ceiling(max(ypred[,2],na.rm=T))*1.2),ylim=c(floor(min(pred[,3],na.rm=T)),ceiling(max(ypred[,3],na.rm=T))*1.2),xaxs="i",yaxs="i")
            graphics::abline(ypred.lm,lty=5,col='red',lwd=2)
            ycor=round(summary(ypred.lm)$r.sq, 3)
            graphics::text(max(ypred[,2],na.rm=T)*1,max(ypred[,3],na.rm=T)*1, paste("R^2=", ycor), col= "forestgreen", cex = 1, pos=3)
            grDevices::dev.off()
        }else{
            print("There is no observed phenotype for predicted phenotype")
        }
    }
    return(list(cor=cor,ycor=ycor))
}#end of `FarmCPU.Pred`
`FarmCPU.Prior` <-
function(GM,P=NULL,Prior=NULL,kinship.algorithm="FARM-CPU"){
    #Object: Set prior on existing p value
    #Input: GM - m by 3  matrix for SNP name, chromosome and BP
    #Input: Prior - s by 4  matrix for SNP name, chromosome, BP and Pvalue
    #Input: P - m by 1 matrix containing probability
    #Requirement: P and GM are in the same order, Prior is part of GM except P value
    #Output: P - m by 1 matrix containing probability
    #Authors: Zhiwu Zhang
    # Last update: March 10, 2013
    ##############################################################################
    #print("FarmCPU.Prior Started")
    #print("dimension of GM")
    #print(dim(GM))
    
    if(is.null(Prior)& kinship.algorithm!="FARM-CPU")return(P)
    if(is.null(Prior)& is.null(P))return(P)
    
    #get prior position
    if(!is.null(Prior)) index=match(Prior[,1],GM[,1],nomatch = 0)
    
    #if(is.null(P)) P=runif(nrow(GM)) #set random p value if not provided (This is not helpful)
    #print("debug set prior  a")
    
    #Get product with prior if provided
    if(!is.null(Prior) & !is.null(P)  )P[index]=P[index]*Prior[,4]
    
    #print("debug set prior   b")
    return(P)
}#The function FarmCPU.Prior ends here

#'
#'
#' FarmCPU
#' 
#' @description 
#' FarmCPU: GWAS and GS by using FarmCPU method
#' 
#' 
#' @param Y = NULL, a data.frame of phenotype data, first column is sample name, second column is the trait.
#' @param GD = NULL,
#' @param GM = NULL,
#' @param CV = NULL,
#' @param GP = NULL,
#' @param Yt = NULL,
#' @param DPP = 1000000,
#' @param kinship.algorithm = "FARM-CPU",
#' @param file.output = TRUE,
#' @param cutOff = 0.01,
#' @param method.GLM = "FarmCPU.LM",
#' @param method.sub = "reward",
#' @param method.sub.final = "reward",
#' @param method.bin = "static",
#' @param bin.size = c(5e5,5e6,5e7),
#' @param bin.selection = seq(10,100,10),
#' @param memo = NULL,
#' @param Prior = NULL,
#' @param ncpus = 1,
#' @param maxLoop = 10,
#' @param threshold.output = .01,
#' @param WS = c(1e0,1e3,1e4,1e5,1e6,1e7),
#' @param alpha = c(.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),
#' @param maxOut = 100,
#' @param QTN.position = NULL,
#' @param converge = 1,
#' @param iteration.output = FALSE,
#' @param acceleration = 0,
#' @param model = "A",
#' @param MAF.calculate = FALSE,
#' @param plot.style = "FarmCPU",
#' @param p.threshold = NA,
#' @param QTN.threshold = 0.01,
#' @param maf.threshold = 0.03,
#' @param ycor = NULL,
#' @param bound = NULL
#' 
#' 
#' @return 
#' A list.
#' 
#' 
#' @author Xiaolei Liu and Zhiwu Zhang
#' 
#' 
#' @examples 
#' \dontrun{
#' myPhenoFile <- system.file("extdata", "mdp_traits.txt.gz", package = "GAPIT3")
#' myPhenotypes <- read.table(myPhenoFile, header = TRUE)
#' myFarmCPU <- FarmCPU(myPhenotypes[, 1:2])
#' }
#' 
#' 
#' @export
`FarmCPU` <- function(Y = NULL,
                      GD = NULL,
                      GM = NULL,
                      CV = NULL,
                      GP = NULL,
                      Yt = NULL,
                      DPP = 1000000,
                      kinship.algorithm = "FARM-CPU",
                      file.output = TRUE,
                      cutOff = 0.01,
                      method.GLM = "FarmCPU.LM",
                      method.sub = "reward",
                      method.sub.final = "reward",
                      method.bin = "static",
                      bin.size = c(5e5,5e6,5e7),
                      bin.selection = seq(10,100,10),
                      memo = NULL,
                      Prior = NULL,
                      ncpus = 1,
                      maxLoop = 10,
                      threshold.output = .01,
                      WS = c(1e0,1e3,1e4,1e5,1e6,1e7),
                      alpha = c(.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),
                      maxOut = 100,
                      QTN.position = NULL,
                      converge = 1,
                      iteration.output = FALSE,
                      acceleration = 0,
                      model = "A",
                      MAF.calculate = FALSE,
                      plot.style = "FarmCPU",
                      p.threshold = NA,
                      QTN.threshold = 0.01,
                      maf.threshold = 0.03,
                      ycor = NULL,
                      bound = NULL){
    #Object: GWAS and GS by using FarmCPU method
    #Input: Y,GD,GM,CV
    #Input: GD - n by m +1 dataframe or n by m big.matrix
    #Input: GD - n by m matrix. This is Genotype Data Pure (GD). THERE IS NOT COLUMN FOR TAXA.
    #Requirement: Y, GD and CV have same taxa order. GD and GM have the same order on SNP
    #Requirement: Y can have missing data. CV, GD and GM can not. Non-variable markers are allowed
    #Output: GWAS,GPS,Pred
    #Authors: Xiaolei Liu and Zhiwu Zhang
    # Date  start: Febuary 24, 2013
    # Last update: April 2, 2013
    ##############################################################################################
    #print("FarmCPU Started")
    #print(date())
    #print("Memory used at begining of BUS")
    #print(memory.size())
    #print(dim(GD))
    #print(dim(GM))
    print("--------------------- Welcome to FarmCPU ----------------------------")
    echo = TRUE
    FarmCPU.Version = FarmCPU.0000()
    print("FarmCPU Started...")
    if(ncol(Y)>2) stop("FarmCPU only accept single phenotype, please specify a column, like myY[,c(1,3)]")
    #Set orientation
    #Strategy: the number of rows in GD and GM are the same if GD has SNP as row
    nm = nrow(GM)
    ny = nrow(Y)
    ngd1 = nrow(GD)
    ngd2 = ncol(GD)
    if(!is.null(CV)){
        CV = as.matrix(CV)
        npc = ncol(CV)
    }else{
        npc = 0
    }
    ngd1 = abs(ngd1-nm)
    ngd2 = abs(ngd2-nm)
    orientation = "col"
    theSNP = 2
    ns = nrow(GD)
    if(min(ngd1,ngd2)==0){
        orientation = "row"
        theSNP = 1
        ns = ncol(GD)
    }
    
    #acceleration
    ac = NULL
    if(acceleration!=0) ac = rep(1.0,nm)
    
    #Handler of non numeric chr
    #GM[,2]=as.numeric(GM[,2])
    
    #Handler 0 bp
    index = which(GM[,3]==0 )
    if(length(index)>0){
        #print("Warning: there is 0 bp which was set to 1")
        #print(length(index))
        GM[index,3] = 1      #This is problematic
    }
    
    #handler of multiple CPU on big.matrix
    if(ncpus>1 & bigmemory::is.big.matrix(GD)){
        #print("Multiple CPUs are not avaiable for big.matrix. ")
        #print("The big.matrix will be converted to regular matrix which takes more memmory")
        #stop("Import the genotype as regula R matrix or set single CPU")
    }
    
    #print("number of CPU required")
    #print(ncpus)
    if(ncpus>1) snowfall::sfInit(parallel = ncpus>1, cpus = ncpus)
    
    P = GP
    
    if(!is.null(GP))P = GP[,4] #get the p value
    
    #print("maxLoop")
    #print(maxLoop)
    gc()
    #print(memory.size())
    #print(date())
    #print(is(GD))
    #print(dim(GD))
    
    #handler of GD with taxa column
    if(ncol(GD)>nm & orientation=="col"){
        #print("GD has taxa column")
        if(bigmemory::is.big.matrix(GD)){
            #retain as bi.matrix
            GD = bigmemory::deepcopy(GD,rows = 1:nrow(GD),cols = 2:ncol(GD))  #This cause problem with multi cpu
        }else{
            GD = as.matrix(GD[,-1])
        }
    }#end of if(ncol...
    
    #Change to regula matrix for multiple CPUs
    if(ncpus>1)  GD = as.matrix(GD)
    
    #print("after remove taxa in GD")
    gc()
    #print(memory.size())
    #print(date())
    #print(is(GD))
    #print(dim(GD))
    
    if(model=="A"){
        shift = 0
    }else if(model=="AD"){
        shift = 1
    }else {
        print("Please choose 'A' model or 'AD' model")
    }
    #print("bin.selection")
    #print(bin.selection)
    
    #calculating MAF
    if(MAF.calculate==FALSE){
        MAF = NA
    }else{
        MAF = apply(GD,theSNP,mean)
        MAF = matrix(MAF,nrow = 1)
        MAF = apply(MAF,2,function(x) min(1-x/2,x/2))
    }
    
    for (trait in 2: ncol(Y))  {
        name.of.trait = colnames(Y)[trait]
        #print(paste("Processing trait: ",name.of.trait,sep=""))
        if(!is.null(memo)) name.of.trait = paste(memo,".",name.of.trait,sep = "")
        
        #===============================================================================
        #handler of missing phenotype (keep raw Y,CV and GD)
        #print(date())
        #print("Memory used before processing missing")
        #print(memory.size())
        
        #index for missing phenotype
        index = 1:nm
        seqTaxa = which(!is.na(Y[,trait]))
        if(MAF.calculate==TRUE){
            if(is.na(maf.threshold)){
                if(length(seqTaxa)<=100) maf.threshold = 0.05
                #if(length(seqTaxa)>100&&length(seqTaxa)<=500) maf.threshold=0.01
                #if(length(seqTaxa)>300&&length(seqTaxa)<=500) maf.threshold=0.05
                #if(length(seqTaxa)>500&&length(seqTaxa)<=1000) maf.threshold=0.01
                if(length(seqTaxa)>100) maf.threshold = 0
            }else{
                maf.threshold = maf.threshold
            }
            mafindex = (1:nm)[MAF>=maf.threshold]
            MAF = MAF[mafindex]
            index = mafindex
            GM = GM[index,]
            nm = length(index)
        }
        #predict = !(length(seqTaxa)==nrow(Y))#judge whether there is NA in phenotype
        predict = !is.null(Yt)#judge whether there is two phenotypes
        PredictYt = NULL
        ypred = NULL
        #print(length(seqTaxa))
        #print(nrow(Y))
        #print("predict")
        #print(predict)
        Y1 = Y[seqTaxa,]
        #if(is.numeric(CV)){CV1=CV[seqTaxa]
        #}else{
        #    CV1=CV[seqTaxa,]}
        CV1 = CV[seqTaxa,]
        
        #print(head(CV1))
        if(length(seqTaxa)<1) stop("FarmCPU stoped as no data in Y")
        
        #print("Extract genotype for phenotyped taxa")
        #print(memory.size())
        #print(is(GD))
        #print(dim(GD))
        #print(length(seqTaxa))
        #print(length(index))
        
        #GD based on big.matrix and orientation
        if(orientation=="col"){
            if(bigmemory::is.big.matrix(GD)){
                GD1 = bigmemory::deepcopy(GD,rows = seqTaxa,cols = index)
            }else{
                GD1 = GD[seqTaxa,index]
            }
        }else{
            if(bigmemory::is.big.matrix(GD)){
                GD1 = bigmemory::deepcopy(GD,rows = index,cols = seqTaxa)
            }else{
                GD1 = GD[index,seqTaxa]
            }
        }# end of if orientation
        
        #prepare the data for predict NA in phenotype
        if(predict){
            seqTaxa2 = which(is.na(Y[,trait]))
            
            #seqTaxa2=which(is.na(Yt[,trait]))
            #Y2=Yt[seqTaxa2,]
            PredictYt = Yt[seqTaxa2,]
            if(is.numeric(CV)){CV2 = CV[seqTaxa2]
            }else{
                CV2 = CV[seqTaxa2,]}
            
            #GD based on big.matrix and orientation
            if(orientation=="col"){
                if(bigmemory::is.big.matrix(GD)){
                    GD2 = bigmemory::deepcopy(GD,rows = seqTaxa2,cols = index)
                }else{
                    GD2 = GD[seqTaxa2,index]
                }
            }else{
                if(bigmemory::is.big.matrix(GD)){
                    GD2 = bigmemory::deepcopy(GD,rows = index,cols = seqTaxa2)
                }else{
                    GD2 = GD[index,seqTaxa2]
                }
            }# end of if orientation
        }
        #print("dim(GD2)")
        #print(dim(GD2))
        #Step 1: preliminary screening
        #print(date())
        #print("Memory used before 1st GLM")
        #print(memory.size())
        
        theLoop = 0
        theConverge = 0
        seqQTN.save = c(0)
        seqQTN.pre = c(-1)
        
        isDone = FALSE
        name.of.trait2 = name.of.trait
        
        
        #while(theLoop<maxLoop & !converge ) {
        while(!isDone) {
            theLoop = theLoop+1
            print(paste("Current loop: ",theLoop," out of maximum of ", maxLoop, sep = ""))
            #print(date())
            
            spacer = "0"
            if(theLoop>9)spacer = ""
            if(iteration.output) name.of.trait2 = paste("Iteration_",spacer,theLoop,".",name.of.trait,sep = "")
            if(method.bin=="NONE")maxLoop = 1 #force to exit for GLM model
            
            #Step 2a: Set prior
            #print("Memory used before Prior")
            #print(memory.size())
            
            myPrior = FarmCPU.Prior(GM = GM,P = P,Prior = Prior,kinship.algorithm = kinship.algorithm)
            #Step 2b: Set bins
            
            #print(myPrior[1:5])
            
            #print("Memory used before Bin")
            #print(memory.size())
            #print(date())
            
            if(theLoop<=2){
                myBin = FarmCPU.BIN(Y = Y1[,c(1,trait)],GDP = GD1,GM = GM,CV = CV1,orientation = orientation,P = myPrior,method = method.bin,b = bin.size,s = bin.selection,theLoop = theLoop,bound = bound)
            }else{
                myBin = FarmCPU.BIN(Y = Y1[,c(1,trait)],GDP = GD1,GM = GM,CV = theCV,orientation = orientation,P = myPrior,method = method.bin,b = bin.size,s = bin.selection,theLoop = theLoop)
            }
            
            #Step 2c: Remove bin dependency
            #print(date())
            #print("Memory used before Remove")
            #print(memory.size())
            
            #Remove QTNs in LD
            seqQTN = myBin$seqQTN
            ve.save = myBin$ve.save
            vg.save = myBin$vg.save
            #print(seqQTN)
            #if(theLoop==2&&is.null(seqQTN)){maxLoop=2}#force to exit for GLM model while seqQTN=NULL and h2=0
            if(theLoop==2){
                #print(head(P))
                #print(min(P,na.rm=TRUE))
                if(!is.na(p.threshold)){
                    if(min(myPrior,na.rm = TRUE)>p.threshold){
                        seqQTN = NULL
                        print("Top snps have little effect, set seqQTN to NULL!")
                        #print("**********FarmCPU ACCOMPLISHED**********")
                    }
                }else{
                    if(min(myPrior,na.rm = TRUE)>0.01/nm){
                        seqQTN = NULL
                        print("Top snps have little effect, set seqQTN to NULL!")
                        #print("**********FarmCPU ACCOMPLISHED**********")
                    }
                }
            }
            
            #when FARM-CPU can not work, make a new QQ plot and manhatthan plot
            if(theLoop==2&&is.null(seqQTN)){
                #Report
                GWAS = cbind(GM,P,MAF,myGLM$B)
                #if(isDone | iteration.output){
                gc()
                pred = myGLM$Pred
                #print(pred)
                if(!is.null(pred)) pred = cbind(Y1,myGLM$Pred) #Need to be consistant to CMLM
                #print(pred)
                p.GLM = GWAS[,4]
                p.GLM.log = -log10(stats::quantile(p.GLM,na.rm = TRUE,0.05))
                #set.seed(666)
                #bonf.log=-log10(quantile(runif(nm),0.05))
                bonf.log = 1.3
                bonf.compare = p.GLM.log/bonf.log
                p.FARMCPU.log = -log10(p.GLM)/bonf.compare
                GWAS[,4] = 10^(-p.FARMCPU.log)
                GWAS[,4][which(GWAS[,4]>1)] = 1
                #colnames(GWAS)=c(colnames(GM),"P.value","maf","nobs","Rsquare.of.Model.without.SNP","Rsquare.of.Model.with.SNP","FDR_Adjusted_P-values")
                colnames(GWAS) = c(colnames(GM),"P.value","maf","effect")
                
                Vp = stats::var(Y1[,2],na.rm = TRUE)
                
                #print("Calling Report..")
                if(file.output){
                    if(npc!=0){
                        betapc = cbind(c(1:npc),myGLM$betapc)
                        colnames(betapc) = c("CV","Effect")
                        utils::write.csv(betapc,paste("FarmCPU.",name.of.trait2,".CVeffect.csv",sep = ""),quote = F,row.names = FALSE)
                    }
                    GAPIT.Report(name.of.trait = name.of.trait2,GWAS = GWAS,pred = NULL,ypred = ypred,tvalue = NULL,stderr = stderr,Vp = Vp,DPP = DPP,cutOff = cutOff,threshold.output = threshold.output,MAF = MAF,seqQTN = QTN.position,MAF.calculate = MAF.calculate,plot.style = plot.style)
                    myPower = GAPIT.Power(WS = WS, alpha = alpha, maxOut = maxOut,seqQTN = QTN.position,GM = GM,GWAS = GWAS,MaxBP = 1e10)
                }
                #} #enf of is done
                break
            }#force to exit for GLM model while seqQTN=NULL and h2=0
            
            #print("debug seqQTN")
            #print(seqQTN)
            #print(seqQTN.save)
            if(!is.null(seqQTN.save)&&theLoop>1){
              #browser()
              # if(seqQTN.save!=0 & seqQTN.save!=-1 & !is.null(seqQTN) ) seqQTN = union(seqQTN,seqQTN.save) #Force previous QTNs in the model
              if(all(seqQTN.save!=0 & seqQTN.save!=-1 & !is.null(seqQTN) ) ) seqQTN = union(seqQTN,seqQTN.save) #Force previous QTNs in the model
                      #print("**********POSSIBLE QTNs combined**********")
            }
            #if(!is.null(seqQTN.save)){
            #if(theLoop>=4 && !is.null(seqQTN.save) && (length(intersect(seqQTN.pre,seqQTN))/length(union(seqQTN.pre,seqQTN)))==1){
            #if(seqQTN.save!=0 & seqQTN.save!=-1 & !is.null(seqQTN) )
            #{seqQTN=union(seqQTN,seqQTN.save) #Force previous QTNs in the model
            #}
            if(theLoop!=1){
                seqQTN.p = myPrior[seqQTN]
                if(theLoop==2){
                    #index.p=seqQTN.p<0.01/nm
                    index.p = seqQTN.p<QTN.threshold
                    if(!is.na(p.threshold)){
                        #index.p=seqQTN.p<p.threshold
                        index.p = seqQTN.p<QTN.threshold
                    }
                    seqQTN.p = seqQTN.p[index.p]
                    seqQTN = seqQTN[index.p]
                    seqQTN.p = seqQTN.p[!is.na(seqQTN)]
                    seqQTN = seqQTN[!is.na(seqQTN)]
                }else{
                    #print("seqQTN.save")
                    #print(seqQTN.save)
                    #print("seqQTN")
                    #print(seqQTN)
                    #print(length(seqQTN.save))
                    #print(seqQTN.p[1:length(seqQTN.save)]<1)
                    #print(seqQTN.p[(length(seqQTN.save)+1):length(seqQTN)]<0.01/nm)
                    
                    #index.p=seqQTN.p<(0.01/nm)
                    index.p = seqQTN.p<QTN.threshold
                    if(!is.na(p.threshold)){
                        #index.p=seqQTN.p<p.threshold
                        index.p = seqQTN.p<QTN.threshold
                    }
                    index.p[seqQTN%in%seqQTN.save]=TRUE
                    #print(index.p)
                    seqQTN.p=seqQTN.p[index.p]
                    seqQTN=seqQTN[index.p]
                    seqQTN.p=seqQTN.p[!is.na(seqQTN)]
                    seqQTN=seqQTN[!is.na(seqQTN)]
                }
            }
            
            myRemove=FarmCPU.Remove(GDP=GD1,GM=GM,seqQTN=seqQTN,seqQTN.p=seqQTN.p,orientation=orientation,threshold=.7)
            
            #Recoding QTNs history
            seqQTN=myRemove$seqQTN
            
            #if(length(setdiff(seqQTN,seqQTN.save))==0 & length(intersect(seqQTN,seqQTN.save))>0   ) converge=TRUE
            theConverge=length(intersect(seqQTN,seqQTN.save))/length(union(seqQTN,seqQTN.save))
            circle=(length(union(seqQTN,seqQTN.pre))==length(intersect(seqQTN,seqQTN.pre))  )
            
            #handler of initial status
            if(is.null(seqQTN.pre)){circle=FALSE
            }else{
                if(seqQTN.pre[1]==0) circle=FALSE
                if(seqQTN.pre[1]==-1) circle=FALSE
            }
            
            #print("circle objective")
            print("seqQTN")
            print(seqQTN)
            print("scanning...")
            if(theLoop==maxLoop){
                print(paste("Total number of possible QTNs in the model is: ", length(seqQTN),sep=""))
            }
            #print(seqQTN.save)
            #print(seqQTN.pre)
            #print(circle)
            
            #print(converge)
            #print("converge current")
            #print(theConverge)
            
            isDone=((theLoop>=maxLoop) | (theConverge>=converge)  |circle )
            
            seqQTN.pre=seqQTN.save
            seqQTN.save=seqQTN
            
            #myRemove=FarmCPU.Remove(GD=GD1,GM=GM,seqQTN=seqQTN,orientation=orientation,threshold=.7)
            #Step 3: Screen with bins
            rm(myBin)
            gc()
            #print(date())
            #print("Memory used before 2nd GLM")
            #print(memory.size())
            
            theCV=CV1
            if(!is.null(myRemove$bin)){
                if(theLoop==1){
                    theCV=cbind(CV1,myRemove$bin)
                }else{
                    #print("remove PCs since 2nd iteration")
                    theCV=cbind(CV1,myRemove$bin)
                    #theCV=myRemove$bin
                }
            }
            myGLM=FarmCPU.GLM(Y = Y1[,c(1,trait)],
                              GDP = GD1,
                              GM = GM,
                              CV = theCV,
                              orientation = orientation,
                              package = method.GLM,
                              ncpus = ncpus,
                              model = model,
                              seqQTN = seqQTN,
                              npc = npc)
            #Step 4: Background unit substitution
            #print(date())
            #print("Memory used before SUB")
            #print(memory.size())
            
            #print("After calling SUB")
            #How about having reward during the process and mean at end?
            if(!isDone){
                myGLM=FarmCPU.SUB(GM=GM,GLM=myGLM,QTN=GM[myRemove$seqQTN,],method=method.sub,model=model)
            }else{
                myGLM=FarmCPU.SUB(GM=GM,GLM=myGLM,QTN=GM[myRemove$seqQTN,],method=method.sub.final,model=model)
            }
            #print(date())
            P=myGLM$P[,ncol(myGLM$P)-shift]
            
            #acceleration
            if(!is.null(ac)){
                # ac = FarmCPU.Accelerate(ac = ac, 
                #                         QTN = myRemove$seqQTN, 
                #                         acceleration = acceleration)
                # The function 'FarmCPU.Accelerate()' does not exist.
                P=P/ac
            }
            #print("Acceleration in bus")
            index=which(ac>1)
            #print(cbind(index,ac[index],P[index]))
            #if P value is 0
            #if(min(P,na.rm=TRUE)==0) break
            P[P==0] <- min(P[P!=0],na.rm=TRUE)*0.01
            #Report
            if(isDone | iteration.output){
                #print("Report assemmbling...")
                #-------------------------------------------------------------------------------
                #Assemble result for report
                gc()
                pred=myGLM$Pred
                PredictY=NULL
                if(!is.null(theCV)&&predict){
                    #Statistics on the reduced model without marker
                    beta <- myGLM$betapred
                    #w=seqQTN
                    if(orientation=="row"){
                        predw=rbind(1,t(CV1),GD2[seqQTN,])
                    }else{
                        predw=cbind(1,CV1,GD2[,seqQTN])
                    }
                    #ypred=predw%*%beta
                    #if(!is.null(theCV)){
                    #nf=length(theCV)/length(seqTaxa2)
                    #theCV=matrix(as.numeric(as.matrix(theCV)),length(seqTaxa2),nf)
                    #predw=cbind(rep(1,length(seqTaxa2)),theCV)#add overall mean indicator
                    #}else{
                    #predw=rep(1,length(seqTaxa2))
                    #}
                    #print(dim(predw))
                    #print(predw)
                    #print(length(beta))
                    #print(beta)
                    PredictY=predw%*%beta
                    #print(PredictY)
                    #PredictYt[seqTaxa2,]=PredictY
                }
                if(!is.null(pred)) pred=cbind(Y1,myGLM$Pred) #Need to be consistant to CMLM
                if(!is.null(PredictY)) ypred=cbind(PredictYt,PredictY) #Need to be consistant to CMLM
                #P=myGLM$P[,ncol(myGLM$P)-shift]
                #myGLM$P is in order of estimate, tvalue, stderr and pvalue
                #nf=ncol(myGLM$P)/4
                #tvalue=myGLM$P[,nf*2-shift]
                #stderr=myGLM$P[,3*nf-shift]
                #print("MAF might cause problem")
                #print(length(MAF))
                #GWAS=cbind(GM,P,MAF,NA,NA,NA,NA)
                GWAS=cbind(GM,P,MAF,myGLM$B)
                #colnames(GWAS)=c(colnames(GM),"P.value","maf","nobs","Rsquare.of.Model.without.SNP","Rsquare.of.Model.with.SNP","FDR_Adjusted_P-values")
                colnames(GWAS)=c(colnames(GM),"P.value","maf","effect")
                Vp = stats::var(Y1[,2],na.rm=TRUE)
                
                if(!is.null(ypred)){
                    yindex=!is.na(ypred[,2])
                    ypredrna=ypred[,2][yindex]
                    ypred.lm = stats::lm(ypred[,3][yindex]~ypredrna)
                    ycor=round(summary(ypred.lm)$r.sq, 3)
                    #print(ycor)
                }
                
                
                #print("Calling Report..")
                if(file.output){
                    if(theLoop==1&&is.null(CV)){
                        
                        if(npc!=0){
                            betapc=cbind(c(1:npc),myGLM$betapc)
                            colnames(betapc)=c("CV","Effect")
                            utils::write.csv(betapc,paste("FarmCPU.",name.of.trait2,".CVeffect.csv",sep=""),quote=F,row.names=FALSE)
                        }
                        
                        GAPIT.Report(name.of.trait=name.of.trait2,GWAS=GWAS,pred=NULL,ypred=NULL,tvalue=NULL,stderr=stderr,Vp=Vp,DPP=DPP,cutOff=cutOff,threshold.output=threshold.output,MAF=MAF,seqQTN=QTN.position,MAF.calculate=MAF.calculate,plot.style=plot.style)
                        
                    }else{
                        if(npc!=0){
                            betapc=cbind(c(1:npc),myGLM$betapc)
                            colnames(betapc)=c("CV","Effect")
                            utils::write.csv(betapc,paste("FarmCPU.",name.of.trait2,".CVeffect.csv",sep=""),quote=F,row.names=FALSE)
                        }
                        
                        GAPIT.Report(name.of.trait=name.of.trait2,GWAS=GWAS,pred=NULL,ypred=ypred,tvalue=NULL,stderr=stderr,Vp=Vp,DPP=DPP,cutOff=cutOff,threshold.output=threshold.output,MAF=MAF,seqQTN=QTN.position,MAF.calculate=MAF.calculate,plot.style=plot.style)
                    }
                }
                #Evaluate Power vs FDR and type I error
                #print("Calling Power..")
                myPower=GAPIT.Power(WS=WS, alpha=alpha, maxOut=maxOut,seqQTN=QTN.position,GM=GM,GWAS=GWAS,MaxBP=1e10)
            } #enf of is done
            #if(length(seqQTN)==1) maxLoop=3
        } #end of while loop
        print("**********FarmCPU ACCOMPLISHED SUCCESSFULLY**********")
        #print(name.of.trait)
        #print("-----------------------------------------------------------------------")
        #===============================================================================
    }# end of loop on trait
    
    if(ncpus>1)snowfall::sfStop()
    gc()
    if(ncol(Y)==2) {
        # return (list(GWAS=GWAS,GPS=NULL,Pred=pred,compression=NULL,kinship.optimum=NULL,kinship=NULL,ycor=ycor,FDR=myPower$FDR,Power=myPower$Power,Power.Alpha=myPower$Power.Alpha,alpha=myPower$alpha,betapc=myGLM$betapc,seqQTN=seqQTN))
        return (list(GWAS=GWAS,GPS=NULL,Pred=pred,compression=NULL,kinship.optimum=NULL,kinship=NULL,ycor=ycor,betapc=myGLM$betapc,seqQTN=seqQTN))

    }else{
        return (list(GWAS=NULL,GPS=NULL,Pred=NULL,compression=NULL,kinship.optimum=NULL,kinship=NULL))
    }
    
}#The FarmCPU function ends here



`FarmCPU.Remove` <-
function(GDP=NULL,GM=NULL,seqQTN=NULL,seqQTN.p=NULL,orientation="col",threshold=.99){
    #Objective: Remove bins that are highly correlated
    #Input: GDP - n by m+1 matrix. The first colum is taxa name. The rest are m genotype
    #Input: GM - m by 3  matrix for SNP name, chromosome and BP
    #Input: seqQTN - s by 1 vecter for index of QTN on GM (+1 for GDP column wise)
    #Requirement: GDP and GM have the same order on SNP
    #Output: bin - n by s0 matrix of genotype
    #Output: binmap - s0 by 3 matrix for map of bin
    #Output: seqQTN - s0 by 1 vecter for index of QTN on GM (+1 for GDP column wise)
    #Relationship: bin=GDP[,c(seqQTN)], binmap=GM[seqQTN,], s0<=s
    #Authors: Zhiwu Zhang
    # Last update: March 4, 2013
    ##############################################################################
    #print("FarmCPU.Remove Started")
    #print(date())
    
    if(is.null(seqQTN))return(list(bin=NULL,binmap=NULL,seqQTN=NULL))
    #remove seqQTN with unsignificant p values
    #index.p=seqQTN.p<0.01
    #seqQTN.p=seqQTN.p[index.p]
    #seqQTN=seqQTN[index.p]
    #sort seqQTN using p values
    seqQTN=seqQTN[order(seqQTN.p)]
    
    hugeNum=10e10
    n=length(seqQTN)
    #print("Number of bins and GDP")
    #print(n)
    #print(dim(GDP))
    #print(seqQTN)
    
    #fielter bins by physical location
    
    binmap=GM[seqQTN,]
    
    #print("binmap")
    #print(binmap)
    
    cb=as.numeric(binmap[,2])*hugeNum+as.numeric(binmap[,3])#create ID for chromosome and bp
    cb.unique=unique(cb)
    
    #print("debuge")
    #print(cb)
    #print(cb.unique)
    
    index=match(cb.unique,cb,nomatch = 0)
    seqQTN=seqQTN[index]
    
    #print("Number of bins after chr and bp fillter")
    n=length(seqQTN) #update n
    #print(n)
    #print(date())
    
    #Set sample
    ratio=.1
    maxNum=100000
    if(orientation=="col"){
        s=nrow(GDP) #sample size
        m=ncol(GDP) #number of markers
    }else{
        m=nrow(GDP) #sample size
        s=ncol(GDP) #number of markers
    }
    
    #print("Determine number of samples")
    #print(date())
    #sampled=floor(ratio*s)
    sampled=s
    if(sampled>maxNum)sampled=maxNum
    
    #print("Number of individuals sampled to test dependency of bins")
    #print(sampled)
    
    #index=sample(s,sampled)
    index=1:sampled
    
    #print("Get the samples")
    #print(date())
    
    #This section has problem of turning big.matrix to R matrix
    #It is OK as x is small
    if(orientation=="col"){
        if(bigmemory::is.big.matrix(GDP)){
            x=as.matrix(bigmemory::deepcopy(GDP,rows=index,cols=seqQTN) )
        }else{
            x=GDP[index,seqQTN]
        }
    }else{
        if(bigmemory::is.big.matrix(GDP)){
            x=t(as.matrix(bigmemory::deepcopy(GDP,rows=seqQTN,cols=index) ))
        }else{
            x=t(GDP[seqQTN,index] )
        }
    }# end of if orientation
    
    #print("Calculating r")
    #print(date())
    #print("matrix x")
    #print(is(x))
    #print(dim(x))
    #print(length(x))
    
    #x=x[,order(seqQTN.p)]
    #print("x")
    #print(head(x))
    r = stats::cor(as.matrix(x))
    #print("r")
    #print(r)
    #print("indexing r")
    #print(date())
    index=abs(r)>threshold
    
    #print("index")
    #print(index)
    #print("Fancy algorithm")
    #print(date())
    #print("dimension of r")
    #print(dim(r))
    b=r*0
    b[index]=1
    c=1-b
    #print("for loop")
    #print(date())
    
    #for(i in 1:(n-1)){
    #  for (j in (i+1):n){
    #    b[j,j]=b[j,j]*c[i,j]
    #  }
    #}
    
    #The above are replaced by following
    c[lower.tri(c)]=1
    diag(c)=1
    bd <- apply(c,2,prod)
    
    #print("Positioning...")
    #print(date())
    
    #position=diag(b)==1
    position=(bd==1)
    seqQTN=seqQTN[position]
    #============================end of optimum============================================
    seqQTN=seqQTN[!is.na(seqQTN)]
    
    #print("Extract bin genotype data")
    #print(date())
    
    #This section has problem of turning big.matrix to R matrix
    
    if(orientation=="col"){
        if(bigmemory::is.big.matrix(GDP)){
            bin=as.matrix(bigmemory::deepcopy(GDP,cols=seqQTN) )
        }else{
            bin=GDP[,seqQTN]
        }
    }else{
        if(bigmemory::is.big.matrix(GDP)){
            bin=t(as.matrix(bigmemory::deepcopy(GDP,rows=seqQTN,) ))
        }else{
            bin=t(GDP[seqQTN,] )
        }
    }# end of if orientation
    
    
    #print("Get bin map")
    #print(date())
    
    binmap=GM[seqQTN,]
    
    #print("Number of bins left:")
    #print(length(seqQTN))
    #print("FarmCPU.Remove accomplished successfully!")
    
    return(list(bin=bin,binmap=binmap,seqQTN=seqQTN))
}#The function FarmCPU.Remove ends here
`FarmCPU.SUB` <-
function(GM=NULL,GLM=NULL,QTN=NULL,method="mean",useapply=TRUE,model="A"){
    #Input: FarmCPU.GLM object
    #Input: QTN - s by 3  matrix for SNP name, chromosome and BP
    #Input: method - options are "penalty", "reward","mean","median",and "onsite"
    #Requirement: P has row name of SNP. s<=t. covariates of QTNs are next to SNP
    #Output: GLM with the last column of P updated by the substituded p values
    #Authors: Xiaolei Liu and Zhiwu Zhang
    # Last update: Febuary 26, 2013
    ##############################################################################
    if(is.null(GLM$P)) return(NULL)  #P is required
    if(is.null(QTN)) return(NULL)  #QTN is required
    #print("FarmCPU.SUB Started")
    #print("dimension of QTN")
    #print(dim(QTN))
    #print(length(QTN))
    
    #print("debug")
    #print(QTN)
    #print(GLM)
    #position=match(QTN[,1], rownames(GLM$P), nomatch = 0)
    position=match(QTN[,1], GM[,1], nomatch = 0)
    #position=(1:nrow(GM))[GM[,1]%in%QTN[,1]]
    nqtn=length(position)
    #print("Position of QTN  on GM")
    #print(length(position))
    #print(position)
    #get position of QTNs (last nqtn columns from the second last)
    if(model=="A"){
        index=(ncol(GLM$P)-nqtn):(ncol(GLM$P)-1)
        spot=ncol(GLM$P)
    }else{
        index=(ncol(GLM$P)-nqtn-1):(ncol(GLM$P)-2)
        spot=ncol(GLM$P)-1
    }
    
    #print("Position of P value of QTN")
    #print(index)
    
    #print("Position of P value of marker")
    #print(spot)
    
    #print('ok')
    #print(ncol(GLM$P))
    #print(nqtn)
    #print((ncol(GLM$P)-nqtn))
    #print((ncol(GLM$P)-1))
    #print(min(GLM$P[,index],na.rm=TRUE))
    #print(GLM$P[position,spot])
    if(ncol(GLM$P)!=1){
        if(length(index)>1){
            if(method=="penalty") P.QTN=apply(GLM$P[,index],2,max,na.rm=TRUE)
            if(method=="reward") P.QTN=apply(GLM$P[,index],2,min,na.rm=TRUE)
            if(method=="mean") P.QTN=apply(GLM$P[,index],2,mean,na.rm=TRUE)
            if(method=="median") P.QTN=apply(GLM$P[, index], 2, stats::median, na.rm=TRUE)
            if(method=="onsite") P.QTN=GLM$P0[(length(GLM$P0)-nqtn+1):length(GLM$P0)]
        }else{
            if(method=="penalty") P.QTN=max(GLM$P[,index],na.rm=TRUE)
            if(method=="reward") P.QTN=min(GLM$P[,index],na.rm=TRUE)
            if(method=="mean") P.QTN=mean(GLM$P[,index],na.rm=TRUE)
            if(method=="median") P.QTN = stats::median(GLM$P[,index], stats::median, na.rm=TRUE)
            if(method=="onsite") P.QTN=GLM$P0[(length(GLM$P0)-nqtn+1):length(GLM$P0)]
        }
        
        #replace SNP pvalues with QTN pvalue
        #print("Substituting...")
        GLM$P[position,spot]=P.QTN
        #print(position)
        #print(GLM$betapred)
        GLM$B[position,]=GLM$betapred
    }
    #write.table(P,file="debuger.csv",sep=",")
    return(GLM)
}#The function FarmCPU.SUB ends here

`FarmCPU.P.Threshold` <-
function(GD=NULL,GM=NULL,Y=NULL,trait="",theRep=100){
    #Input: GD - Genotype
    #Input: GM - SNP name, chromosome and BP
    #Input: Y - phenotype, 2 columns
    #Input: trait - name of the trait
    #Input: theRep - number of replicates
    #Output: get minimum p value of each permutation and the recommend p.threshold used for FarmCPU model
    #Authors: Xiaolei Liu
    # Last update: April 6, 2015
    ##############################################################################
    
    #theRep=theRep
    #trait=trait
    if(is.null(GD))return(NULL)
    if(is.null(GM))return(NULL)
    if(is.null(Y))return(NULL)
    set.seed(12345)
    i=1
    for(i in 1:theRep){
        index=1:nrow(Y)
        index.shuffle=sample(index,length(index),replace=F)
        Y.shuffle=Y
        Y.shuffle[,2]=Y.shuffle[index.shuffle,2]
        
        #GWAS with FarmCPU...
        myFarmCPU=FarmCPU(
            Y=Y.shuffle[,c(1,2)],#Phenotype
            GD=GD,#Genotype
            GM=GM,#Map information
            file.output=FALSE,
            method.bin="optimum", #options are "static" and "optimum", default is static and this gives the fastest speed. If you want to use random model to optimize possible QTNs selection, use method.bin="optimum"
            maxLoop=1,#maxLoop is used to set the maximum iterations you want
            iteration.output=TRUE,#iteration.output=TRUE means to output results of every iteration
        )
        
        pvalue=min(myFarmCPU$GWAS[,4],na.rm=T)
        if(i==1){
            pvalue.final=pvalue
        }else{
            pvalue.final=c(pvalue.final,pvalue)
        }
    }#end of theRep
    
    utils::write.table(pvalue.final,paste("FarmCPU.p.threshold.optimize.",trait,".txt",sep=""),sep="\t",col.names=F,quote=F,row.names=F)
    
    print("The p.threshold of this data set should be:")
    print(sort(pvalue.final)[ceiling(theRep*0.05)])
    
}#end of `FarmCPU.P.Threshold`


`FarmCPU.Burger` <-
function(Y=NULL,CV=NULL,GK=NULL){
    #Object: To calculate likelihood, variances and ratio, revised by Xiaolei based on GAPIT.Burger function from GAPIT package
    #Straitegy: NA
    #Output: P value
    #intput:
    #Y: phenotype with columns of taxa,Y1,Y2...
    #CV: covariate variables with columns of taxa,v1,v2...
    #GK: Genotype data in numerical format, taxa goes to row and snp go to columns. the first column is taxa (same as GAPIT.bread)
    #Authors: Xiaolei Liu ,Jiabo Wang and Zhiwu Zhang
    #Last update: Dec 21, 2016
    ##############################################################################################
    #print("FarmCPU.Burger in progress...")
    
    if(!is.null(CV)){
        CV=as.matrix(CV)#change CV to a matrix when it is a vector xiaolei changed here
        theCV=as.matrix(cbind(matrix(1,nrow(CV),1),CV)) ###########for FarmCPU
    }else{
        theCV=matrix(1,nrow(Y),1)
    }
    
    #handler of single column GK
    n=nrow(GK)
    m=ncol(GK)
    if(m>2){
        theGK=as.matrix(GK)#GK is pure genotype matrix
    }else{
        theGK=matrix(GK,n,1)
    }
    
    myFaSTREML=GAPIT.get.LL(pheno=matrix(Y[,-1],nrow(Y),1),geno=NULL,snp.pool=theGK,X0=theCV)
    REMLs=-2*myFaSTREML$LL
    delta=myFaSTREML$delta
    vg=myFaSTREML$vg
    ve=myFaSTREML$ve
    
    #print("FarmCPU.Burger succeed!")
    return (list(REMLs=REMLs,vg=vg,ve=ve,delta=delta))
} #end of FarmCPU.Burger
#=============================================================================================



