`GAPIT.get.LL` <-
cmpfun(function(pheno,geno=NULL,snp.pool,X0=NULL){
    # evaluation of the maximum likelihood
    #Input: ys, xs, vg, delta, Z, X0, snp.pool
    #Output: LL
    #Authors: Qishan Wang, Feng Tian and Zhiwu Zhang
    #Last update: April 16, 2012
    ################################################################################
    #print("GAPIT.get.LL started")
    #print("dimension of pheno, snpool and X0")
    #print(dim(pheno))
    #print(length(pheno))
    #print(dim(snp.pool))
    #print(length(snp.pool))
    #print(dim(X0))
    #print(length(X0))
    
    y=pheno
    p=0
    deltaExpStart = -5
    deltaExpEnd = 5
    snp.pool=snp.pool[,]
    if(!is.null(snp.pool)&&var(snp.pool)==0){
        deltaExpStart = 100
        deltaExpEnd = deltaExpStart
        #print("deltaExp change here")
    }
    if(is.null(X0)) {
        X0 = matrix(1, nrow(snp.pool), 1)
    }
    #snp.test=as.numeric(geno[,1])
    #X <- cbind(X0, snp.test)
    X=X0
    
    #########SVD of X
    K.X.svd= svd(snp.pool,LINPACK=TRUE)######rivised by Jiabo Wang 2016.1.8
    # snp.pool=NA problem occurred
    #####rivised 2012.4.15 by qishan wang
    d=K.X.svd$d
    d=d[d>1e-08]
    d=d^2
    U1=K.X.svd$u
    U1=U1[,1:length(d)] ##rivised 2012.4.15 by qishan wang
    
    #handler of single snp
    if(is.null(dim(U1))) U1=matrix(U1,ncol=1)

    
    ###################
    n=nrow(U1)
    #I= diag(1,nrow(U1)) #xiaolei removed, this costs lots of memory
    
    U1TX=crossprod(U1,X)
    U1TY=crossprod(U1,y)
    yU1TY<- y-U1%*%U1TY
    XU1TX<- X-U1%*%U1TX  ### i is out of bracket
    #xiaolei rewrite following 4 lines
    IU = -tcrossprod(U1,U1)
    diag(IU) = rep(1,n) + diag(IU)
    #IUU=(I-tcrossprod(U1,U1))
    IUX=crossprod(IU,X )
    IUY=crossprod(IU,y)
    
    #Iteration on the range of delta (-5 to 5 in glog scale)
    for (m in seq(deltaExpStart,deltaExpEnd,by=0.1))
    {
        p=p+1
        delta<- exp(m)
        
        #----------------------------calculate beta-------------------------------------
        #######get beta compnents 1
        beta1=0
        for(i in 1:length(d)){
            one=matrix(U1TX[i,], nrow=1)
            beta=crossprod(one,(one/(d[i]+delta)))  #This is not real beta, confusing
            beta1= beta1+beta
        }
        
        #######get beta components 2
        beta2=0
        for(i in 1:nrow(U1)){
            one=matrix(IUX[i,], nrow=1)
            dim(one)
            beta=crossprod(one,one)
            beta2= beta2+beta
        }
        beta2<-beta2/delta
        
        #######get b3
        beta3=0
        for(i in 1:length(d)){
            one1=matrix(U1TX[i,], nrow=1)
            one2=matrix(U1TY[i,], nrow=1)
            beta=crossprod(one1,(one2/(d[i]+delta)))  #This is not real beta, confusing
            beta3= beta3+beta
        }
        
        ###########get beta4
        beta4=0
        for(i in 1:nrow(U1)){
            one1=matrix(IUX[i,], nrow=1)
            one2=matrix(IUY[i,], nrow=1)
            beta=crossprod(one1,one2)       #This is not real beta, confusing
            beta4= beta4+beta
        }
        beta4<-beta4/delta
        
        #######get final beta
        #zw1=solve(beta1+beta2)
        zw1 <- try(solve(beta1+beta2),silent=TRUE)
        if(inherits(zw1, "try-error")){
            zw1 <- ginv(beta1+beta2)
        }
        
        #zw1=ginv(beta1+beta2)
        zw2=(beta3+beta4)
        beta=crossprod(zw1,zw2)  #This is the real beta
        
        #----------------------------calculate LL---------------------------------------
        ####part 1
        part11<-n*log(2*3.14)
        part12<-0
        for(i in 1:length(d)){
            part12_pre=log(d[i]+delta)
            part12= part12+part12_pre
        }
        part13<- (nrow(U1)-length(d))*log(delta)
        part1<- -1/2*(part11+part12+part13)
        
        ######  part2
        part21<-nrow(U1)
        ######part221
        
        part221=0
        for(i in 1:length(d)){
            one1=matrix(U1TX[i,], nrow=1)
            one2=matrix(U1TY[i,], nrow=1)
            part221_pre=(one2-one1%*%beta)^2/(d[i]+delta) ###### beta contain covariate and snp %*%
            part221= part221+part221_pre
        }
        
        ######part222
        part222=0
        
        for(i in 1:n){
            one1=matrix(XU1TX[i,], nrow=1)
            one2=matrix(yU1TY[i,], nrow=1)
            part222_pre=((one2-one1%*%beta)^2)/delta
            part222= part222+part222_pre
        }
        part22<-n*log((1/n)*(part221+part222))
        part2<- -1/2*(part21+part22)
        
        ################# likihood
        LL<-part1+part2
        part1<-0
        part2<-0
        
        #-----------------------Save the optimum---------------------------------------
        if(p==1){
            beta.save=beta
            delta.save=delta
            LL.save=LL
        }else{
            if(LL>LL.save){
                beta.save=beta
                delta.save=delta
                LL.save=LL
            }
        }
        
    } # end of Iteration on the range of delta (-5 to 5 in glog scale)
    
    #--------------------update with the optimum------------------------------------
    beta=beta.save
    delta=delta.save
    LL=LL.save
    names(delta)=NULL
    names(LL)=NULL
    
    #--------------------calculating Va and Vem-------------------------------------
    #sigma_a1
    #U1TX=crossprod(U1,X)#xiaolei removed, it is re-calculated
    #U1TY=crossprod(U1,y)#xiaolei removed, it is re-calculated
    sigma_a1=0
    for(i in 1:length(d)){
        one1=matrix(U1TX[i,], nrow=1)
        one2=matrix(U1TY[i,], nrow=1)
        sigma_a1_pre=(one2-one1%*%beta)^2/(d[i]+delta)
        sigma_a1= sigma_a1+sigma_a1_pre
    }
    
    ### sigma_a2
    #xiaolei removed following 3 lines
    #IU=I-tcrossprod(U1,U1)    #This needs to be done only once
    #IUX=crossprod(IU,X)
    #IUY=crossprod(IU,y)
    sigma_a2=0
    
    for(i in 1:nrow(U1)){
        one1=matrix(IUX[i,], nrow=1)
        one2=matrix(IUY[i,], nrow=1)
        sigma_a2_pre<-(one2-one1%*%beta)^2
        sigma_a2= sigma_a2+sigma_a2_pre
    }
    
    sigma_a2<-sigma_a2/delta
    sigma_a<- 1/n*(sigma_a1+sigma_a2)
    sigma_e<-delta*sigma_a
    
    return(list(beta=beta, delta=delta, LL=LL, vg=sigma_a,ve=sigma_e))
}
)#end of cmpfun(
#=============================================================================================

