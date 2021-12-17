


#' Blink
#'
#' @description 
#' Blink
#' 
#' @param Y = NULL, data.frame of phenotypic data, column one is sample (taxa), column two is a trait, subsequent columns are other traits.
#' @param GD = NULL, data.frame of genetic data in 'numerical' format, samples in rows, variants in columns.
#' @param GM = NULL, Genetic Map data.frame to provide genomic coordinates for GD
#' @param QTN.position = NULL,
#' @param CV = NULL, Covariates
#' @param DPP = 100000000,
#' @param kinship.algorithm = "FARM-CPU",
#' @param file.output = TRUE,
#' @param cutOff = 0.01,
#' @param method.GLM = "FarmCPU.LM",
#' @param Prior = NULL,
#' @param ncpus = 1,
#' @param maxLoop = 10,
#' @param LD = 0.7,
#' @param threshold.output = .0001,
#' @param alpha = c(.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),
#' @param WS = c(1e0,1e3,1e4,1e5,1e6,1e7),
#' @param GP = NULL,
#' @param FDRcut = FALSE,
#' @param maxOut = 10,
#' @param converge = 1,
#' @param iteration.output = FALSE,
#' @param acceleration = 0,
#' @param threshold = NA,
#' @param model = "A",
#' @param MAF.calculate = FALSE,
#' @param plot.style = "FarmCPU",
#' @param p.threshold = NA,
#' @param maf.threshold = 0,
#' @param bound = FALSE,
#' @param method.sub = "reward",
#' @param method.sub.final = "reward",
#' @param stepwise = FALSE,
#' @param BIC.method = "naive",
#' @param LD.wise = FALSE,
#' @param time.cal = FALSE,
#' @param Prediction  =  F
#' 
#' 
#' @return 
#' list(GWAS=GWAS,myGLM=myGLM,PEV = PEV,seqQTN=seqQTN)
#'
#' @export
`Blink` <- function(Y = NULL,
                    QTN.position = NULL,
                    GD = NULL,
                    GM = NULL,
                    CV = NULL,
                    DPP = 100000000,
                    kinship.algorithm = "FARM-CPU",
                    file.output = TRUE,
                    cutOff = 0.01,
                    method.GLM = "FarmCPU.LM",
                    Prior = NULL,
                    ncpus = 1,
                    maxLoop = 10,
                    LD = 0.7,
                    threshold.output = .0001,
                    alpha = c(.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),
                    WS = c(1e0,1e3,1e4,1e5,1e6,1e7),
                    GP = NULL,
                    FDRcut = FALSE,
                    maxOut = 10,
                    converge = 1,
                    iteration.output = FALSE,
                    acceleration = 0,
                    threshold = NA,
                    model = "A",
                    MAF.calculate = FALSE,
                    plot.style = "FarmCPU",
                    p.threshold = NA,
                    maf.threshold = 0,
                    bound = FALSE,
                    method.sub = "reward",
                    method.sub.final = "reward",
                    stepwise = FALSE,
                    BIC.method = "naive",
                    LD.wise = FALSE,
                    time.cal = FALSE,
                    Prediction  =  FALSE){
  # Jiabo modified the Blink GS codes in 2020.9
  print("----------------------Welcome to Blink----------------------")
  time.start=proc.time()
  nm=nrow(GM)
  ny=nrow(Y)
  if(is.na(threshold)){
    threshold = floor(ny / log(ny))
  }
  ngd = nrow(GD)
  orientation = "col"
  theSNP = 2
  ns = nrow(GD)
  seqQTN = NULL
  if(nm == ngd){
    orientation = "row"
    theSNP = 1
    ns = ncol(GD)
  }
  if(maf.threshold > 0) {
    MAF.calculate = TRUE
  }

  if(MAF.calculate==FALSE){
    MAF=NA
  }else{
    MAF=apply(GD,theSNP,mean)
    MAF=matrix(MAF,nrow=1)
    MAF=apply(X = MAF, 
              MARGIN = 2,
              function(x){ min(1 - x/2, x/2) }
              )
  }
  MAF.index = 1:nm
  if(maf.threshold > 0) {
    MAF.index = MAF > maf.threshold
    MAF = MAF[MAF.index]
  }
  ac=NULL
  if(acceleration != 0){
    ac = rep(1.0, nm)
  }
  index = which(GM[,3] == 0 )
  if(length(index) > 0){
    GM[index,3]=1
  }

  P = GP
  gc()
  if(ncol(GD) > nm & orientation == "col"){
    if( bigmemory::is.big.matrix(GD) ){
      GD = bigmemory::deepcopy(GD, 
                               rows=1:nrow(GD), 
                               cols=2:ncol(GD))
    }else{
      GD=as.matrix(GD[,-1])
    }
  }
  # GD=as.matrix(GD)
  gc()
  shift = 0
  for(trait in 2:2){
    name.of.trait = colnames(Y)[trait]
    index = MAF.index
    seqTaxa = which(!is.na(Y[,trait]))
    Y1 = Y[seqTaxa,]
    if(!is.null(CV)){
        CV1 = CV[seqTaxa,] #Thanks for jloat's suggestion in Jul 23 2021
    }else{
      CV1=NULL
    }
    if(orientation == "col"){
      if(bigmemory::is.big.matrix(GD)){
        GD1=bigmemory::deepcopy(GD,rows=seqTaxa,cols=index)
      }else{
        GD1=GD[seqTaxa,index]
      }
    } else {
      if(bigmemory::is.big.matrix(GD)){
        GD1=bigmemory::deepcopy(GD,rows=index,cols=seqTaxa)
      }else{
        GD1=GD[index,seqTaxa]
        GD1=as.matrix(GD1)
      }
    }
    LD.time = rep(0,maxLoop)
    BIC.time = rep(0,maxLoop)
    GLM.time = rep(0,maxLoop)
    theLoop = 0
    theConverge = 0
    seqQTN.save = c(0)
    isDone = FALSE
    name.of.trait2 = name.of.trait
    while(!isDone) {
      theLoop = theLoop + 1
      print(paste("----------------------Iteration:",theLoop,"----------------------",sep=" "))
      if(iteration.output){
        name.of.trait2 = paste("Iteration_",
                               theLoop,".",
                               name.of.trait,
                               sep="")
      }
      myPrior = FarmCPU.Prior(GM = GM,
                              P = P,
                              Prior = Prior,
                              kinship.algorithm = kinship.algorithm)
      if(!is.null(myPrior)){
        if(theLoop!=1){
          seqQTN.p = myPrior
          if(theLoop == 2){
            bonferroniCutOff = cutOff/nm
            sp = sort(seqQTN.p)
            spd = abs(cutOff - sp * nm/cutOff)
            index_fdr = grep(min(spd), spd)[1]
            FDRcutoff = cutOff * index_fdr/nm
            if(FDRcut){
              index.p = seqQTN.p < (FDRcutoff)
            }else{
              index.p = seqQTN.p < (bonferroniCutOff)
            }
            if(!is.na(p.threshold)){
              index.p = seqQTN.p < p.threshold
            }
            index.p[ is.na(index.p) ] = FALSE
            seqQTN.selected = as.numeric( which( index.p ) )
            }else{
              index.p = seqQTN.p < (1/nm)
              if(!is.na(p.threshold)){
                index.p=seqQTN.p<p.threshold
              }
              index.p[is.na(index.p)] = FALSE
              seqQTN.selected = as.numeric(which(index.p))
            }

            Porder = order(myPrior[seqQTN.selected], 
                           na.last = TRUE, 
                           decreasing = FALSE)
            t1 = proc.time()
            if(length(Porder) > 1){
              
              if(LD.wise & (length(Porder) > threshold)){
                max_Porder = max(Porder)
                if(max_Porder > 10000) max_Porder = 10000
                step_bin = 10
                Porder_new = rep(Porder[1],threshold)
                Po = 1
                for( po in 2:max_Porder){
                  if(min(abs(seqQTN.selected[Porder[po]]-seqQTN.selected[Porder_new]))>step_bin){
                    Po = Po + 1
                    Porder_new[Po]=Porder[po]
                    if (Po >=threshold) break
                  }
                }
                Porder=Porder_new
              }
              
              if(bigmemory::is.big.matrix(GD1)){
                if(orientation=="col"){
                  GDnew=bigmemory::deepcopy(GD1,cols=seqQTN.selected)
                  GDneo=bigmemory::deepcopy(GDnew,cols=Porder)
                }else{
                  GDnew=bigmemory::deepcopy(GD1,rows=seqQTN.selected)
                  GDneo=bigmemory::deepcopy(GDnew,rows=Porder)
                }
              } else {
                if(orientation=="col"){
                  GDnew=GD1[,seqQTN.selected]
                  GDneo=GDnew[,Porder]
                }else{
                  GDnew=GD1[seqQTN.selected,]
                  GDneo=GDnew[Porder,]
                }
              }
              
              print("LD remove is working....")
              print("Number SNPs for LD remove:")
              print(length(Porder))
              Psort=Blink.LDRemove(Porder=Porder,GDneo=GDneo,bound=bound,LD=LD,model=model,orientation=orientation)
              seqQTN.can=seqQTN.selected[Psort]
              t2=proc.time()
              print("Model selection based on BIC is working....")
              print("Number of SNPs for BIC selection:")
              print(length(seqQTN.can))
              myBIC = Blink.BICselection(Psort = seqQTN.can,
                                         GD = GD1,
                                         Y = Y1,
                                         orientation = orientation,
                                         BIC.method = BIC.method)
              seqQTN = myBIC$seqQTN
              #if(theLoop==6) print(seqQTN)
              t3 = proc.time()
              LD.time[theLoop] = as.numeric(t2)[3] - as.numeric(t1)[3]
              BIC.time[theLoop] = as.numeric(t3)[3] - as.numeric(t2)[3]
            }else if(length(Porder) == 1){
              print("LD remove is working....")
              print("Model selection based on BIC is working....")
              seqQTN=seqQTN.selected
            }else{
              seqQTN=NULL
            }
          }
        }else{
          seqQTN=NULL
        }

        if(theLoop==2){
          if(!is.na(p.threshold)){
            if(min(myPrior,na.rm=TRUE)>p.threshold){
              seqQTN=NULL
                print("Top snps have little effect, set seqQTN to NULL!")
            }
          }else{
            sp=sort(seqQTN.p)
            spd=abs(cutOff-sp*nm/cutOff)
            index_fdr=grep(min(spd),spd)[1]
            FDRcutoff=cutOff*index_fdr/nm

            if(FDRcut){
              index.cutoff=FDRcutoff
            }else{
              index.cutoff=bonferroniCutOff
            }
            # index.p=seqQTN.p<(FDRcutoff)
            if(min(myPrior,na.rm=TRUE)>index.cutoff){
              seqQTN=NULL
              print("Top snps have little effect, set seqQTN to NULL!")
            }
          }
        }
      
      if(theLoop==2&&is.null(seqQTN)|length(seqQTN)==0&&theLoop==2){
          print(paste("seqQTN is:",seqQTN,",stop here",sep=""))
          if(!isDone | iteration.output){
          gc()
              p.GLM=GWAS[,4]
                p.GLM.log=-log10(stats::quantile(p.GLM,na.rm=TRUE,0.05))
                bonf.log=1.3
                bonf.compare=p.GLM.log/bonf.log
                p.FARMCPU.log=-log10(p.GLM)/bonf.compare
              GWAS[,4]=10^(-p.FARMCPU.log)
                GWAS[,4][which(GWAS[,4]>1)]=1
                colnames(GWAS)=c(colnames(GM),"P.value","maf","nobs","Rsquare.of.Model.without.SNP","Rsquare.of.Model.with.SNP","FDR_Adjusted_P-values")
                Vp=stats::var(Y1[,2],na.rm=TRUE)
                # if(file.output) GAPIT.Report(name.of.trait=name.of.trait2,GWAS=GWAS,pred=NULL,tvalue=NULL,stderr=stderr,Vp=Vp,DPP=DPP,cutOff=cutOff,threshold.output=threshold.output,MAF=MAF,seqQTN=QTN.position,MAF.calculate=MAF.calculate,plot.style=plot.style)
                    # myPower=GAPIT.Power(WS=WS, alpha=alpha, maxOut=maxOut,seqQTN=QTN.position,GM=GM,GWAS=GWAS,MaxBP=1e10)
            }
              break
      }
      if(theLoop>1){
        if(all(seqQTN.save!=0 & seqQTN.save!=-1 & !is.null(seqQTN))){
          seqQTN=union(seqQTN,seqQTN.save)
        }
      }
      if(theLoop>2 ){
        if( length(Porder)>1){
          BIC=Blink.BICselection(Psort=seqQTN,GD=GD1,Y=Y1,orientation=orientation,BIC.method=BIC.method)
          seqQTN = BIC$seqQTN
        }
      }

      print("seqQTN:")
      print(seqQTN)
      theConverge=length(intersect(seqQTN,seqQTN.save))/length(union(seqQTN,seqQTN.save))
      isDone=((theLoop>=maxLoop)|(theConverge>=converge))
      if(!is.null(seqQTN)) seqQTN.save=seqQTN
      gc()
      if(!is.null(seqQTN)){
        if(orientation=="col"){
          theCV=cbind(CV1,GD1[,seqQTN])
        }else{
          if(length(seqQTN)>1){
            theCV1=t(GD1[seqQTN,])
            theCV=cbind(CV1,theCV1)
          }else{
            theCV=cbind(CV1,GD1[seqQTN,])
          }
        }
      }else{
        theCV=CV1
      }
      t4=proc.time()
      if(!is.null(theCV)) theCV=as.matrix(theCV)
      #if(theLoop==4) write.table(theCV,"CV.txt",col.names=F,row.names=F,quote=F,sep="\t")
      myGLM=FarmCPU.LM(y=Y1[,trait],GDP=GD1,w=theCV,orientation=orientation)
        #print(dim(myGLM$P))
        #myGLM=Blink.LM(y=Y1[,trait],GDP=GD1,w=theCV,orientation=orientation)
        #print(dim(myGLM$P))

        if(!isDone){
          myGLM=Blink.SUB(GM=GM,GLM=myGLM,QTN=GM[seqQTN,],method=method.sub,model=model)
        }else{
          #save(myGLM,file="myGLM_last.Rdata")
            myGLM=Blink.SUB(GM=GM,GLM=myGLM,QTN=GM[seqQTN,],method=method.sub,model=model)
        }
      t5=proc.time()
      GLM.time[theLoop]=as.numeric(t5)[3]-as.numeric(t4)[3]
      P=myGLM$P[,ncol(myGLM$P)-shift]
      index=which(ac>1)
      P[P==0] <- min(P[P!=0],na.rm=TRUE)*0.01
      P[is.na(P)] =1
      gc()
      nf=ncol(myGLM$P)/4
      tvalue=myGLM$P[,nf*2-shift]
      stderr=myGLM$P[,3*nf-shift]
      GWAS=cbind(GM[MAF.index,],P,MAF,NA,NA,NA,NA)
      colnames(GWAS)=c(colnames(GM),"P.value","maf","nobs","Rsquare.of.Model.without.SNP","Rsquare.of.Model.with.SNP","FDR_Adjusted_P-values")
      Vp=stats::var(Y1[,2],na.rm=TRUE)
      # if(file.output){
      #   if(theLoop==1&&is.null(CV)){
      #     GAPIT.Report(name.of.trait=name.of.trait2,GWAS=GWAS,pred=NULL,ypred=NULL,tvalue=NULL,stderr=stderr,Vp=Vp,DPP=DPP,cutOff=cutOff,threshold.output=threshold.output,MAF=MAF,seqQTN=QTN.position,MAF.calculate=MAF.calculate,plot.style=plot.style)
      #   }else{
      #     GAPIT.Report(name.of.trait=name.of.trait2,GWAS=GWAS,pred=NULL,ypred=NULL,tvalue=NULL,stderr=stderr,Vp=Vp,DPP=DPP,cutOff=cutOff,threshold.output=threshold.output,MAF=MAF,seqQTN=QTN.position,MAF.calculate=MAF.calculate,plot.style=plot.style)
      #   }
      # }# end of file.out
    } #end of theLoop
    PEV=NULL
    if(Prediction){
      YP = cbind(Y[,1],Y[,trait])
      p.rank = order(GWAS[,4],na.last = T,decreasing=F)
      if(sum(GWAS[p.rank,4]<p.threshold)>0){
        seqQTN = p.rank[GWAS[p.rank,4] < p.threshold]
      }else{
        seqQTN = p.rank[1]
      }
      if(!bigmemory::is.big.matrix(GD)){
        if(orientation=="col"){
          GDpred = GD[,seqQTN]
        }else{
          if(length(seqQTN)>1){
            GDpred = t(GD[seqQTN,])
          }else{
            GDpred = as.matrix(GD[seqQTN,])
          }
        }
      }else{
        if(orientation=="col"){
          GDpred = bigmemory::deepcopy(GD,cols=seqQTN)
        }else{
          GDpred = bigmemory::deepcopy(GD,rows=seqQTN)
        }
      }
      PEV = Blink.Pred(Y = YP,
                       GD = GDpred,
                       CV = CV,
                       orientation = orientation
                       )
    }
    if(time.cal){
      print("LD.time(sec):")
      print(LD.time[1:theLoop])
      print("BIC.time(sec):")
      print(BIC.time[1:theLoop])
      print("GLM.time(sec):")
      print(GLM.time[1:theLoop])
    }
    time.end=proc.time()
    time.all=as.numeric(time.end)[3]-as.numeric(time.start)[3]
    print(paste("-------------Blink finished successfully in",round(time.all,2),"seconds!-----------------"))
  # print(proc.time())
    # write.table(GWAS,paste(name.of.trait2,"_GWAS.txt",sep=""),sep="\t",col.names=T,row.names=F)
  }#end of phenotype
  return(list(GWAS=GWAS,myGLM=myGLM,PEV = PEV,seqQTN=seqQTN))
}#  end of function Blink


# GWAS4_blinkc=result
# blinkc=merge(GWAS4[,c(1,4)],GWAS4_blinkc[,c(1,5)],by.x="SNP",by.y="taxa")




# cor(blinkc[,2],blinkc[,3])
# cor(farmcpu[,2],farmcpu[,3])


# `Blink.BICselection` <-  function(Psort = NULL,
#                                   CV = NULL,
#                                   GD = NULL,
#                                   Y = Y1,
#                                   orientation = NULL,
#                                   BIC.method = "even"){

`Blink.BICselection` <-  function(Y, 
                                  Psort = NULL,
                                  CV = NULL,
                                  GD = NULL,
                                  orientation = NULL,
                                  BIC.method = "even"){
#Objects: fixed model selection using BIC
#Input:Y,GD,Psort
#   BIC.method: Naive: detect all SNPs of Psort
#         even: detect some SNPs, step=floor(sqrt(m))+1
#         fixed: detect the SNPs by fixed steps. Default is 20
#         log: detect the SNPs by log(10,N) transform
#         ln: detect the SNPs by ln transform
#Output: seqQTN: SNP position
#Author: Yao Zhou
#Last update: 01/05/2016, modified 03/31/2016
  GD = as.matrix(GD)
  n=nrow(Y)
  threshold=floor(n/log(n))
  if(threshold < length(Psort)){
    seqQTN=Psort[1:threshold]
  }else{
    seqQTN=Psort
  }
  y=Y[,2]

  s=0
  a=0
  m=length(seqQTN)
  pmatrix=matrix(1,m,m)

  if(BIC.method=="naive"){
    position=seq(1:m)

  } else if(BIC.method=="even"){
    step.length=floor(sqrt(m))+1
    step=floor(m/step.length)
    if ((m-step*step.length)>=(0.5*step.length)) {
          step=step+1
      }
    if (step.length>m) {
      step.length=m
            step=1
      }
    position=seq(step,m,step)
    if(position[length(position)]<m) position=append(position,m)
  } else if(BIC.method=="lg"){
    if(m==1){
      position =c(1)
    }else{
      le=seq(1:m)
      step=le/log10(le)
      for(i in 2:m){
        le[i]=le[i-1]+step[i]
        le=round(le)
        if(le[i]>m){
          position=le[1:i]
          break
        }
      }
    }

  } else if(BIC.method=="ln"){
    if(m==1){
      position =c(1)
    }else{
      le=seq(1:m)
      step=le/log(le)

      for(i in 2:m){
        le[i]=le[i-1]+step[i]
        le=round(le)
        if(le[i]>m){
          position=le[1:i]
          break
        }
      }
    }
  } else if(BIC.method=="fixed"){
    if(m>20){
      position=floor(seq(1,m,m/20))
    }else{
      position=seq(1:m)
    }
  }else{
    print("please choose one method for BIC")
#    break
  }

  BICv=rep(NA,length(position))
  if(is.null(CV)){
    w=as.matrix(rep(1,n))
    ww=n
    ncov=2
  }else{
    CV=as.matrix(CV)
    w=cbind(1,CV)
    ww=crossprod(w)
    ncov=ncol(ww)+1
  }
  wwi=MASS::ginv(ww)


  pos.pre=0
  k=0
  for(pos in position){
    if(pos>m) pos=m
    if(orientation=="col"){
      x=GD[,seqQTN[(pos.pre+1):pos]]
    }else{
      x=GD[seqQTN[(pos.pre+1):pos],]
      if(is.matrix(x)){
        x=t(x)
      }else{
        x=as.matrix(x)
      }
    }
    k=k+1
    pos.pre=pos
    x=as.matrix(x)

    if(k==1){
      ww=crossprod(w,w)
    }else{
      WW=matrix(0,(nwc+nxc),(nwc+nxc))
      WW[1:nwc,1:nwc]=ww
      WW[1:nwc,(nwc+1):(nwc+nxc)]=xw
      WW[(nwc+1):(nxc+nwc),1:nwc]=wx
      WW[(nwc+1):(nwc+nxc),(nwc+1):(nwc+nxc)]=xx
      ww=WW
    }
    nwc = ncol(w)
    nxc = ncol(x)
        iXX = matrix(0,(nwc+nxc),(nwc+nxc))
    xx = crossprod(x,x)
    xw = crossprod(x,w)
    wx = crossprod(w,x)
    t1 = wwi %*% wx
    t2 = xx - xw %*% t1
    if (!is.null(t2)){
    M22 = MASS::ginv(t2)
    t3=xw %*% wwi
    M21=-M22 %*% t3
    M12=-t1 %*% M22
    M11=wwi + t1 %*% M22 %*% t3
    iXX[1:nwc,1:nwc]=M11
    iXX[(nwc+1):(nwc+nxc),(nwc+1):(nwc+nxc)]=M22
    iXX[(nwc+1):(nwc+nxc),1:nwc]=M21
    iXX[1:nwc,(nwc+1):(nwc+nxc)]=M12
    w=cbind(w,x)
    wy=crossprod(w,y)
    wwi=iXX
    beta=wwi %*% wy
    yp= w %*% beta
    ve=as.numeric(stats::var(yp-y))
    RSS= (yp-y)^2
    n2LL=n*log(2*pi)+n*log(ve)+2*sum(RSS/(2*ve))
    # BICv[k]=n2LL+2*(nwc+nxc-1)*log(n)
    BICv[k]=n2LL+(nwc+nxc-1)*log(n)
    df=(n-pos-1)
    MSE=sum(RSS)/df
    se=sqrt(diag(iXX)*MSE)
    tvalue=beta/se
        pvalue <- 2 * stats::pt(abs(tvalue), df,lower.tail = FALSE)
        pmatrix[1:pos,pos]=pvalue[ncov:length(pvalue)]
    }
  }
  seqQTN=Psort[1:position[which(BICv==min(BICv,na.rm=T))]]
  pvalue=as.numeric(pmatrix[1:length(seqQTN),length(seqQTN)])
  return(list(seqQTN=seqQTN,pvalue=pvalue,BIC=BICv))
}


`Blink.LDRemoveBlock`<-function(GDneo=NULL,LD=NULL,Porder=NULL,bound=FALSE,model="A",orientation=NULL){
#`Blink.LDRemove`<-function(GDneo=NULL,LD=NULL,Porder=NULL,bound=FALSE,model="A",orientation=NULL){
#Objects: Calculate LD and remove the correlated SNPs
#Authors: Yao Zhou
#Last Update:  03/03/16
  if (model=="D"){
    GDneo=1-abs(GDneo-1)
  }

  GDneo=as.matrix(GDneo)
  if(min(ncol(GDneo),nrow(GDneo))<201) bound=FALSE
  if(orientation=="col"){
    n=nrow(GDneo)
    if(bound){
      GDneo=GDneo[sample(n,200,replace=F),]
    }
  }else{
    n=ncol(GDneo)
    if(bound){
      GDneo=GDneo[,sample(n,200,replace=F)]
    }
    GDneo=t(GDneo)
  }
  # cat("ncol(GDneo) is",ncol(GDneo),"\n")
  corr = stats::cor(GDneo)
  corr[is.na(corr)]=1
  corr[abs(corr)<=LD]=0
  corr[abs(corr)>LD]=1
  Psort=as.numeric(matrix(1,1,ncol(corr)))
  # print(ncol(corr))
  for(i in 2:ncol(corr)){
    p.a=Psort[1:(i-1)]
    p.b=as.numeric(corr[1:(i-1),i])
    index=(p.a==p.b)
    index[(p.a==0)&(p.b==0)]=FALSE
    if(sum(index)!=0) Psort[i]=0
  }
  seqQTN=Porder[Psort==1]
  return(seqQTN)
}

`Blink.LDRemove`<-function(GDneo=NULL,LD=0.7,Porder=NULL,bound=FALSE,model="A",orientation="row",block=1000,LD.num =50){
#Objects: LD remove, especially length(Porder)>10000
#Authors: Yao Zhou
#Last update: 08/15/2016
  GDneo = as.matrix(GDneo)
  if (orientation == "row") {
    SNP.index = apply(GDneo,1, stats::sd)!=0
    GDneo = GDneo[SNP.index,]
  } else {
    SNP.index = apply(GDneo, 2, stats::sd) != 0
    GDneo = GDneo[, SNP.index]
  }
  Porder = Porder[SNP.index]
  l = block
  seqQTN=NULL
  lp=length(Porder)
  k=ceiling(lp/l)
  GDneo=as.matrix(GDneo)
  if(min(ncol(GDneo),nrow(GDneo))<201) bound=FALSE
  if(orientation=="col"){
    n=nrow(GDneo)
    if(bound){
      GDneo=GDneo[sample(n,200,replace=F),]
    }
  }else{
    n=ncol(GDneo)
    if(bound){
      GDneo=GDneo[,sample(n,200,replace=F)]
    }
    GDneo=t(GDneo)
  }
  for(i in 1:k){
    bottom=(i-1)*l+1
    up=l*i
    if(up>lp) up = lp
    Porderb=Porder[bottom:up]

    index = seq(bottom:up)
    GDneob = GDneo[,index]
    # cat("i is ",i,"\n")
    # print(length(index))
    seqQTNs = Blink.LDRemoveBlock(GDneo=GDneob,LD=LD,Porder=Porderb,orientation="col",model=model)
    # print(seqQTN)
    seqQTN = append(seqQTN,seqQTNs)
    if(k >1){
      index1 = which(Porder %in% seqQTN)
      Porderb = Porder[index1]
      GDneob = GDneo[,index1]
      if(length(index1)>1){
        seqQTN = Blink.LDRemoveBlock(GDneo=GDneob,LD=LD,Porder=Porderb,orientation="col",model=model)
      }else{
        seqQTN = Porderb
      }

    }
    if(LD.num < length(seqQTN)) break
  }
  rm(GDneob,Porderb)
  return(seqQTN)
}

`Blink.LM` <-function(y,w=NULL,GDP,orientation="col"){
    N=length(y) #Total number of taxa, including missing ones
    direction=2
    if(orientation=="row"){
    GDP=t(GDP)
    }
    ntest=ncol(GDP)
    if(orientation=="row"){
        B <- matrix(NA,nrow=nrow(GDP),ncol=ncol(w)+1)
    }else{
        B <- matrix(NA,nrow=ncol(GDP),ncol=ncol(w)+1)
    }
    print(dim(B))
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
    
    n=N
    k=1 #number of genetic effect: 1 and 2 for A and AD respectively
    
    q1=(q0+1) # vecter index for the posistion of genetic effect (a)
    q2=(q0+1):(q0+2) # vecter index for the posistion of genetic effect (a and d)
    df=n-q0-k #residual df (this should be varied based on validating d)
    
    iXX=matrix(0,q0+k,q0+k) #Reserve the maximum size of inverse of LHS
    
    ww=crossprod(w,w)
    wy=crossprod(w,y)
    yy=crossprod(y,y)
    # wwi=solve(ww) Revised by Jiabo on 2021.3.4
    wwi <- try(solve(ww),silent=TRUE)
     if(inherits(wwi, "try-error")){
      print("!!!!!")
     wwi <- MASS::ginv(ww)
     }
    #Statistics on the reduced model without marker
    rhs=wy
    gc()
    y=as.matrix(y)
  gy=crossprod(GDP,y)
  gw=crossprod(w,GDP)
  bw=crossprod(gw,wwi)
  lbw=ncol(bw)
  P <- matrix(NA,nrow=ncol(GDP),ncol=4*(nf+1))
  for(i in 1:ntest){ 
    x=GDP[,i]
    xy=gy[i,1]
    xw=gw[,i]
        xx=crossprod(x,x)
     #   B21 <- crossprod(xw, wwi)
        B21=matrix(bw[i,],1,lbw)
      t2=B21%*%xw #I have problem of using crossprod and tcrossprod here
    B22 <- xx - t2
    invB22=1/B22
    NeginvB22B21 <- crossprod(-invB22,B21)
    iXX11 <- wwi + as.numeric(invB22)*crossprod(B21,B21)
        iXX[1:q0,1:q0]=iXX11
      iXX[q1,q1]=invB22
        iXX[q1,1:q0]=NeginvB22B21
        iXX[1:q0,q1]=NeginvB22B21
        rhs=c(wy,xy)
        beta <- crossprod(iXX,rhs)   #both a and d go in
        df=n-q0-1
        ve=(yy-crossprod(beta,rhs))/df
        se=sqrt(diag(iXX)*ve)
        tvalue=beta/se
        pvalue <- 2 * stats::pt(abs(tvalue), df,lower.tail = FALSE)
        if(!is.na(abs(B22[1,1]))){
            if(abs(B22[1,1])<10e-8)pvalue[]=NA}
        P[i,]=c(beta[-1],tvalue[-1],se[-1],pvalue[-1])
        print(length(beta))
    #     B[i,]=beta[length(beta)]
    }
    return(list(P=P,PF=NULL,beta=beta))
} #end of function

`Blink.Pred` <- function(GD = NULL, 
                         Y = NULL,
                         CV = NULL, 
                         orientation = "col"){
## Objects: Prediction using significant pseudo QTNs
## Input: Y, CV and GD
## Output: Predicted Phenotype
## Authors: Yao Zhou
## Last update: 2/6/2017

  if(bigmemory::is.big.matrix(GD)) GD = as.matrix(GD)
  if(orientation =="row"){
    GD = t(GD)
    if(nrow(GD)==1) GD = t(GD)
  }
  
  seqTaxa=which(!is.na(Y[,2]))
  Y1 = Y[seqTaxa,2]
  GD1 = GD[seqTaxa,]
  
  if(is.null(CV)){
    mylm = stats::lm(Y1 ~ GD1)
    PEV = stats::predict(mylm,as.data.frame(GD))
  }else{
    CV1 = CV[seqTaxa,]
    mylm = stats::lm(Y1 ~ CV1 + GD1)
    PEV = stats::predict(mylm,as.data.frame(cbind(CV,GD)))
  }
  return(PEV) 
}


`Blink.SUB` <-
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
    position=match(QTN[,1], GM[,1], nomatch = 0)
    nqtn=length(position)
    if(model=="A"){
        index=(ncol(GLM$P)-nqtn):(ncol(GLM$P)-1)
        spot=ncol(GLM$P)
    }else{
        index=(ncol(GLM$P)-nqtn-1):(ncol(GLM$P)-2)
        spot=ncol(GLM$P)-1
    }
    if(ncol(GLM$P)!=1){
        if(length(index)>1){
            if(method=="penalty") P.QTN=apply(GLM$P[,index],2,max,na.rm=TRUE)
            if(method=="reward") P.QTN=apply(GLM$P[,index],2,min,na.rm=TRUE)
            if(method=="mean") P.QTN=apply(GLM$P[,index],2,mean,na.rm=TRUE)
            if(method=="median") P.QTN = apply(GLM$P[,index], 2, stats::median, na.rm = TRUE)
            if(method=="onsite") P.QTN=GLM$P0[(length(GLM$P0)-nqtn+1):length(GLM$P0)]
        }else{
            if(method=="penalty") P.QTN=max(GLM$P[,index],na.rm=TRUE)
            if(method=="reward") P.QTN=min(GLM$P[,index],na.rm=TRUE)
            if(method=="mean") P.QTN=mean(GLM$P[,index],na.rm=TRUE)
            if(method=="median") P.QTN=stats::median(GLM$P[,index], stats::median,na.rm=TRUE)
            if(method=="onsite") P.QTN=GLM$P0[(length(GLM$P0)-nqtn+1):length(GLM$P0)]
        }
        GLM$P[position,spot]=P.QTN
    }
    return(GLM)
}#The function FarmCPU.SUB ends here


`Blink.cor`<-function(Y,
                      GD,
                      w = NULL,
                      orientation = "row",
                      ms = ms,
#                      n = ny,
#                      m = nm
                      n = nrow(Y),
                      m = nrow(GD)
                      ){
  #Objects: calculate R value with covariates
  #Input: pheontype(nx1), ms is marker size for slicing the genotype, genotype(orientation="row", mxn or orientation="col", nxm,) and covariates(nxp)
  #   n is individual number, m is marker number, p is covariate number
  #Output: abs(r)
  #Author: Yao Zhou
  #Last updated: Jun 28, 2016
  if(!is.matrix(Y)) Y=as.matrix(Y)
  # Orthogonolize phenotype w.r.t. covariates
  {
    if(!is.null(w)){
      w = cbind(1,w)
    }else{
      w = matrix(1,n,1)
    }
    if(!is.matrix(w)) w = as.matrix(w)
    qw = qr(w)
    if( min(abs(diag(qr.R(qw)))) < .Machine$double.eps * m ) {
      stop("Colinear or zero covariates detected");
    }
    w = qr.Q(qw)
    tw=t(w)
    rm(qw)
  } 
  
  # Orthogonolize phenotype w.r.t. covariates
  
  {
    Y = Y - w%*%crossprod(w,Y)
    colsq = colSums(Y^2)
    div = sqrt(colsq)
    Y = Y/div
    rm(colsq,div)
  }
  time.start = proc.time()
  #Orthogonolize genotype w.r.t. covariates
  {
    if(orientation == "row"){
      rabs = matrix(NA,nrow = nrow(GD),ncol = 1)
      ntest = nrow(GD)
      ns = ceiling(ntest/ms)
      for(i in 1:ns){
        bottom=(ms*(i-1)+1)
        if(i<ns){
          up=ms*i
        }else{
          up=ntest
        }
        GDs = GD[bottom:up,]
        GDs = GDs - tcrossprod(GDs,tw)%*%tw
        rowsq= rowSums(GDs^2)
        div = sqrt(rowsq)
        GDs=GDs/div
        rabs[bottom:up,1]=abs(GDs%*%Y)
      }
    }else{
      rabs = matrix(NA,nrow = ncol(GD),ncol = 1)
      ntest = ncol(GD)
      ns = ceiling(ntest/ms)
      for(i in 1:ns){
        bottom=(ms*(i-1)+1)
        if(i<ns){
          up=ms*i
        }else{
          up=ntest
        }
        GDs = GD[,bottom:up]
        GDs = GDs- crossprod(tw,tw%*%GDs)
        colsq= colSums(GDs^2)
        div = sqrt(colsq)
        GDs=t(GDs)/div
        rabs[bottom:up,1] = abs(GDs%*%Y)
      }
    }
    rm(GDs,div)     
  }    
  time.end= proc.time()
  time.all = time.end - time.start
  print(time.all)
  return(rabs)
}


`Blink.ptor`<-function(p,df){
#Objects: transform the p value to r value
#Input: p: p value
#   df: degree of freedom, df = n-ncov-1
#   
#Output: r value
#Author: Yao Zhou
#Last Update: 8/15/2015
  t = stats::qt(0.5*p,df-1,lower.tail = FALSE)
  r=sqrt(t^2/(df-1+t^2))
  return(r)
}

`Blink.rtop`<-function(r,df){
#Objects: transform the p value to r value
#Input: r: r value
#   df: degree of freedom, df = n-ncov-1
#   
#Output: p value
#Author: Yao Zhou
#Last Update: 8/15/2015
  tvalue=sqrt(df-1)*r/sqrt(1-r^2)
  pvalue <- 2 * stats::pt(abs(tvalue), df-1,lower.tail = FALSE)
  return(pvalue)
}

# `BlinkR.SUB` <-function(CV,seqQTN,Y,r,ny,ms,m){
# #Objects: subsitution of r value for covariates
# #Input: Y, nx1 vector, phenotype
# #   w, all covariates without 1 
# #   seq = length(seqQTN), number of SNPs added as covariate
# #Output：r, r value for SNPs added as covariates
# #Author:Yao Zhou
# #Last update: 08/15/2016
#   rsnp=matrix(NA,nsnp,1)
#   ncov=ncol(CV)
#   nf=ncov-nsnp
#   GDP=CV[,(nf+1):ncov]
#   for(i in 1:nsnp){
#     w=GDP[,-i]
#     if(nsnp==1) w=NULL
#     GD=as.matrix(GDP[,i])
#     rsnp[i,1]=Blink.cor(Y=Y,w=w,GD=GD,orientation=orientation,ms=ms,n=ny,m=nm)
#   }
#   rm(GDP, GD, w, ncov, nf)
#     return(rsnp)
# }

