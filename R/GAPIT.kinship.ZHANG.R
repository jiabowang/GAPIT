`GAPIT.kinship.ZHANG` <-
  function(snps,hasInbred=TRUE) {
    # Object: To calculate ZHANG (Zones Harbored Adjustments of Negligent Genetic) relationship
    # Authors: Zhwiu Zhang
    # Last update: october 25, 2014 
    ############################################################################################## 
    print("Calculating ZHANG relationship defined by Zhiwu Zhang...")
    #Remove invariants
    fa=colSums(snps)/(2*nrow(snps))
    index.non=fa>=1| fa<=0
    snps=snps[,!index.non]
    
    het=1-abs(snps-1)
    ind.sum=rowSums(het)
    fi=ind.sum/(2*ncol(snps))
    inbreeding=1-min(fi)
    
    nSNP=ncol(snps)
    nInd=nrow(snps)
    n=nInd 
    snpMean= apply(snps,2,mean)   #get mean for each snp
    print("substracting mean...")
    snps=t(snps)-snpMean    #operation on matrix and vector goes in direction of column
    print("Getting X'X...")
    #K=tcrossprod((snps), (snps))
    K=crossprod((snps), (snps)) 
    if(is.na(K[1,1])) stop ("GAPIT says: Missing data is not allowed for numerical genotype data")
    
    print("Adjusting...")
    #Extract diagonals
    i =1:n
    j=(i-1)*n
    index=i+j
    d=K[index]
    DL=min(d)
    DU=max(d)
    floor=min(K)
    
    
    #Set range between 0 and 2
    top=1+inbreeding
    K=top*(K-floor)/(DU-floor)
    Dmin=top*(DL-floor)/(DU-floor)
    
    #Adjust based on expected minimum diagonal (1)
    if(Dmin<1) {
      print("Adjustment by the minimum diagonal")
      K[index]=(K[index]-Dmin+1)/((top+1-Dmin)*.5)
      K[-index]=K[-index]*(1/Dmin)
    }
    
    #Limiting the maximum offdiagonal to the top
    Omax=max(K[-index])
    if(Omax>top){
      print("Adjustment by the minimum off diagonal")
      K[-index]=K[-index]*(top/Omax)
    }
    
    print("Calculating kinship with Zhang method: done")
    return(K)
  }
#=============================================================================================

