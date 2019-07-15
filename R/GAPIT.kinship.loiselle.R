`GAPIT.kinship.loiselle` <-
function(snps, method="additive", use="all") {
# Object: To calculate the kinship matrix using the method of Loiselle et al. (1995)
# Authors: Alex Lipka and Hyun Min Kang
# Last update: May 31, 2011 
############################################################################################## 
  #Number of SNP types that are 0s
  n0 <- sum(snps==0,na.rm=TRUE)
  #Number of heterozygote SNP types
  nh <- sum(snps==0.5,na.rm=TRUE)
  #Number of SNP types that are 1s
  n1 <- sum(snps==1,na.rm=TRUE)
  #Number of SNP types that are missing
  nNA <- sum(is.na(snps))
  

 
  #Self explanatory
  dim(snps)[1]*dim(snps)[2]
  #stopifnot(n0+nh+n1+nNA == length(snps))

    
  #Note that the two lines in if(method == "dominant") and if(method == "recessive") are found in
  #if(method == "additive").  Worry about this only if you have heterozygotes, which you do not.
  if( method == "dominant" ) {
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
    snps[!is.na(snps) && (snps == 0.5)] <- flags[!is.na(snps) && (snps == 0.5)]
  }
  else if( method == "recessive" ) {
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
    snps[!is.na(snps) && (snps == 0.5)] <- flags[!is.na(snps) && (snps == 0.5)]
  }
  else if( ( method == "additive" ) && ( nh > 0 ) ) {
    dsnps <- snps
    rsnps <- snps
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
    dsnps[!is.na(snps) && (snps==0.5)] <- flags[is.na(snps) && (snps==0.5)]
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
    rsnps[!is.na(snps) && (snps==0.5)] <- flags[is.na(snps) && (snps==0.5)]
    snps <- rbind(dsnps,rsnps)
  }

  #mafs is a (# SNPs)x(# lines) matrix.  The columns of mafs are identical, and the ij^th element is the average
  #allele frequency for the SNP in the i^th row.
  
  #if(use == "all") imputes missing SNP type values with the expected (average) allele frequency.
  if( use == "all" ) {
    mafs <- matrix(rowMeans(snps,na.rm=TRUE),nrow(snps),ncol(snps))
    snps[is.na(snps)] <- mafs[is.na(snps)]
  }
  else if( use == "complete.obs" ) {
    mafs <- matrix(rowMeans(snps,na.rm=TRUE),nrow(snps),ncol(snps))
    snps <- snps[rowSums(is.na(snps))==0,]
  }
  mafs_comp <- 1-mafs
  snps_comp <- 1-snps
  

  n <- ncol(snps)
  K <- matrix(nrow=n,ncol=n)
  diag(K) <- 1
  #Create the k term on page 1422 of Loiselle et al. (1995)

  missing <- rep(NA, dim(snps)[1])  
  for(i in 1:dim(snps)[1]) {
    missing[i] <- sum(is.na(snps[i,]))
  }
  

  for(i in 1:(n-1)) {
    for(j in (i+1):n) {
      Num_First_Term_1 <- (snps[,i]-mafs[,i])*(snps[,j]-mafs[,j])
      Num_First_Term_2 <- (snps_comp[,i]-mafs_comp[,i])*(snps_comp[,j]-mafs_comp[,j])
      First_Term <- sum(Num_First_Term_1)+sum(Num_First_Term_2)

      Num_Second_Term_1 <- mafs[,i]*(1-mafs[,i])
      Num_Second_Term_2 <- mafs_comp[,i]*(1-mafs_comp[,i])
      Num_Second_Term_Bias_Correction <- 1/((2*n)-missing - 1)
      Num_Second_Term <-  Num_Second_Term_1 + Num_Second_Term_2
      Second_Term <- sum(Num_Second_Term*Num_Second_Term_Bias_Correction)

      Third_Term <- sum(Num_Second_Term) 
      
      f <- (First_Term + Second_Term)/Third_Term

      K[i,j] <- f
      if(K[i,j]<0) K[i,j]=0
      
      K[j,i] <- K[i,j]
    }
  }
  return(K)
}
#=============================================================================================

