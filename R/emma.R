emma.kinship <- function(snps, method="additive", use="all") {
  n0 <- sum(snps==0,na.rm=TRUE)
  nh <- sum(snps==0.5,na.rm=TRUE)                                                                                         
  n1 <- sum(snps==1,na.rm=TRUE)
  nNA <- sum(is.na(snps))

  stopifnot(n0+nh+n1+nNA == length(snps))

  if ( method == "dominant" ) {
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
    snps[!is.na(snps) && (snps == 0.5)] <- flags[!is.na(snps) && (snps == 0.5)]
  }
  else if ( method == "recessive" ) {
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
    snps[!is.na(snps) && (snps == 0.5)] <- flags[!is.na(snps) && (snps == 0.5)]
  }
  else if ( ( method == "additive" ) && ( nh > 0 ) ) {
    dsnps <- snps
    rsnps <- snps
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
    dsnps[!is.na(snps) && (snps==0.5)] <- flags[is.na(snps) && (snps==0.5)]
    flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
    rsnps[!is.na(snps) && (snps==0.5)] <- flags[is.na(snps) && (snps==0.5)]
    snps <- rbind(dsnps,rsnps)
  }

  if ( use == "all" ) {
    mafs <- matrix(rowMeans(snps,na.rm=TRUE),nrow(snps),ncol(snps))
    snps[is.na(snps)] <- mafs[is.na(snps)]
  }
  else if ( use == "complete.obs" ) {
    snps <- snps[rowSums(is.na(snps))==0,]
  }

  n <- ncol(snps)
  K <- matrix(nrow=n,ncol=n)
  diag(K) <- 1

  for(i in 1:(n-1)) {
    for(j in (i+1):n) {
      x <- snps[,i]*snps[,j] + (1-snps[,i])*(1-snps[,j])
      K[i,j] <- sum(x,na.rm=TRUE)/sum(!is.na(x))
      K[j,i] <- K[i,j]
    }
  }
  return(K)
}

emma.eigen.L <- function(Z,K,complete=TRUE) {
  if ( is.null(Z) ) {
    return(emma.eigen.L.wo.Z(K))
  }
  else {
    return(emma.eigen.L.w.Z(Z,K,complete))
  }
}

emma.eigen.L.wo.Z <- function(K) {
  eig <- eigen(K,symmetric=TRUE)
  return(list(values=eig$values,vectors=eig$vectors))
}

emma.eigen.L.w.Z <- function(Z,K,complete=TRUE) {
  if ( complete == FALSE ) {
    vids <- colSums(Z)>0
    Z <- Z[,vids]
    K <- K[vids,vids]
  }
  eig <- eigen(K%*%crossprod(Z,Z),symmetric=FALSE,EISPACK=TRUE)
  return(list(values=eig$values,vectors=qr.Q(qr(Z%*%eig$vectors),complete=TRUE)))
}

emma.eigen.R <- function(Z,K,X,complete=TRUE) {
  if ( is.null(Z) ) {
    return(emma.eigen.R.wo.Z(K,X))
  }
  else {
    return(emma.eigen.R.w.Z(Z,K,X,complete))
  }
}

emma.eigen.R.wo.Z <- function(K, X) {
  n <- nrow(X)
  q <- ncol(X)
  S <- diag(n)-X%*%solve(crossprod(X,X))%*%t(X)
  eig <- eigen(S%*%(K+diag(1,n))%*%S,symmetric=TRUE)
  stopifnot(!is.complex(eig$values))
  return(list(values=eig$values[1:(n-q)]-1,vectors=eig$vectors[,1:(n-q)]))
}

emma.eigen.R.w.Z <- function(Z, K, X, complete = TRUE) {
  if ( complete == FALSE ) {
    vids <-  colSums(Z) > 0
    Z <- Z[,vids]
    K <- K[vids,vids]
  }
  n <- nrow(Z)
  t <- ncol(Z)
  q <- ncol(X)
 

  
  SZ <- Z - X%*%solve(crossprod(X,X))%*%crossprod(X,Z)
  eig <- eigen(K%*%crossprod(Z,SZ),symmetric=FALSE,EISPACK=TRUE)
  if ( is.complex(eig$values) ) {
    eig$values <- Re(eig$values)
    eig$vectors <- Re(eig$vectors)    
  }
  qr.X <- qr.Q(qr(X))
  return(list(values=eig$values[1:(t-q)],
              vectors=qr.Q(qr(cbind(SZ%*%eig$vectors[,1:(t-q)],qr.X)),
                complete=TRUE)[,c(1:(t-q),(t+1):n)]))   
}

emma.delta.ML.LL.wo.Z <- function(logdelta, lambda, etas, xi) {
  n <- length(xi)
  delta <- exp(logdelta)
  return( 0.5*(n*(log(n/(2*pi))-1-log(sum((etas*etas)/(lambda+delta))))-sum(log(xi+delta))) )  
}

emma.delta.ML.LL.w.Z <- function(logdelta, lambda, etas.1, xi.1, n, etas.2.sq ) {
  t <- length(xi.1)
  delta <- exp(logdelta)
#  stopifnot(length(lambda) == length(etas.1))
  return( 0.5*(n*(log(n/(2*pi))-1-log(sum(etas.1*etas.1/(lambda+delta))+etas.2.sq/delta))-(sum(log(xi.1+delta))+(n-t)*logdelta)) )
}

emma.delta.ML.dLL.wo.Z <- function(logdelta, lambda, etas, xi) {
  n <- length(xi)
  delta <- exp(logdelta)
  etasq <- etas*etas
  ldelta <- lambda+delta
  return( 0.5*(n*sum(etasq/(ldelta*ldelta))/sum(etasq/ldelta)-sum(1/(xi+delta))) )
}

emma.delta.ML.dLL.w.Z <- function(logdelta, lambda, etas.1, xi.1, n, etas.2.sq ) {
  t <- length(xi.1)
  q <- t-length(lambda)
  delta <- exp(logdelta)
  etasq <- etas.1*etas.1
  ldelta <- lambda+delta
  return( 0.5*(n*(sum(etasq/(ldelta*ldelta))+etas.2.sq/(delta*delta))/(sum(etasq/ldelta)+etas.2.sq/delta)-(sum(1/(xi.1+delta))+(n-t)/delta) ) )
}

emma.delta.REML.LL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <-  exp(logdelta)
  return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas*etas/(lambda+delta))))-sum(log(lambda+delta))) )
}

emma.delta.REML.LL.w.Z <- function(logdelta, lambda, etas.1, n, t, etas.2.sq ) {
  tq <- length(etas.1)
  nq <- n - t + tq
  delta <-  exp(logdelta)
  return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas.1*etas.1/(lambda+delta))+etas.2.sq/delta))-(sum(log(lambda+delta))+(n-t)*logdelta)) ) 
}

emma.delta.REML.dLL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <- exp(logdelta)
  etasq <- etas*etas
  ldelta <- lambda+delta
  return( 0.5*(nq*sum(etasq/(ldelta*ldelta))/sum(etasq/ldelta)-sum(1/ldelta)) )
}

emma.delta.REML.dLL.w.Z <- function(logdelta, lambda, etas.1, n, t1, etas.2.sq ) {
  t <- t1
  tq <- length(etas.1)
  nq <- n - t + tq
  delta <- exp(logdelta)
  etasq <- etas.1*etas.1
  ldelta <- lambda+delta
  return( 0.5*(nq*(sum(etasq/(ldelta*ldelta))+etas.2.sq/(delta*delta))/(sum(etasq/ldelta)+etas.2.sq/delta)-(sum(1/ldelta)+(n-t)/delta)) )
}

emma.MLE <- function(y, X, K, Z=NULL, ngrids=100, llim=-10, ulim=10,
  esp=1e-10, eig.L = NULL, eig.R = NULL)
{
  n <- length(y)
  t <- nrow(K)
  q <- ncol(X)
  
#  stopifnot(nrow(K) == t)
  stopifnot(ncol(K) == t)
  stopifnot(nrow(X) == n)

  if ( det(crossprod(X,X)) == 0 ) {
    warning("X is singular")
    return (list(ML=0,delta=0,ve=0,vg=0))
  }

  if ( is.null(Z) ) {
    if ( is.null(eig.L) ) {
      eig.L <- emma.eigen.L.wo.Z(K)
    }
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.wo.Z(K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    
  
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas <- matrix(eig.R$values,n-q,m) + matrix(delta,n-q,m,byrow=TRUE)
    Xis <- matrix(eig.L$values,n,m) + matrix(delta,n,m,byrow=TRUE)
    Etasq <- matrix(etas*etas,n-q,m)
    LL <- 0.5*(n*(log(n/(2*pi))-1-log(colSums(Etasq/Lambdas)))-colSums(log(Xis)))
    dLL <- 0.5*delta*(n*colSums(Etasq/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(1/Xis))
    
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.ML.LL.wo.Z(llim,eig.R$values,etas,eig.L$values))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.ML.LL.wo.Z(ulim,eig.R$values,etas,eig.L$values))
    }

    for( i in 1:(m-1) )
      {
        if ( ( dLL[i]*dLL[i+1] < 0 ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
        {
          r <- uniroot(emma.delta.ML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas, xi=eig.L$values)
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.ML.LL.wo.Z(r$root,eig.R$values, etas, eig.L$values))
        }
      }
#    optdelta <- exp(optlogdelta)
  }
  else {
    if ( is.null(eig.L) ) {
      eig.L <- emma.eigen.L.w.Z(Z,K)
    }
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.w.Z(Z,K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    etas.1 <- etas[1:(t-q)]
    etas.2 <- etas[(t-q+1):(n-q)]
    etas.2.sq <- sum(etas.2*etas.2)

    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim

    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas <- matrix(eig.R$values,t-q,m) + matrix(delta,t-q,m,byrow=TRUE)
    Xis <- matrix(eig.L$values,t,m) + matrix(delta,t,m,byrow=TRUE)
    Etasq <- matrix(etas.1*etas.1,t-q,m)
    #LL <- 0.5*(n*(log(n/(2*pi))-1-log(colSums(Etasq/Lambdas)+etas.2.sq/delta))-colSums(log(Xis))+(n-t)*log(deltas))
    dLL <- 0.5*delta*(n*(colSums(Etasq/(Lambdas*Lambdas))+etas.2.sq/(delta*delta))/(colSums(Etasq/Lambdas)+etas.2.sq/delta)-(colSums(1/Xis)+(n-t)/delta))
    
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.ML.LL.w.Z(llim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.ML.LL.w.Z(ulim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
    }

    for( i in 1:(m-1) )
      {
        if ( ( dLL[i]*dLL[i+1] < 0 ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
        {
          r <- uniroot(emma.delta.ML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, xi.1=eig.L$values, n=n, etas.2.sq = etas.2.sq )
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.ML.LL.w.Z(r$root,eig.R$values, etas.1, eig.L$values, n, etas.2.sq ))
        }
      }
#    optdelta <- exp(optlogdelta)
  }

  maxdelta <- exp(optlogdelta[which.max(optLL)])
  maxLL <- max(optLL)
  if ( is.null(Z) ) {
    maxva <- sum(etas*etas/(eig.R$values+maxdelta))/n    
  }
  else {
    maxva <- (sum(etas.1*etas.1/(eig.R$values+maxdelta))+etas.2.sq/maxdelta)/n
  }
  maxve <- maxva*maxdelta

  return (list(ML=maxLL,delta=maxdelta,ve=maxve,vg=maxva))
}

emma.REMLE <- function(y, X, K, Z=NULL, ngrids=100, llim=-10, ulim=10,
  esp=1e-10, eig.L = NULL, eig.R = NULL) {
  n <- length(y)
  t <- nrow(K)
  q <- ncol(X)

#  stopifnot(nrow(K) == t)
  stopifnot(ncol(K) == t)
  stopifnot(nrow(X) == n)

  if ( det(crossprod(X,X)) == 0 ) {
    warning("X is singular")
    return (list(REML=0,delta=0,ve=0,vg=0))
  }

  if ( is.null(Z) ) {
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.wo.Z(K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
  
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas <- matrix(eig.R$values,n-q,m) + matrix(delta,n-q,m,byrow=TRUE)
    Etasq <- matrix(etas*etas,n-q,m)
    LL <- 0.5*((n-q)*(log((n-q)/(2*pi))-1-log(colSums(Etasq/Lambdas)))-colSums(log(Lambdas)))
    dLL <- 0.5*delta*((n-q)*colSums(Etasq/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(1/Lambdas))
    
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.REML.LL.wo.Z(llim,eig.R$values,etas))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.REML.LL.wo.Z(ulim,eig.R$values,etas))
    }

    for( i in 1:(m-1) )
      {
        if ( ( dLL[i]*dLL[i+1] < 0 ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
        {
          r <- uniroot(emma.delta.REML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas)
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.REML.LL.wo.Z(r$root,eig.R$values, etas))
        }
      }
#    optdelta <- exp(optlogdelta)
  }
  else {
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.w.Z(Z,K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    etas.1 <- etas[1:(t-q)]
    etas.2 <- etas[(t-q+1):(n-q)]
    etas.2.sq <- sum(etas.2*etas.2)
  
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas <- matrix(eig.R$values,t-q,m) + matrix(delta,t-q,m,byrow=TRUE)
    Etasq <- matrix(etas.1*etas.1,t-q,m)
    dLL <- 0.5*delta*((n-q)*(colSums(Etasq/(Lambdas*Lambdas))+etas.2.sq/(delta*delta))/(colSums(Etasq/Lambdas)+etas.2.sq/delta)-(colSums(1/Lambdas)+(n-t)/delta))
    
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.REML.LL.w.Z(llim,eig.R$values,etas.1,n,t,etas.2.sq))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.REML.LL.w.Z(ulim,eig.R$values,etas.1,n,t,etas.2.sq))
    }

    for( i in 1:(m-1) )
      {
        if ( ( dLL[i]*dLL[i+1] < 0 ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
        {
          r <- uniroot(emma.delta.REML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, n=n, t1=t, etas.2.sq = etas.2.sq )
          optlogdelta <- append(optlogdelta, r$root)
          optLL <- append(optLL, emma.delta.REML.LL.w.Z(r$root,eig.R$values, etas.1, n, t, etas.2.sq ))
        }
      }
#    optdelta <- exp(optlogdelta)
  }  

  maxdelta <- exp(optlogdelta[which.max(optLL)])
  maxLL <- max(optLL)
  if ( is.null(Z) ) {
    maxva <- sum(etas*etas/(eig.R$values+maxdelta))/(n-q)    
  }
  else {
    maxva <- (sum(etas.1*etas.1/(eig.R$values+maxdelta))+etas.2.sq/maxdelta)/(n-q)
  }
  maxve <- maxva*maxdelta

  return (list(REML=maxLL,delta=maxdelta,ve=maxve,vg=maxva))
}

emma.ML.LRT <- function(ys, xs, K, Z=NULL, X0 = NULL, ngrids=100, llim=-10, ulim=10, esp=1e-10, ponly = FALSE) {
  if ( is.null(dim(ys)) || ncol(ys) == 1 ) {
    ys <- matrix(ys,1,length(ys))
  }
  if ( is.null(dim(xs)) || ncol(xs) == 1 ) {
    xs <- matrix(xs,1,length(xs))
  }
  if ( is.null(X0) ) {
    X0 <- matrix(1,ncol(ys),1)
  }  
  
  g <- nrow(ys)
  n <- ncol(ys)
  m <- nrow(xs)
  t <- ncol(xs)
  q0 <- ncol(X0)
  q1 <- q0 + 1

  if ( !ponly ) {
    ML1s <- matrix(nrow=m,ncol=g)
    ML0s <- matrix(nrow=m,ncol=g)
    vgs <- matrix(nrow=m,ncol=g)
    ves <- matrix(nrow=m,ncol=g)
  }
  stats <- matrix(nrow=m,ncol=g)
  ps <- matrix(nrow=m,ncol=g)
  ML0 <- vector(length=g)
  
  stopifnot(nrow(K) == t)
  stopifnot(ncol(K) == t)
  stopifnot(nrow(X0) == n)

  if ( sum(is.na(ys)) == 0 ) {
    eig.L <- emma.eigen.L(Z,K)
    eig.R0 <- emma.eigen.R(Z,K,X0)
      
    for(i in 1:g) {
      ML0[i] <- emma.MLE(ys[i,],X0,K,Z,ngrids,llim,ulim,esp,eig.L,eig.R0)$ML
    }

    x.prev <- vector(length=0)
    
    for(i in 1:m) {
      vids <- !is.na(xs[i,])
      nv <- sum(vids)
      xv <- xs[i,vids]

      if ( ( mean(xv) <= 0 ) || ( mean(xv) >= 1 ) ) {
        if (!ponly) {
          stats[i,] <- rep(NA,g)
          vgs[i,] <- rep(NA,g)
          ves[i,] <- rep(NA,g)
          ML1s[i,] <- rep(NA,g)
          ML0s[i,] <- rep(NA,g)
        }
        ps[i,] = rep(1,g)
      }
      else if ( identical(x.prev, xv) ) {
        if ( !ponly ) {
          stats[i,] <- stats[i-1,]
          vgs[i,] <- vgs[i-1,]
          ves[i,] <- ves[i-1,]
          ML1s[i,] <- ML1s[i-1,]
          ML0s[i,] <- ML0s[i-1,]
        }
        ps[i,] <- ps[i-1,]
      }
      else {
        if ( is.null(Z) ) {
          X <- cbind(X0[vids,,drop=FALSE],xs[i,vids])
          eig.R1 = emma.eigen.R.wo.Z(K[vids,vids],X)
        }
        else {
          vrows <- as.logical(rowSums(Z[,vids]))
          nr <- sum(vrows)
          X <- cbind(X0[vrows,,drop=FALSE],Z[vrows,vids]%*%t(xs[i,vids,drop=FALSE]))
          eig.R1 = emma.eigen.R.w.Z(Z[vrows,vids],K[vids,vids],X)          
        }

        for(j in 1:g) {
          if ( nv == t ) {
            MLE <- emma.MLE(ys[j,],X,K,Z,ngrids,llim,ulim,esp,eig.L,eig.R1)
#            MLE <- emma.MLE(ys[j,],X,K,Z,ngrids,llim,ulim,esp,eig.L,eig.R1)            
            if (!ponly) { 
              ML1s[i,j] <- MLE$ML
              vgs[i,j] <- MLE$vg
              ves[i,j] <- MLE$ve
            }
            stats[i,j] <- 2*(MLE$ML-ML0[j])
            
          }
          else {
            if ( is.null(Z) ) {
              eig.L0 <- emma.eigen.L.wo.Z(K[vids,vids])
              MLE0 <- emma.MLE(ys[j,vids],X0[vids,,drop=FALSE],K[vids,vids],NULL,ngrids,llim,ulim,esp,eig.L0)
              MLE1 <- emma.MLE(ys[j,vids],X,K[vids,vids],NULL,ngrids,llim,ulim,esp,eig.L0)
            }
            else {
              if ( nr == n ) {
                MLE1 <- emma.MLE(ys[j,],X,K,Z,ngrids,llim,ulim,esp,eig.L)
              }
              else {
                eig.L0 <- emma.eigen.L.w.Z(Z[vrows,vids],K[vids,vids])              
                MLE0 <- emma.MLE(ys[j,vrows],X0[vrows,,drop=FALSE],K[vids,vids],Z[vrows,vids],ngrids,llim,ulim,esp,eig.L0)
                MLE1 <- emma.MLE(ys[j,vrows],X,K[vids,vids],Z[vrows,vids],ngrids,llim,ulim,esp,eig.L0)
              }
            }
            if (!ponly) { 
              ML1s[i,j] <- MLE1$ML
              ML0s[i,j] <- MLE0$ML
              vgs[i,j] <- MLE1$vg
              ves[i,j] <- MLE1$ve
            }
            stats[i,j] <- 2*(MLE1$ML-MLE0$ML)
          }
        }
        if ( ( nv == t ) && ( !ponly ) ) {
          ML0s[i,] <- ML0
        }
        ps[i,] <- pchisq(stats[i,],1,lower.tail=FALSE)
      }
    }
  }
  else {
    eig.L <- emma.eigen.L(Z,K)
    eig.R0 <- emma.eigen.R(Z,K,X0)
      
    for(i in 1:g) {
      vrows <- !is.na(ys[i,])      
      if ( is.null(Z) ) {
        ML0[i] <- emma.MLE(ys[i,vrows],X0[vrows,,drop=FALSE],K[vrows,vrows],NULL,ngrids,llim,ulim,esp)$ML
      }
      else {
        vids <- colSums(Z[vrows,]>0)
            
        ML0[i] <- emma.MLE(ys[i,vrows],X0[vrows,,drop=FALSE],K[vids,vids],Z[vrows,vids],ngrids,llim,ulim,esp)$ML        
      }
    }

    x.prev <- vector(length=0)
    
    for(i in 1:m) {
      vids <- !is.na(xs[i,])
      nv <- sum(vids)
      xv <- xs[i,vids]

      if ( ( mean(xv) <= 0 ) || ( mean(xv) >= 1 ) ) {
        if (!ponly) {
          stats[i,] <- rep(NA,g)
          vgs[i,] <- rep(NA,g)
          ves[i,] <- rep(NA,g)
          ML1s[i,] <- rep(NA,g)
          ML0s[,i] <- rep(NA,g)
        }
        ps[i,] = rep(1,g)
      }      
      else if ( identical(x.prev, xv) ) {
        if ( !ponly ) {
          stats[i,] <- stats[i-1,]
          vgs[i,] <- vgs[i-1,]
          ves[i,] <- ves[i-1,]
          ML1s[i,] <- ML1s[i-1,]
        }
        ps[i,] = ps[i-1,]
      }
      else {
        if ( is.null(Z) ) {
          X <- cbind(X0,xs[i,])
          if ( nv == t ) {
            eig.R1 = emma.eigen.R.wo.Z(K,X)
          }          
        }
        else {
          vrows <- as.logical(rowSums(Z[,vids]))
          X <- cbind(X0,Z[,vids,drop=FALSE]%*%t(xs[i,vids,drop=FALSE]))
          if ( nv == t ) {
            eig.R1 = emma.eigen.R.w.Z(Z,K,X)
          }
        }

        for(j in 1:g) {
#          print(j)
          vrows <- !is.na(ys[j,])
          if ( nv == t ) {
            nr <- sum(vrows)
            if ( is.null(Z) ) {
              if ( nr == n ) {
                MLE <- emma.MLE(ys[j,],X,K,NULL,ngrids,llim,ulim,esp,eig.L,eig.R1)                
              }
              else {
                MLE <- emma.MLE(ys[j,vrows],X[vrows,],K[vrows,vrows],NULL,ngrids,llim,ulim,esp)
              }
            }
            else {
              if ( nr == n ) {
                MLE <- emma.MLE(ys[j,],X,K,Z,ngrids,llim,ulim,esp,eig.L,eig.R1)                
              }
              else {
                vtids <- as.logical(colSums(Z[vrows,,drop=FALSE]))
                MLE <- emma.MLE(ys[j,vrows],X[vrows,],K[vtids,vtids],Z[vrows,vtids],ngrids,llim,ulim,esp)
              }
            }
            
            if (!ponly) { 
              ML1s[i,j] <- MLE$ML
              vgs[i,j] <- MLE$vg
              ves[i,j] <- MLE$ve
            }
            stats[i,j] <- 2*(MLE$ML-ML0[j])
          }
          else {
            if ( is.null(Z) ) {
              vtids <- vrows & vids
              eig.L0 <- emma.eigen.L(NULL,K[vtids,vtids])
              MLE0 <- emma.MLE(ys[j,vtids],X0[vtids,,drop=FALSE],K[vtids,vtids],NULL,ngrids,llim,ulim,esp,eig.L0)
              MLE1 <- emma.MLE(ys[j,vtids],X[vtids,],K[vtids,vtids],NULL,ngrids,llim,ulim,esp,eig.L0)
            }
            else {
              vtids <- as.logical(colSums(Z[vrows,])) & vids
              vtrows <- vrows & as.logical(rowSums(Z[,vids]))
              eig.L0 <- emma.eigen.L(Z[vtrows,vtids],K[vtids,vtids])
              MLE0 <- emma.MLE(ys[j,vtrows],X0[vtrows,,drop=FALSE],K[vtids,vtids],Z[vtrows,vtids],ngrids,llim,ulim,esp,eig.L0)
              MLE1 <- emma.MLE(ys[j,vtrows],X[vtrows,],K[vtids,vtids],Z[vtrows,vtids],ngrids,llim,ulim,esp,eig.L0)
            }
            if (!ponly) { 
              ML1s[i,j] <- MLE1$ML
              vgs[i,j] <- MLE1$vg
              ves[i,j] <- MLE1$ve
              ML0s[i,j] <- MLE0$ML
            }
            stats[i,j] <- 2*(MLE1$ML-MLE0$ML)
          }
        }
        if ( ( nv == t ) && ( !ponly ) ) {
          ML0s[i,] <- ML0
        }
        ps[i,] <- pchisq(stats[i,],1,lower.tail=FALSE)
      }
    }    
  }
  if ( ponly ) {
    return (ps)
  }
  else {
    return (list(ps=ps,ML1s=ML1s,ML0s=ML0s,stats=stats,vgs=vgs,ves=ves))
  }  
}

emma.REML.t <- function(ys, xs, K, Z=NULL, X0 = NULL, ngrids=100, llim=-10, ulim=10, esp=1e-10, ponly = FALSE) {
  if ( is.null(dim(ys)) || ncol(ys) == 1 ) {
    ys <- matrix(ys,1,length(ys))
  }
  if ( is.null(dim(xs)) || ncol(xs) == 1 ) {
    xs <- matrix(xs,1,length(xs))
  }
  if ( is.null(X0) ) {
    X0 <- matrix(1,ncol(ys),1)
  }
  
  g <- nrow(ys)
  n <- ncol(ys)
  m <- nrow(xs)
  t <- ncol(xs)
  q0 <- ncol(X0)
  q1 <- q0 + 1
  
  stopifnot(nrow(K) == t)
  stopifnot(ncol(K) == t)
  stopifnot(nrow(X0) == n)

  if ( !ponly ) {
    REMLs <- matrix(nrow=m,ncol=g)
    vgs <- matrix(nrow=m,ncol=g)
    ves <- matrix(nrow=m,ncol=g)
  }
  dfs <- matrix(nrow=m,ncol=g)
  stats <- matrix(nrow=m,ncol=g)
  ps <- matrix(nrow=m,ncol=g)

  if ( sum(is.na(ys)) == 0 ) {
    eig.L <- emma.eigen.L(Z,K)

    x.prev <- vector(length=0)

    for(i in 1:m) {
      vids <- !is.na(xs[i,])
      nv <- sum(vids)
      xv <- xs[i,vids]

      if ( ( mean(xv) <= 0 ) || ( mean(xv) >= 1 ) ) {
        if ( !ponly ) {
          vgs[i,] <- rep(NA,g)
          ves[i,] <- rep(NA,g)
          dfs[i,] <- rep(NA,g)
          REMLs[i,] <- rep(NA,g)
          stats[i,] <- rep(NA,g)
        }
        ps[i,] = rep(1,g)
        
      }
      else if ( identical(x.prev, xv) ) {
        if ( !ponly ) {
          vgs[i,] <- vgs[i-1,]
          ves[i,] <- ves[i-1,]
          dfs[i,] <- dfs[i-1,]
          REMLs[i,] <- REMLs[i-1,]
          stats[i,] <- stats[i-1,]
        }
        ps[i,] <- ps[i-1,]
      }
      else {
        if ( is.null(Z) ) {
          X <- cbind(X0[vids,,drop=FALSE],xs[i,vids])
          eig.R1 = emma.eigen.R.wo.Z(K[vids,vids],X)
        }
        else {
          vrows <- as.logical(rowSums(Z[,vids]))              
          X <- cbind(X0[vrows,,drop=FALSE],Z[vrows,vids,drop=FALSE]%*%t(xs[i,vids,drop=FALSE]))
          eig.R1 = emma.eigen.R.w.Z(Z[vrows,vids],K[vids,vids],X)
        }
        
        for(j in 1:g) {
          if ( nv == t ) {
            REMLE <- emma.REMLE(ys[j,],X,K,Z,ngrids,llim,ulim,esp,eig.R1)
            if ( is.null(Z) ) {
              U <- eig.L$vectors * matrix(sqrt(1/(eig.L$values+REMLE$delta)),t,t,byrow=TRUE)
              dfs[i,j] <- nv - q1
            }
            else {
              U <- eig.L$vectors * matrix(c(sqrt(1/(eig.L$values+REMLE$delta)),rep(sqrt(1/REMLE$delta),n-t)),n,n,byrow=TRUE)
              dfs[i,j] <- n - q1
            }
            yt <- crossprod(U,ys[j,])
            Xt <- crossprod(U,X)
            iXX <- solve(crossprod(Xt,Xt))
            beta <- iXX%*%crossprod(Xt,yt)
            
            if ( !ponly ) {
              vgs[i,j] <- REMLE$vg
              ves[i,j] <- REMLE$ve
              REMLs[i,j] <- REMLE$REML
            }
            stats[i,j] <- beta[q1]/sqrt(iXX[q1,q1]*REMLE$vg)
          }
          else {
            if ( is.null(Z) ) {
              eig.L0 <- emma.eigen.L.wo.Z(K[vids,vids])
              nr <- sum(vids)
              yv <- ys[j,vids]
              REMLE <- emma.REMLE(yv,X,K[vids,vids,drop=FALSE],NULL,ngrids,llim,ulim,esp,eig.R1)
              U <- eig.L0$vectors * matrix(sqrt(1/(eig.L0$values+REMLE$delta)),nr,nr,byrow=TRUE)
              dfs[i,j] <- nr - q1
            }
            else {
              eig.L0 <- emma.eigen.L.w.Z(Z[vrows,vids,drop=FALSE],K[vids,vids])              
              yv <- ys[j,vrows]
              nr <- sum(vrows)
              tv <- sum(vids)
              REMLE <- emma.REMLE(yv,X,K[vids,vids,drop=FALSE],Z[vrows,vids,drop=FALSE],ngrids,llim,ulim,esp,eig.R1)
              U <- eig.L0$vectors * matrix(c(sqrt(1/(eig.L0$values+REMLE$delta)),rep(sqrt(1/REMLE$delta),nr-tv)),nr,nr,byrow=TRUE)
              dfs[i,j] <- nr - q1
            }
            yt <- crossprod(U,yv)
            Xt <- crossprod(U,X)
            iXX <- solve(crossprod(Xt,Xt))
            beta <- iXX%*%crossprod(Xt,yt)
            if (!ponly) {
              vgs[i,j] <- REMLE$vg
              ves[i,j] <- REMLE$ve
              REMLs[i,j] <- REMLE$REML
            }
            stats[i,j] <- beta[q1]/sqrt(iXX[q1,q1]*REMLE$vg)
          }
        }
        ps[i,] <- 2*pt(abs(stats[i,]),dfs[i,],lower.tail=FALSE)
      }
    }
  }
  else {
    eig.L <- emma.eigen.L(Z,K)
    eig.R0 <- emma.eigen.R(Z,K,X0)
      
    x.prev <- vector(length=0)
    
    for(i in 1:m) {
      vids <- !is.na(xs[i,])
      nv <- sum(vids)
      xv <- xs[i,vids]

      if ( ( mean(xv) <= 0 ) || ( mean(xv) >= 1 ) ) {
        if (!ponly) {
          vgs[i,] <- rep(NA,g)
          ves[i,] <- rep(NA,g)
          REMLs[i,] <- rep(NA,g)
          dfs[i,] <- rep(NA,g)
        }
        ps[i,] = rep(1,g)
      }      
      else if ( identical(x.prev, xv) ) {
        if ( !ponly ) {
          stats[i,] <- stats[i-1,]
          vgs[i,] <- vgs[i-1,]
          ves[i,] <- ves[i-1,]
          REMLs[i,] <- REMLs[i-1,]
          dfs[i,] <- dfs[i-1,]
        }
        ps[i,] = ps[i-1,]
      }
      else {
        if ( is.null(Z) ) {
          X <- cbind(X0,xs[i,])
          if ( nv == t ) {
            eig.R1 = emma.eigen.R.wo.Z(K,X)
          }
        }
        else {
          vrows <- as.logical(rowSums(Z[,vids,drop=FALSE]))
          X <- cbind(X0,Z[,vids,drop=FALSE]%*%t(xs[i,vids,drop=FALSE]))
          if ( nv == t ) {
            eig.R1 = emma.eigen.R.w.Z(Z,K,X)
          }          
        }

        for(j in 1:g) {
          vrows <- !is.na(ys[j,])
          if ( nv == t ) {
            yv <- ys[j,vrows]
            nr <- sum(vrows)
            if ( is.null(Z) ) {
              if ( nr == n ) {
                REMLE <- emma.REMLE(yv,X,K,NULL,ngrids,llim,ulim,esp,eig.R1)
                U <- eig.L$vectors * matrix(sqrt(1/(eig.L$values+REMLE$delta)),n,n,byrow=TRUE)                
              }
              else {
                eig.L0 <- emma.eigen.L.wo.Z(K[vrows,vrows,drop=FALSE])
                REMLE <- emma.REMLE(yv,X[vrows,,drop=FALSE],K[vrows,vrows,drop=FALSE],NULL,ngrids,llim,ulim,esp)
                U <- eig.L0$vectors * matrix(sqrt(1/(eig.L0$values+REMLE$delta)),nr,nr,byrow=TRUE)
              }
              dfs[i,j] <- nr-q1
            }
            else {
              if ( nr == n ) {
                REMLE <- emma.REMLE(yv,X,K,Z,ngrids,llim,ulim,esp,eig.R1)
                U <- eig.L$vectors * matrix(c(sqrt(1/(eig.L$values+REMLE$delta)),rep(sqrt(1/REMLE$delta),n-t)),n,n,byrow=TRUE)                
              }
              else {
                vtids <- as.logical(colSums(Z[vrows,,drop=FALSE]))
                eig.L0 <- emma.eigen.L.w.Z(Z[vrows,vtids,drop=FALSE],K[vtids,vtids,drop=FALSE])
                REMLE <- emma.REMLE(yv,X[vrows,,drop=FALSE],K[vtids,vtids,drop=FALSE],Z[vrows,vtids,drop=FALSE],ngrids,llim,ulim,esp)
                U <- eig.L0$vectors * matrix(c(sqrt(1/(eig.L0$values+REMLE$delta)),rep(sqrt(1/REMLE$delta),nr-sum(vtids))),nr,nr,byrow=TRUE)
              }
              dfs[i,j] <- nr-q1
            }

            yt <- crossprod(U,yv)
            Xt <- crossprod(U,X[vrows,,drop=FALSE])
            iXX <- solve(crossprod(Xt,Xt))
            beta <- iXX%*%crossprod(Xt,yt)
            if ( !ponly ) {
              vgs[i,j] <- REMLE$vg
              ves[i,j] <- REMLE$ve
              REMLs[i,j] <- REMLE$REML
            }
            stats[i,j] <- beta[q1]/sqrt(iXX[q1,q1]*REMLE$vg)
          }
          else {
            if ( is.null(Z) ) {
              vtids <- vrows & vids
              eig.L0 <- emma.eigen.L.wo.Z(K[vtids,vtids,drop=FALSE])
              yv <- ys[j,vtids]
              nr <- sum(vtids)
              REMLE <- emma.REMLE(yv,X[vtids,,drop=FALSE],K[vtids,vtids,drop=FALSE],NULL,ngrids,llim,ulim,esp)
              U <- eig.L0$vectors * matrix(sqrt(1/(eig.L0$values+REMLE$delta)),nr,nr,byrow=TRUE)
              Xt <- crossprod(U,X[vtids,,drop=FALSE])
              dfs[i,j] <- nr-q1
            }
            else {
              vtids <- as.logical(colSums(Z[vrows,,drop=FALSE])) & vids
              vtrows <- vrows & as.logical(rowSums(Z[,vids,drop=FALSE]))
              eig.L0 <- emma.eigen.L.w.Z(Z[vtrows,vtids,drop=FALSE],K[vtids,vtids,drop=FALSE])
              yv <- ys[j,vtrows]
              nr <- sum(vtrows)
              REMLE <- emma.REMLE(yv,X[vtrows,,drop=FALSE],K[vtids,vtids,drop=FALSE],Z[vtrows,vtids,drop=FALSE],ngrids,llim,ulim,esp)
              U <- eig.L0$vectors * matrix(c(sqrt(1/(eig.L0$values+REMLE$delta)),rep(sqrt(1/REMLE$delta),nr-sum(vtids))),nr,nr,byrow=TRUE)
              Xt <- crossprod(U,X[vtrows,,drop=FALSE])
              dfs[i,j] <- nr-q1
            }
            yt <- crossprod(U,yv)
            iXX <- solve(crossprod(Xt,Xt))
            beta <- iXX%*%crossprod(Xt,yt)
            if ( !ponly ) {
              vgs[i,j] <- REMLE$vg
              ves[i,j] <- REMLE$ve
              REMLs[i,j] <- REMLE$REML
            }
            stats[i,j] <- beta[q1]/sqrt(iXX[q1,q1]*REMLE$vg)
            
          }
        }
        ps[i,] <- 2*pt(abs(stats[i,]),dfs[i,],lower.tail=FALSE)        
      }
    }    
  }
  if ( ponly ) {
    return (ps)
  }
  else {
    return (list(ps=ps,REMLs=REMLs,stats=stats,dfs=dfs,vgs=vgs,ves=ves))
  }
}

