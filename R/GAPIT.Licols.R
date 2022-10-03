
`GAPIT.Licols` <- function(X,tol=1e-10){
# Extract a linearly independent set of columns of a given matrix X
#
#    [Xsub,idx]=licols(X,tol)
#
## Input:
#  X: The given input matrix
#  tol: A rank estimation tolerance. Default=1e-10
#
## Output:
#  Xsub: The extracted columns of X
#  idx:  The indices (into X) of the extracted columns 
#Authors: Jiabo Wang
#Writer:  Li Chen and Jiabo Wang
# Last update: MAY 12, 2022 
##############################################################################################
  if (all(X==0)){     # X is a zero matrix        
    idx <- Xsub <- c()
  }else{              # X is not a zero matrix 
    qr_res <- qr(X,LAPACK=F) # QR decomposition 
    Q <- qr.Q(qr_res)
    R <- qr.R(qr_res)
    E <- qr_res$pivot
    
    if (is.vector(R) == 0){
      diagr <- abs(diag(R))
    }else{
      diagr <- abs(R[1])
    }
    
    # Rank estimation
    r <- sum(diagr >= tol * diagr[1])
    
    idx <- sort(E[1:r])
    Xsub <- X[,idx,drop=FALSE]
  }
  res <- vector("list")
  res$Xsub <- Xsub
  res$idx <- idx
  return(res)
}
