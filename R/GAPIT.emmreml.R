`emmreml`<-function (y, X, Z, K, varbetahat = FALSE, varuhat = FALSE, PEVuhat = FALSE, 
    test = FALSE) 
{
    q = dim(X)[2]
    n = length(y)
    spI <- diag(n)
    try.so=try(solve(crossprod(as.matrix(t(X)))),silent=TRUE)
    try.xx=TRUE
    if(inherits(try.so, "try-error"))
        {try.xx=FALSE}
    if(try.xx)
    {S <- spI - tcrossprod(X %*% solve(crossprod(X)), X)
    }else{
     S <- spI - tcrossprod(X %*% MASS::ginv(crossprod(X)), X)
    }
    ZK <- Z %*% K
    offset <- log(n)
    ZKZt <- tcrossprod(ZK, Z)
    ZKZtandoffset <- ZKZt + offset * spI
    SZKZtSandoffset <- {
        S %*% ZKZtandoffset
    } %*% S
    svdSZKZtSandspI <- eigen(SZKZtSandoffset, symmetric = TRUE)
    Ur <- svdSZKZtSandspI$vectors[, 1:(n - q)]
    lambda <- svdSZKZtSandspI$values[1:(n - q)] - offset
    eta <- crossprod(Ur, y)
    minimfunc <- function(delta) {
        (n - q) * log(sum(eta^2/{
            lambda + delta
        })) + sum(log(lambda + delta))
    }
    optimout <- stats::optimize(minimfunc, lower = 9^(-9), upper = 9^9, 
        tol = 1e-06)
    deltahat <- optimout$minimum
    # Hinvhat <- solve(ZKZt + deltahat * spI)
    if(try.xx)
    {Hinvhat <- solve(ZKZt + deltahat * spI)
    }else{
     Hinvhat <- MASS::ginv(ZKZt + deltahat * spI)
    }
    XtHinvhat <- crossprod(X, Hinvhat)

    # betahat <- solve(XtHinvhat %*% X, XtHinvhat %*% y)
    if(try.xx)
    {betahat <- solve(XtHinvhat %*% X, XtHinvhat %*% y)
    }else{
     betahat <- MASS::ginv(XtHinvhat %*% X, XtHinvhat %*% y)
    }
    
    ehat <- (y - {
        X %*% betahat
    })
    Hinvhatehat <- Hinvhat %*% ehat
    sigmausqhat <- sum(eta^2/{
        lambda + deltahat
    })/(n - q)
    Vinv <- (1/sigmausqhat) * Hinvhat
    sigmaesqhat <- deltahat * sigmausqhat
    uhat <- crossprod(ZK, Hinvhatehat)
    df <- n - q
    loglik <- -0.5 * (optimout$objective + df + df * log(2 * 
        pi/df))
    if (varuhat) {
        # P <- Vinv - Vinv %*% X %*% solve(crossprod(X, Vinv %*% 
        #     X), crossprod(X, Vinv))
        if(try.xx)
    {P <- Vinv - Vinv %*% X %*% solve(crossprod(X, Vinv %*% 
            X), crossprod(X, Vinv))
    }else{
     P <- Vinv - Vinv %*% X %*% MASS::ginv(crossprod(X, Vinv %*% 
            X), crossprod(X, Vinv))
    }
        varuhat <- sigmausqhat^2 * crossprod(ZK, P) %*% ZK
    }
    if (PEVuhat) {
        if (!exists("P")) {
            # P <- Vinv - Vinv %*% X %*% solve(crossprod(X, Vinv %*% 
            #     X), crossprod(X, Vinv))
            if(try.xx)
    {P <- Vinv - Vinv %*% X %*% solve(crossprod(X, Vinv %*% 
            X), crossprod(X, Vinv))
    }else{
     P <- Vinv - Vinv %*% X %*% MASS::ginv(crossprod(X, Vinv %*% 
            X), crossprod(X, Vinv))
    }
        }
        PEVuhat <- sigmausqhat * K - varuhat
    }
    if (varbetahat) {
        # varbetahat <- solve(crossprod(X, Vinv %*% X))
        if(try.xx)
    {varbetahat <- solve(crossprod(X, Vinv %*% X))
    }else{
     varbetahat <- MASS::ginv(crossprod(X, Vinv %*% X))
    }
    }
    if (test) {
        Xsqtestu <- uhat^2/diag(varuhat)
        puhat <- stats::pchisq(Xsqtestu, df = 1, lower.tail = F, log.p = F)
        p.adjust.M <- stats::p.adjust.methods
        p.adjuhat <- sapply(p.adjust.M, function(meth) stats::p.adjust(puhat, 
            meth))
        Xsqtestbeta <- betahat^2/diag(varbetahat)
        pbetahat <- stats::pchisq(Xsqtestbeta, df = 1, lower.tail = F, 
            log.p = F)
        p.adjbetahat <- sapply(p.adjust.M, function(meth) stats::p.adjust(pbetahat, 
            meth))
    }
    if (!exists("Xsqtestbeta")) {
        Xsqtestbeta <- c()
    }
    if (!exists("pvalbeta")) {
        pvalbeta <- c()
    }
    if (!exists("Xsqtestu")) {
        Xsqtestu <- c()
    }
    if (!exists("p.adjuhat")) {
        p.adjuhat <- c()
    }
    if (!exists("p.adjbetahat")) {
        p.adjbetahat <- c()
    }
    if (!exists("varuhat")) {
        varuhat <- c()
    }
    if (!exists("varbeta")) {
        varubeta <- c()
    }
    if (!exists("PEVuhat")) {
        PEVuhat <- c()
    }
    return(list(Vu = sigmausqhat, Ve = sigmaesqhat, betahat = betahat, 
        uhat = uhat, Xsqtestbeta = Xsqtestbeta, pvalbeta = p.adjbetahat, 
        Xsqtestu = Xsqtestu, pvalu = p.adjuhat, varuhat = diag(varuhat), 
        varbetahat = diag(varbetahat), PEVuhat = diag(PEVuhat), 
        loglik = loglik))
}
