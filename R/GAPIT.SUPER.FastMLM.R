`GAPIT.SUPER.FastMLM` <-
function(ys, xs, vg, delta, Z = NULL, X0 = NULL, snp.pool=NULL,LD=0.01,method="FaST") {
#Input: ys, xs, vg, delta, Z, X0, snp.pool
#Output: GWAS
#Authors: Qishan Wang, Feng Tian and Zhiwu Zhang
#Last update: April 16, 2012
################################################################################
#print("GAPIT.SUPER.FastMLM started")
#print("dimension of ys,xs,X0 and snp.pool")
#print(length(ys))
#print(dim(xs))
#print(dim(X0))
#print(dim(snp.pool))
#print((LD))


#Set data to the require format
ys=unlist(ys)
if(is.null(dim(ys)) || ncol(ys) == 1)  ys <- matrix(ys, 1, length(ys))
if(is.null(dim(xs)) || ncol(xs) == 1)  xs <- matrix(xs, 1, length(xs))
if(is.null(X0))  X0 <- matrix(1, nrow(snp.pool), 1)

#Exract data size
g <- nrow(ys)
n <- nrow(xs)   #####  generaol nrow(xs)=nrow(U1) rivised by qishan 2012.4.16
m <- ncol(xs)
t <- nrow(xs)
q0 <- ncol(X0)
q1 <- q0 + 1

#Allocate space
dfs <- matrix(nrow = m, ncol = g)
stats <- matrix(nrow = m, ncol = g)
ps <- matrix(nrow = m, ncol = g)
betavalue <- matrix(nrow = m, ncol = g)
####################
if(method=="SUPER"){
 LDsqr=sqrt(LD)  
##################  
 
#Iteration on trait (j) and SNP (i)
for(j in 1:g)
{
 
for (i in 1:m)
{
  if((i >0)&(floor(i/500)==i/500))  print(paste("SNP: ",i," ",sep=""))


  #No variation on the SNP
  if(min(xs[,i])==max(xs[,i]))
  {
    dfs[i,j] <- n - q1
    betavalue[i,j]=0
    stats[i,j] <- 0
  }
  #The SNP has variation
  if(min(xs[,i])!=max(xs[,i]))
  {
      #SUPER
      snp.corr = stats::cor(xs[,i],snp.pool)
      index.k=which( abs(snp.corr)<=LDsqr )
      #handler of snp correlated with all QTNs
      if(length(index.k)<2){
       index.k=1:length(snp.corr) #keep going to have them all
       #print("warning: there is a snp correlated with all QTNs")
      }   
      K.X= snp.pool[,index.k]
      ####################
      K.X.svd= svd(K.X) ###start 2012.4.16 by qishan
  
       d=K.X.svd$d
       d=d[d>1e-8]
       d=d^2

       U1=K.X.svd$u   
       U1=U1[,1:length(d)]  ### end 2012.4.16 by qishan
 
       n<-nrow(U1)

      
       I= diag(1,nrow(U1))
      
      ################ get iXX
         X <- cbind(X0, xs[,i]) ####marker by column
         U <- U1*matrix(sqrt(1/(d + delta)), nrow(U1), length(d), byrow = TRUE) 
         Xt <- crossprod(U, X) 
         XX1<- crossprod(Xt, Xt)
         XX2<- crossprod((I-tcrossprod(U1,U1))%*%X,(I-tcrossprod(U1,U1))%*%X)/delta
         #iXX<-solve(XX1+XX2) 
         
           iXX <- try(solve(XX1+XX2),silent=T)
     if(inherits(iXX, "try-error")){
     iXX <- MASS::ginv(XX1+XX2)
     }
      #################  end get ixx
      ################   begin get beta
      ################
    #######get beta compnents 1
#U1TX=t(U1)%*%X
U1TX=crossprod(U1,X)
beta1=0
for(ii in 1:length(d)){
one=matrix(U1TX[ii,], nrow=1)
dim(one)
#beta=t(one)%*%one/(d[ii]+delta)
beta=crossprod(one,one)/(d[ii]+delta)
beta1= beta1+beta
}

#######get beta components 2
#IUX=(I-U1%*%t(U1))%*%X
IUX=(I-tcrossprod(U1,U1))%*%X
beta2=0
for(ii in 1:nrow(U1)){
one=matrix(IUX[ii,], nrow=1)
dim(one)
beta=t(one)%*%one
beta2= beta2+beta
}
beta2<-beta2/delta

#######get b3
#U1TY=t(U1)%*%ys[j,]
U1TY=crossprod(U1,ys[j,])
beta3=0
for(ii in 1:length(d)){
one1=matrix(U1TX[ii,], nrow=1)
one2=matrix(U1TY[ii,], nrow=1)
beta=crossprod(one1,one2)/(d[ii]+delta)
beta3= beta3+beta
}

###########get beta4
#IUY=(I-U1%*%t(U1))%*%ys[j,]
IUY=(I-tcrossprod(U1,U1))%*%ys[j,]
beta4=0
for(ii in 1:nrow(U1)){
one1=matrix(IUX[ii,], nrow=1)
one2=matrix(IUY[ii,], nrow=1)
#beta=t(one1)%*%one2
beta=crossprod(one1,one2)
beta4= beta4+beta
}
beta4<-beta4/delta

#######get final beta
beta = MASS::ginv(beta1+beta2)%*%(beta3+beta4)
   
      ##############
      ################    end get beta

    betavalue[i,j]=beta[q1,1]
    stats[i,j] <- beta[q1,1]/sqrt(iXX[q1, q1] * vg)
    dfs[i,j] <- n - q1
  } #end of SNP variation stutus detection
} #loop for markers

#print("Calculating p-values...")
ps[,j] <- 2 * stats::pt(abs(stats[,j]), dfs[,j],  lower.tail = FALSE)
} #end of loop on traits

return(list(beta=betavalue, ps = ps, stats = stats, dfs = dfs,effect=betavalue))
} #Enf of SUPERMLM

#######################
if(method=="FaST"){
 K.X.svd= svd(snp.pool) ###start 2012.4.16 by qishan
  
       d=K.X.svd$d
       d=d[d>1e-8]
       d=d^2

       U1=K.X.svd$u   
       U1=U1[,1:length(d)]  ### end 2012.4.16 by qishan
 
       n<-nrow(U1)
       I= diag(1,nrow(U1))
   U <- U1*matrix(sqrt(1/(d + delta)), nrow(U1), length(d), byrow = TRUE) 
################## 
 
#Iteration on trait (j) and SNP (i)
for(j in 1:g)
{
 
for (i in 1:m)
{
  if((i >0)&(floor(i/500)==i/500))  print(paste("SNP: ",i," ",sep=""))


  #No variation on the SNP
  if(min(xs[,i])==max(xs[,i]))
  {
    dfs[i,j] <- n - q1
    betavalue[i,j]=0
    stats[i,j] <- 0
  }
  #The SNP has variation
  if(min(xs[,i])!=max(xs[,i]))
  {
      #SUPER
      
      ####################
      K.X.svd= svd(snp.pool) ###start 2012.4.16 by qishan
  
       d=K.X.svd$d
       d=d[d>1e-8]
       d=d^2

       U1=K.X.svd$u   
       U1=U1[,1:length(d)]  ### end 2012.4.16 by qishan
 
       n<-nrow(U1)
       I= diag(1,nrow(U1))
      
      ################ get iXX
         X <- cbind(X0, xs[,i]) ####marker by column
         U <- U1*matrix(sqrt(1/(d + delta)), nrow(U1), length(d), byrow = TRUE) 
         Xt <- crossprod(U, X) 
         XX1<- crossprod(Xt, Xt)
         XX2<- crossprod((I-tcrossprod(U1,U1))%*%X,(I-tcrossprod(U1,U1))%*%X)/delta
                iXX <- try(solve(XX1+XX2),silent=T)
     if(inherits(iXX, "try-error")){
     iXX <- MASS::ginv(XX1+XX2)
     }
      #################  end get ixx
      ################   begin get beta
    #######get beta compnents 1
#U1TX=t(U1)%*%X
U1TX=crossprod(U1,X)
beta1=0
for(ii in 1:length(d)){
one=matrix(U1TX[ii,], nrow=1)
dim(one)
beta=crossprod(one,one)/(d[ii]+delta)
beta1= beta1+beta
}

#######get beta components 2
IUX=(I-tcrossprod(U1,U1))%*%X
beta2=0
for(ii in 1:nrow(U1)){
one=matrix(IUX[ii,], nrow=1)
dim(one)
beta=crossprod(one,one)
beta2= beta2+beta
}
beta2<-beta2/delta

#######get b3
#U1TY=t(U1)%*%ys[j,]
U1TY=crossprod(U1,ys[j,])
beta3=0
for(ii in 1:length(d)){
one1=matrix(U1TX[ii,], nrow=1)
one2=matrix(U1TY[ii,], nrow=1)
#beta=t(one1)%*%one2/(d[ii]+delta)
beta=crossprod(one1,one2)/(d[ii]+delta)
beta3= beta3+beta
}

###########get beta4
#IUY=(I-U1%*%t(U1))%*%ys[j,]
IUY=(I-tcrossprod(U1,U1))%*%ys[j,]
beta4=0
for(ii in 1:nrow(U1)){
one1=matrix(IUX[ii,], nrow=1)
one2=matrix(IUY[ii,], nrow=1)
#beta=t(one1)%*%one2
beta=crossprod(one1,one2)
beta4= beta4+beta
}
beta4<-beta4/delta

#######get final beta
beta = MASS::ginv(beta1+beta2)%*%(beta3+beta4)
   
      ##############
      ################    end get beta

    betavalue[i,j]=beta[q1,1]
    stats[i,j] <- beta[q1,1]/sqrt(iXX[q1, q1] * vg)
    dfs[i,j] <- n - q1
  } #end of SNP variation stutus detection
} #loop for markers

#print("Calculating p-values...")
ps[,j] <- 2 * stats::pt(abs(stats[,j]), dfs[,j],  lower.tail = FALSE)
} #end of loop on traits

return(list(beta=betavalue, ps = ps, stats = stats, dfs = dfs,effect=betavalue))
} #Enf of FastMLM

}####end function
#=============================================================================================

