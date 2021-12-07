`GAPIT.ROC` <-
function(t=NULL,se=NULL,Vp=1,trait="",plot.style="rainbow"){
    #Object: To make table and plot for ROC (power vs FDR)
    #Input: t and se are the vectors of t value and their standard error
    #Input: Vp is phenotypic variance and trait is name of the phenotype
    #Output: A table and plot
    #Requirment: error df is same for all SMP or large
    #Authors: Zhiwu Zhang
    # Last update: Feb 11, 2013
    ##############################################################################################
#print("GAPIT.ROC start")
#print("Length of t se and Vp")
#print(length(t))
# aa=read.csv("GAPIT.MLM.MLM.V1.Df.tValue.StdErr.csv",head=T)
# #print(length(se))
# t=aa$t.Value
# se=aa$std.Error
# Vp=var(mySim$Y[,2])
# trait="V1"
#print((Vp))
if(length(t)==length(t[is.na(t)]) ){
#print("NA t, No ROC plot")
return(NULL)
}
    
    #test
    #n=1000
    #trait="test"
    #t=rnorm(n)
    #se=sqrt(abs(rnorm(n))  )
    #Vp=10
    
    #Remove NAs
    index=is.na(t)
    t=t[!index]
    se=se[!index]
    #print(head(cbind(t,se)))
    #Configration
    FDR=c(0,.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1)
    coefficient=c(0,0.01,.02,.05,.1,.2,.3)
    
    #Power holder
    nf=length(FDR)
    nc=length(coefficient)
    power=matrix(NA,nf,nc)
    
    #Handler of matrix format
    if(!is.null(dim(t))) t=t[,1]
    if(!is.null(dim(se))) se=se[,1]
    
    n=length(t)
    
    #Discard negative
    t=abs(t)
    #print("@@@@@@@@@@@@@@")
    #sort t and se
    position=order(t,decreasing = TRUE)
    t=t[position]
    se=se[position]
    EFFECT=coefficient*sqrt(Vp)
    newbit=matrix(1/se,n,1)%*%EFFECT   #n by nc matrix
    tnew=newbit+t  #n by nc matrix
    
    for (i in 1:nf){
        fdr=FDR[i]
        cutpoint=floor(n*fdr)
        cutoff=t[cutpoint]
        
        
        for (j in 1:nc){
            effect= EFFECT[j]
            singnificant=tnew[,j]>cutoff
            count=length(t[singnificant])
            power[i,j]=count/n
            
        } #end of for on fdr
    } #end of for on effect
    
    #output
    rownames(power)=FDR
    tkk<-c(.3,.2,.1,.05,.02,0.01,0)
    tc1<-c(0,0.25,0.5,0.75,1.0)
    #colnames(power)=paste("QTN=",coefficient,sep="")
    colnames(power)=paste("QTN=",tkk,sep="")

    if(plot.style=="FarmCPU"){
    utils::write.table(power,file=paste("FarmCPU.",trait,".ROC.csv",sep=""),quote = TRUE, sep = ",", row.names = TRUE,col.names = NA)
    }
    if(plot.style=="rainbow"){
        utils::write.table(power,file=paste("GAPIT.",trait,".ROC.csv",sep=""),quote = TRUE, sep = ",", row.names = TRUE,col.names = NA)
    }
    FDR_log<-FDR/10
    #palette(c("black","red","blue","brown", "orange","cyan", "green",rainbow(nc)))
    if(plot.style=="FarmCPU"){
    grDevices::pdf(paste("FarmCPU.", trait,".ROC.pdf" ,sep = ""), width = 5,height=5)
    graphics::par(mar = c(5,6,5,3))
    }
    if(plot.style=="rainbow"){
        grDevices::pdf(paste("GAPIT.", trait,".ROC.pdf" ,sep = ""), width = 7,height=7)
        graphics::par(mar = c(5,5,5,3))
    }
  
 grDevices::palette(c("black","red","blue","brown", "orange","cyan", "green", grDevices::rainbow(nc)))
    plot(FDR_log,power[,1],log="x",type="o",yaxt="n",lwd=3,col=1,xlab="Type I error",ylab="Power",main = trait,cex.axis=1.3, cex.lab=1.3)
    graphics::axis(side=2,at=tc1,labels=tc1,cex.lab=1.3,cex.axis=1.3)
    for(i in 2:nc){
        graphics::lines(power[,i]~FDR_log, lwd=3,type="o",pch=i,col=i)
    }
    #legend("bottomright", colnames(power), pch = c(1:nc), lty = c(1,2),col=c(1:nc))
   graphics::legend("bottomright", colnames(power), pch = c(nc:1), lty = c(1,2),col=c(nc:1),lwd=2,bty="n")
    grDevices::palette("default")      # reset back to the default
    #print("@@@@@@@@@@@@@@")
    #print(power)
    grDevices::dev.off()
print("ROC completed!")
    
}   #GAPIT.ROC ends here
#=============================================================================================

