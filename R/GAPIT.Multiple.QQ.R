`GAPIT.Multiple.QQ`<-function(mapP,DPP=5000,cutOff=0.01, wd=2, ratio=1,
    allpch=NULL,memo=NULL)
    #Object: Make a Multiple QQ Plot
    #Output: pdfs of the Multiple Manhattan Plot
    #Authors: Jiabo Wang
    # Last update: JUN 22, 2022
    ##############################################################################################
{
 Nenviron=ncol(mapP)-3
 allpch0=c(0,1,2,5,6)
 add.pch=c("+","*","-","#","<",">","^","$","=","|","?",as.character(1:9),letters[1:26],LETTERS[1:26]) 
 n.vals=ceiling(Nenviron/length(allpch0))-1
 # s=size-wd/ratio/2
 # DPP=500
 P=mapP[,-c(1:3)]
 values=apply(P,1,min)
 values=-log10(values)
 cut0=ceiling(-log10(0.01/length(values))/2)
 rv=runif(length(values))
 values=values+rv*(values-5+cut0)
 index=values>cut0
 print(table(index))
 # if(nrow(mapP)>DPP) mapP=mapP[index,]
 
 P=mapP[,-c(1:3)]
 taxa=colnames(mapP)[-c(1:3)]
 if(!is.null(memo))taxa=paste(taxa,"_",memo,sep="")
 if(Nenviron>5)
 {
  yourpch=c(rep(allpch0,n.vals),allpch0[1:(Nenviron-length(allpch0)*n.vals)])
 }else{
  yourpch=allpch0[1:Nenviron]
 }


 NN=nrow(P)
 themax.y0=max(-log10(P))
 pdf(paste("GAPIT.QQ.Mutiple.Plot.",memo,".symphysic",".pdf" ,sep = ""), width = 30,height=18)
 par(mfrow=c(1,1))
 par(mar = c(5,5,5,1))
 par(cex=1.8)
 themax.y02=ceiling((ceiling(themax.y0/4))*4)
 p_value_quantiles0 <- (1:NN)/(NN+1)
 log.Quantiles0 <- -log10(p_value_quantiles0)

 N1=NN
        ## create the confidence intervals
 c95 <- rep(NA,N1)
 c05 <- rep(NA,N1)
 for(j in 1:N1)
 {
    i=ceiling((10^-log.Quantiles0[j])*NN)
    if(i==0)i=1
    c95[j] <- stats::qbeta(0.95,i,NN-i+1)
    c05[j] <- stats::qbeta(0.05,i,NN-i+1)
            #print(c(j,i,c95[j],c05[j]))
 }

 plot(NULL, xlim = c(0,max(log.Quantiles0)), ylim = c(0,themax.y0), 
 	type="l",lty=5, lwd = 2, las=1,
 	ylab=expression(Observed~~-log[10](italic(p))), xlab=expression(Expected~~-log[10](italic(p))),
 	col="gray")
 index=length(c95):1        
 graphics::polygon(c(log.Quantiles0[index],log.Quantiles0),c(-log10(c05)[index],-log10(c95)),col='gray',border=NA)
 graphics::abline(a = 0, b = 1, col = "red",lwd=2)
 step.vals0=NULL  
 for(i in 1:Nenviron)
 {
 	P.values=P[,i]
    P.values=P.values[P.values>0]
    P.values=P.values[P.values<=1]
    
    P.values <- P.values[order(P.values)]
    step.vals=ceiling(i/length(allpch0))-1
    mypch=allpch0[i-step.vals*length(allpch0)]

    #Set up the p-value quantiles
    #print("Setting p_value_quantiles...")
    p_value_quantiles <- (1:NN)/(NN+1)
    log.P.values <- -log10(P.values)
    log.Quantiles <- -log10(p_value_quantiles)
    # log.P.values2=apply(cbind(log.P.values,log.Quantiles),1,mean)
    # print(max(log.P.values))
    # log.P.values[5:NN]=log.P.values2[5:NN]
    if(nrow(mapP)>DPP) log.P.values=log.P.values[index]
    if(nrow(mapP)>DPP) log.Quantiles=log.Quantiles[index]

    # log.P.values=log.P.values[order(log.P.values)]
    # print(themax.y0)
    # print(max(log.P.values2))
    # print(max(log.P.values))
        par(new=T)
        plot(log.Quantiles, log.P.values, xlim = c(0,max(log.Quantiles0)), 
        	ylim = c(0,themax.y0), cex.axis=1, cex.lab=1.3, axes=FALSE, 
        	lty = 1,  lwd = 2, col = step.vals+1 ,xlab ="",
        	ylab ="", pch=mypch,
        	)
    #     if(step.vals!=0)
    # {
    #   points(log.Quantiles, log.P.values,pch=add.pch[step.vals],col="Blue",cex=1,cex.main=4)
    # }
    step.vals0=append(step.vals0,step.vals)

 }# end of Nenviron

graphics::legend("topleft",taxa,pch=yourpch,col=step.vals0+1,pt.lwd=3,text.font=6,box.col=NA)

grDevices::dev.off()



}# end of function



