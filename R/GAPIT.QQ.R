`GAPIT.QQ` <-
function(P.values, plot.type = "log_P_values", name.of.trait = "Trait",DPP=50000,plot.style="rainbow"){
    #Object: Make a QQ-Plot of the P-values
    #Options for plot.type = "log_P_values" and "P_values"
    #Output: A pdf of the QQ-plot
    #Authors: Alex Lipka and Zhiwu Zhang
    # Last update: May 9, 2011
    ##############################################################################################
    # Sort the data by the raw P-values
    #print("Sorting p values")
    #print(paste("Number of P values: ",length(P.values)))
    #remove NAs and keep the ones between between 0 and 1
    P.values=P.values[!is.na(P.values)]
    P.values=P.values[P.values>0]
    P.values=P.values[P.values<=1]
    
    if(length(P.values[P.values>0])<1) return(NULL)
    N=length(P.values)
    DPP=round(DPP/4) #Reduce to 1/4 for QQ plot
    P.values <- P.values[order(P.values)]
    
    #Set up the p-value quantiles
    #print("Setting p_value_quantiles...")
    p_value_quantiles <- (1:length(P.values))/(length(P.values)+1)
    
    
    if(plot.type == "log_P_values")
    {
        log.P.values <- -log10(P.values)
        log.Quantiles <- -log10(p_value_quantiles)
        
        index=GAPIT.Pruning(log.P.values,DPP=DPP)
        log.P.values=log.P.values[index ]
        log.Quantiles=log.Quantiles[index]
        
        if(plot.style=="FarmCPU"){
        pdf(paste("FarmCPU.", name.of.trait,".QQ-Plot.pdf" ,sep = ""),width = 5,height=5)
        par(mar = c(5,6,5,3))
        }
        if(plot.style=="rainbow"){
            pdf(paste("GAPIT.", name.of.trait,".QQ-Plot.pdf" ,sep = ""),width = 5,height=5)
            par(mar = c(5,6,5,3))
        }
        #Add conficence interval
        N1=length(log.Quantiles)
        ## create the confidence intervals
        c95 <- rep(NA,N1)
        c05 <- rep(NA,N1)
        for(j in 1:N1){
            i=ceiling((10^-log.Quantiles[j])*N)
            if(i==0)i=1
            c95[j] <- qbeta(0.95,i,N-i+1)
            c05[j] <- qbeta(0.05,i,N-i+1)
            #print(c(j,i,c95[j],c05[j]))
        }
        
        #CI Lines
        #plot(log.Quantiles, -log10(c05), xlim = c(0,max(log.Quantiles)), ylim = c(0,max(log.P.values)), type="l",lty=5, axes=FALSE, xlab="", ylab="",col="black")
        #par(new=T)
        #plot(log.Quantiles, -log10(c95), xlim = c(0,max(log.Quantiles)), ylim = c(0,max(log.P.values)), type="l",lty=5, axes=FALSE, xlab="", ylab="",col="black")
        
        #CI shade
        plot(NULL, xlim = c(0,max(log.Quantiles)), ylim = c(0,max(log.P.values)), type="l",lty=5, lwd = 2, axes=FALSE, xlab="", ylab="",col="gray")
        index=length(c95):1
        polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col='gray',border=NA)
        
        #Diagonal line
        abline(a = 0, b = 1, col = "red",lwd=2)
        
        #data
        par(new=T)
        if(plot.style=="FarmCPU"){
            plot(log.Quantiles, log.P.values, cex.axis=1.1, cex.lab=1.3, lty = 1,  lwd = 2, col = "Black" ,bty='l', xlab =expression(Expected~~-log[10](italic(p))), ylab = expression(Observed~~-log[10](italic(p))), main = paste(name.of.trait,sep=""),pch=20)
        }
        if(plot.style=="rainbow"){
            plot(log.Quantiles, log.P.values, xlim = c(0,max(log.Quantiles)), ylim = c(0,max(log.P.values)), cex.axis=1.1, cex.lab=1.3, lty = 1,  lwd = 2, col = "Blue" ,xlab =expression(Expected~~-log[10](italic(p))),ylab = expression(Observed~~-log[10](italic(p))), main = paste(name.of.trait,sep=""))
        }
        
        dev.off()
    }
    
    
    if(plot.type == "P_values")
    {
        pdf(paste("QQ-Plot_", name.of.trait,".pdf" ,sep = ""))
        par(mar = c(5,5,5,5))
        qqplot(p_value_quantiles, P.values, xlim = c(0,1),
        ylim = c(0,1), type = "l" , xlab = "Uniform[0,1] Theoretical Quantiles", 
        lty = 1, lwd = 1, ylab = "Quantiles of P-values from GWAS", col = "Blue",
        main = paste(name.of.trait,sep=" "))
        abline(a = 0, b = 1, col = "red")
        dev.off()   
    }
    #print("GAPIT.QQ  accomplished successfully!")
}
#=============================================================================================

