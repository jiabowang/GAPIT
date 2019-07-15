`GAPIT.MAF` <-
function(MAF=NULL,P=NULL,E=NULL,trait="",threshold.output=.1,plot.style="rainbow"){
    #Object: To display probability and effect over MAF
    #Input: MAF vector of MAF
    #Input: P vector of P values
    #Output: A table and plot
    #Requirment: NA
    #Authors: Zhiwu Zhang
    # Start  date: April 5, 2013
    # Last update: Oct 27, 2015 by Jiabo Wang add notice for P<0.1 is empty
    ##############################################################################################
    #print("MAF plot started")
    #print(threshold.output)
    #Remove NAs and under threshold
    index= which(P<threshold.output & !is.na(MAF))
    MAF=MAF[index]
    #E=E[index]
    P=P[index]
LP=-log10(P) 
LPC=round(LP*10,digits = 0)+20
ncolors=max(LPC,na.rm=T)
if(ncolors > 1024) {ncolors=1024}
if(ncolors==-Inf) 
{
print("There are no significant gene by this method(<0.1)")
}else{
#print("MAF plot started 0001")
#print(length(P))
#print(ncolors)
#palette(rainbow(ncolors))
#palette(gray(seq(.9,0,len = ncolors)))
#print("MAF plot started 0001b")
pdf(paste("GAPIT.", trait,".MAF.pdf" ,sep = ""), width = 5,height=5) 
par(mar = c(5,6,5,3))
theColor=heat.colors(ncolors, alpha = 1)
palette(rev(theColor))
plot(MAF,LP,type="p",lty = 1,lwd=2,col=LPC,xlab="MAF",ylab =expression(Probability~~-log[10](italic(p))),main = trait, cex.axis=1.1, cex.lab=1.3)
#for(i in 2:nc){
#lines(power[,i]~FDR, lwd=2,type="o",pch=i,col=i)
#}
#legend("bottomright", colnames(power), pch = c(1:nc), lty = c(1,2),col=c(1:nc))
palette("default")      # reset back to the default
dev.off()
}
}   #GAPIT.MAF ends here
#=============================================================================================

