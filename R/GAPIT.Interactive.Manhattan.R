`GAPIT.Interactive.Manhattan`<-
function(GWAS=NULL,MAF.threshold=seq(0,0.5,.1),cutOff=0.01,DPP=50000,X_fre=NULL,plot.type=c("m","q"),name.of.trait = "Trait"
  )
{   
    if(is.null(GWAS)) stop("Please add GWAS result in here!!!")
 
    MP=GWAS[,2:4]
    #print(head(GWAS))
    GWAS=GWAS[order(GWAS[,3]),]
    GWAS=GWAS[order(GWAS[,2]),]
    #print(GWAS[GWAS[,4]==min(GWAS[,4]),2])

    taxa=as.character(GWAS[,1])
    numMarker=nrow(GWAS)
    bonferroniCutOff01=-log10(0.01/numMarker)
    bonferroniCutOff05=-log10(0.05/numMarker)
    # deal with P value to log
    Ps=as.numeric(as.vector(GWAS[,4]))
    logPs <-  -log10(Ps)
    logPs[is.na(logPs)]=0
    
    y.lim <- ceiling(max(GWAS[,4]))
    chrom_total=as.numeric(as.character((GWAS[,2])))
    #print(GWAS[GWAS[,4]==min(GWAS[,4]),])
    #print("!!!!")
    #print(chrom_total[logPs==max(logPs)])
    POS=as.numeric(as.vector(GWAS[,3]))
    #print(head(POS))
    chm.to.analyze <- as.numeric(as.character(unique(GWAS[,2])))
    chm.to.analyze=chm.to.analyze[order(as.numeric(as.character(chm.to.analyze)))]
    #chm.to.analyze = factor(sort(chm.to.analyze))

    numCHR= length(chm.to.analyze)
    print(chm.to.analyze)

    ticks=NULL
    lastbase=0
    
        #change base position to accumulatives (ticks)
        for (i in chm.to.analyze)
        {
            index=(chrom_total==i)
            ticks <- c(ticks, lastbase+mean(POS[index]))
            POS[index]=POS[index]+lastbase
            lastbase=max(POS[index])
        }

        x0 <- POS
        y0 <- as.numeric(logPs)
        z0 <- chrom_total
        posi0<-as.numeric(as.vector(GWAS$Position))
        maf0 <- as.numeric(as.vector(GWAS$maf))
        effect0<- as.numeric(as.vector(GWAS$effect))
        #print(head(z0))
        position=order(y0,decreasing = TRUE)
        index0=GAPIT.Pruning(y0[position],DPP=DPP)
        index=position[index0]
        #order by P value
        x=x0[index]
        y=y0[index]
        z=z0[index]
        taxa=taxa[index]
        posi=posi0[index]
        maf=maf0[index]
        effect=effect0[index]
    
        plot.color=rep(c(  '#EC5f67',    '#FAC863',  '#99C794',    '#6699CC',  '#C594C5'),ceiling(numCHR/5))
    

if(c("m")%in%plot.type)
{

  Position=x
  P_value=y
  z[z<10]=paste("0",z[z<10],sep="")
  zz=paste("Chr_",z,sep="")
  #print(zz)
#  if(!require(plotly)) install.packages("plotly")
  #print("!!!!!")
  #print(head(Position))
 library(plotly)
  p <- plotly::plot_ly(
    type = 'scatter',
    x = ~Position,
    y = ~P_value,
    colorscale='Viridis',
    reversescale =T,
    #symbol="circle",
    text = ~paste("SNP: ", taxa, "<br>Posi: ", posi,"<br>MAF: ", round(maf,2),"<br>Effect: ",round(effect,2)),
    color = ~as.character(zz)
    )%>%
   plotly::add_trace(y=bonferroniCutOff01,name = 'CutOff-0.01',color=I("red"),mode="line",width=1.4,text="")%>%
   plotly::add_trace(y=bonferroniCutOff05,name = 'CutOff-0.05',color=I("red"),mode="line",line=list(width=1.4,dash='dot'),text="")%>%
   layout(title = "Interactive.Manhattan.Plot",
                  #showticklabels = FALSE,
                  #legend = list(orientation = 'h'),
                  xaxis = list(title = "Chromsome",zeroline = FALSE,showticklabels = FALSE),
                  yaxis = list (title = "-Log10(p)"))
   # plotly::add_trace(p,y=bonferroniCutOff01,name = 'CutOff-0.01',color=I("red"),mode="line",width=1.4,text="")%>%
   # plotly::add_trace(p,y=bonferroniCutOff05,name = 'CutOff-0.05',color=I("red"),mode="line",line=list(width=1.4,dash='dot'),text="")%>%
   # plotly::layout(title = "Interactive.Manhattan.Plot",
   #       #showticklabels = FALSE,
   #       #legend = list(orientation = 'h'),
   #       xaxis = list(title = "Chromsome",zeroline = FALSE,showticklabels = FALSE),
   #       yaxis = list (title = "-Log10(p)"))

    htmltools::save_html(p, paste("Interactive.Manhattan.",name.of.trait,".html",sep=""))
}


################ for QQ plot
if(c("q")%in%plot.type)
{
        P.values=y
        p_value_quantiles <- (1:length(P.values))/(length(P.values)+1)
        log.P.values <- P.values
        log.Quantiles <- -log10(p_value_quantiles)
        
        index=GAPIT.Pruning(log.P.values,DPP=DPP)
        log.P.values=log.P.values[index ]
        log.Quantiles=log.Quantiles[index]
        N=length(P.values)
        N1=length(log.Quantiles)
        ## create the confidence intervals
        c95 <- rep(NA,N1)
        c05 <- rep(NA,N1)
        for(j in 1:N1){
            i=ceiling((10^-log.Quantiles[j])*N)
            if(i==0)i=1
            c95[j] <- stats::qbeta(0.95,i,N-i+1)
            c05[j] <- stats::qbeta(0.05,i,N-i+1)
            #print(c(j,i,c95[j],c05[j]))
        }
        
        #CI shade
        #plot3d(NULL, xlim = c(0,max(log.Quantiles)), zlim = c(0,max(log.P.values)), type="l",lty=5, lwd = 2, axes=FALSE, xlab="", ylab="",col="gray")
        index=length(c95):1
        zz=paste("Chr_",z,sep="")
        Expected=log.Quantiles
        Observed=log.P.values
        #abline(a = 0, b = 1, col = "red",lwd=2)
        qp <- plotly::plot_ly(
    type = 'scatter',
    x = ~Expected,
    y = ~Observed,
    text = ~paste("SNP: ", taxa,"<br>Chr: ",zz,"<br>Posi: ", posi, "<br>MAF: ", round(maf,2),"<br>Effect: ",round(effect,2)),
    #size=2*y/max(y),
    name = "SNP",
    opacity=0.5,
    ) %>% plotly::add_lines(x=log.Quantiles,y=log.Quantiles,color=I("red"), 
    mode = 'lines',name="Diag",text="")%>%
          plotly::layout(title = "Interactive.QQ.Plot",
                         xaxis = list(title = "Expected -Log10(p)"),
                         yaxis = list (title = "Observed -Log10(p)"),
                         #showticklabels = FALSE,
                         showlegend = FALSE)
    htmltools::save_html(qp, paste("Interactive.QQ ",name.of.trait,".html",sep=""))

    # plotly::layout(title = "Interactive.QQ.Plot",
    #     xaxis = list(title = "Expected -Log10(p)"),
    #      yaxis = list (title = "Observed -Log10(p)"),
    #      #showticklabels = FALSE,
    #      showlegend = FALSE)
    #     htmltools::save_html(qp, paste("Interactive.QQ ",name.of.trait,".html",sep=""))


}   
print("GAPIT.Interactive.Plot has done !!!")

}#end of GAPIT.Interactive.Manhattan
#=============================================================================================





