`GAPIT.Manhattan` <-
function(GI.MP = NULL,GD=NULL,name.of.trait = "Trait",plot.type = "Genomewise",width0=18,height0=5.75,
DPP=50000,cutOff=0.01,band=5,seqQTN=NULL,plot.style="Oceanic",CG=NULL,plot.bin=10^9,chor_taxa=NULL){
    #Object: Make a Manhattan Plot
    #Options for plot.type = "Separate_Graph_for_Each_Chromosome" and "Same_Graph_for_Each_Chromosome"
    #Output: A pdf of the Manhattan Plot
    #Authors: Alex Lipka, Zhiwu Zhang, Meng Li and Jiabo Wang
    # Last update: Oct 10, 2016
    #Add r2 between candidata SNP and other markers in on choromosome
    ##############################################################################################
    #print("Manhattan ploting...")
    
    #print(cutOff)
    #do nothing if null input
    if(is.null(GI.MP)) return
  #Handler of lable paosition only indicated by negatie position
    position.only=F
    if(!is.null(seqQTN)){
      if(seqQTN[1]<0){
        seqQTN=-seqQTN
        position.only=T
      }      
    }  
    borrowSlot=4
    GI.MP[,borrowSlot]=0 
    GI.MP[,5]=1:(nrow(GI.MP))
    GI.MP=matrix(as.numeric(as.matrix(GI.MP) ) ,nrow(GI.MP),ncol(GI.MP))
    GI.MP=GI.MP[order(GI.MP[,2]),]
    GI.MP=GI.MP[order(GI.MP[,1]),]
    #Inicial as 0   
    if(!is.null(seqQTN))GI.MP[seqQTN,borrowSlot]=1   
    if(!is.null(GD))
    {  if(ncol(GD)!=nrow(GI.MP))print("GD does not mach GM in Manhattan !!!")
    }
    #Remove all SNPs that do not have a choromosome, bp position and p value(NA)
    GI.MP <- GI.MP[!is.na(GI.MP[,1]),]
    GI.MP <- GI.MP[!is.na(GI.MP[,2]),]
    if(!is.null(GD)) GD=GD[,!is.na(GI.MP[,3])]
    GI.MP <- GI.MP[!is.na(GI.MP[,3]),]
    
    #Retain SNPs that have P values between 0 and 1 (not na etc)
    if(!is.null(GD)) GD=GD[,GI.MP[,3]>0]
    GI.MP <- GI.MP[GI.MP[,3]>0,]
    if(!is.null(GD)) GD=GD[,GI.MP[,3]<=1]
    GI.MP <- GI.MP[GI.MP[,3]<=1,]
    
    #Remove chr 0 and 99
    GI.MP <- GI.MP[GI.MP[,1]!=0,]
    numMarker=nrow(GI.MP)
    #print(numMarker)
    bonferroniCutOff=-log10(cutOff/numMarker)
    sp=sort(GI.MP[,3])
    spd=abs(cutOff-sp*numMarker/cutOff)
    index_fdr=grep(min(spd),spd)[1]
    FDRcutoff=-log10(cutOff*index_fdr/numMarker)
    #Replace P the -log10 of the P-values
    if(!is.null(GD))
    {  if(ncol(GD)!=nrow(GI.MP))
    {print("GD does not match GM in Manhattan !!!")
    return
    }}
    #print(ncol(GD))
    #print(nrow(GI.MP))
    GI.MP[,3] <-  -log10(GI.MP[,3])
    index_GI=GI.MP[,3]>0
    GI.MP <- GI.MP[index_GI,]
    if(!is.null(GD)) GD=GD[,index_GI]
    
    GI.MP[,5]=1:(nrow(GI.MP))
    y.lim <- ceiling(max(GI.MP[,3]))
    chm.to.analyze <- unique(GI.MP[,1])
    chm.to.analyze=chm.to.analyze[order(chm.to.analyze)]
    numCHR= length(chm.to.analyze)
    #GI.MP[,5]=1:(nrow(GI.MP))
     bin.mp=GI.MP[,1:3]
     bin.mp[,3]=0 # for r2
     bin.mp[,1]=as.numeric(as.vector(GI.MP[,2]))+as.numeric(as.vector(GI.MP[,1]))*(10^(max(GI.MP[,1])+1))
     
     
     #as.numeric(as.vector(GP[,3]))+as.numeric(as.vector(GP[,2]))*MaxBP
     #print(head(bin.mp))
     bin.mp[,2]=floor(bin.mp[,1]/plot.bin)
     if(!is.null(GD)) X=GD

     #print(head(bin.mp))
        #Chromosomewise plot
    if(plot.type == "Chromosomewise"&!is.null(GD))
    {
        #print("Manhattan ploting Chromosomewise")
        GI.MP=cbind(GI.MP,bin.mp)
        pdf(paste("GAPIT.", name.of.trait,".Manhattan.Plot.Chromosomewise.pdf" ,sep = ""), width = 10)
            #par(mar = c(5,5,4,3), lab = c(8,5,7))
        layout(matrix(c(1,1,2,1,1,1,1,1,1),3,3,byrow=TRUE), c(2,1), c(1,1), TRUE)
        for(i in 1:numCHR)
        {
            #Extract SBP on this chromosome
            subset=GI.MP[GI.MP[,1]==chm.to.analyze[i],,drop=FALSE]
            # print(head(subset))
            if(nrow(subset)==0)next #thanks to lvclark to fix it
            subset[,1]=1:(nrow(subset))
            #sub.bin.mp=bin.mp[GI.MP[,1]==chm.to.analyze[i],]
            #subset=cbind(subset,sub.bin.mp)
            sig.mp=subset[subset[,3]>bonferroniCutOff,,drop=FALSE]
            sig.index=subset[,3]>bonferroniCutOff ### index of significont SNP          
            num.row=nrow(sig.mp)
            if(!is.null(dim(sig.mp)))sig.mp=sig.mp[!duplicated(sig.mp[,7]),]
            num.row=nrow(sig.mp)
            if(is.null(dim(sig.mp))) num.row=1
            bin.set=NULL
            r2_color=matrix(0,nrow(subset),2)
            #r2_color
            print(paste("select ",num.row," candidate significont markers in ",i," chromosome ",sep="") )
            #print(sig.mp)
            if(length(unique(sig.index))==2)
            {
                for(j in 1:num.row)
                {   sig.mp=matrix(sig.mp,num.row,8)
                    bin.store=subset[which(subset[,7]==sig.mp[j,7]),]
                    if(is.null(dim(bin.store)))
                      {subset[which(subset[,7]==sig.mp[j,7]),8]=1
                          next
                      }
                    bin.index=unique(bin.store[,5])
                    subGD=X[,bin.store[,5]]
                    #print(dim(bin.store))
                    if(is.null(CG))candidata=bin.store[bin.store[,3]==max(bin.store[,3]),5]
                    if(length(candidata)!=1)candidata=candidata[1]
                    
                    for (k in 1:ncol(subGD))
                    {
                        r2=cor(X[,candidata],subGD[,k])^2
                        #print(r2)
                        bin.store[k,8]=r2
                        
                    }
                    subset[bin.store[,1],8]=bin.store[,8]
                    #print()
                }###end for each sig.mp
         
            }###end if empty of sig.mp
            rm(sig.mp,num.row)
            y.lim <- ceiling(max(subset[,3]))+1  #set upper for each chr
            if(length(subset)>3){
                x <- as.numeric(subset[,2])/10^(6)
                y <- as.numeric(subset[,3])
            }else{
                x <- as.numeric(subset[2])/10^(6)
                y <- as.numeric(subset[3])
            }
            
            ##print(paste("befor prune: chr: ",i, "length: ",length(x),"max p",max(y), "min p",min(y), "max x",max(x), "Min x",min(x)))
            n_col=10
            r2_color[,2]=subset[,8]
            do_color=colorRampPalette(c("orangeRed", "blue"))(n_col)
            #Prune most non important SNPs off the plots
            order=order(y,decreasing = TRUE)
            y=y[order]
            x=x[order]
            r2_color=r2_color[order,,drop=FALSE]
            index=GAPIT.Pruning(y,DPP=round(DPP/numCHR))
            x=x[index]
            y=y[index]
            r2_color=r2_color[index,,drop=FALSE]
            r2_color[which(r2_color[,2]<=0.2),2]=do_color[n_col]
            r2_color[which(r2_color[,2]<=0.4&r2_color[,2]>0.2),2]=do_color[n_col*0.8]
            r2_color[which(r2_color[,2]<=0.6&r2_color[,2]>0.4),2]=do_color[n_col*0.6]
            r2_color[which(r2_color[,2]<=0.8&r2_color[,2]>0.6),2]=do_color[n_col*0.4]
            r2_color[which(r2_color[,2]<=1&r2_color[,2]>0.8),2]=do_color[n_col/n_col]
            par(mar=c(0,0,0,0))
            par(mar=c(5,5,2,1),cex=0.8)

            plot(y~x,type="p", ylim=c(0,y.lim), xlim = c(min(x), max(x)),las=1,
            col = r2_color[,2], xlab = expression(Base~Pairs~(x10^-6)),
            ylab = "-Log Base 10 p-value", main =           paste("Chromosome",chm.to.analyze[i],sep=" "),
            cex.lab=1.6,pch=21,bg=r2_color[,2])
            
            abline(h=bonferroniCutOff,col="forestgreen")
            abline(h=FDRcutoff,col="forestgreen",lty=2)
            par(mar=c(15,5,6,5),cex=0.5)
            
            barplot(matrix(rep(1,times=n_col),n_col,1),beside=T,col=do_color,border=do_color,axes=FALSE,)
        #legend(x=10,y=2,legend=expression(R^"2"),,lty=0,cex=1.3,bty="n",bg=par("bg"))
            axis(3,seq(11,1,by=-2),seq(0,1,by=0.2),las=1)

        }# end plot.type == "Chromosomewise"&!is.null(GD)
        dev.off()
        
        print("manhattan plot on chromosome finished")
    } #Chromosomewise plot
    
    
    #Genomewise plot
    if(plot.type == "Genomewise")
    {
        #print("Manhattan ploting Genomewise")
        #Set corlos for chromosomes
        #nchr=max(chm.to.analyze)
        nchr=length(chm.to.analyze)

    #Set color schem            
        ncycle=ceiling(nchr/band)
        ncolor=band*ncycle
        #palette(rainbow(ncolor+1))
        cycle1=seq(1,nchr,by= ncycle)
        thecolor=cycle1
        for(i in 2:ncycle){thecolor=c(thecolor,cycle1+(i-1))}
        col.Rainbow=rainbow(ncolor+1)[thecolor]         
          col.FarmCPU=rep(c("#CC6600","deepskyblue","orange","forestgreen","indianred3"),ceiling(numCHR/5))
          col.Rushville=rep(c("orangered","navyblue"),ceiling(numCHR/2))    
            col.Congress=rep(c("deepskyblue3","firebrick"),ceiling(numCHR/2))
            col.Ocean=rep(c("steelblue4","cyan3"),ceiling(numCHR/2))        
            col.PLINK=rep(c("gray10","gray70"),ceiling(numCHR/2))       
            col.Beach=rep(c("turquoise4","indianred3","darkolivegreen3","red","aquamarine3","darkgoldenrod"),ceiling(numCHR/5))
            #col.Oceanic=rep(c( '#EC5f67',  '#F99157',  '#FAC863',  '#99C794',  '#5FB3B3',  '#6699CC',  '#C594C5',  '#AB7967'),ceiling(numCHR/8))
            #col.Oceanic=rep(c( '#EC5f67',      '#FAC863',  '#99C794',      '#6699CC',  '#C594C5',  '#AB7967'),ceiling(numCHR/6))
            col.Oceanic=rep(c(  '#EC5f67',      '#FAC863',  '#99C794',      '#6699CC',  '#C594C5'),ceiling(numCHR/5))
            col.cougars=rep(c(  '#990000',      'dimgray'),ceiling(numCHR/2))
        
        if(plot.style=="Rainbow")plot.color= col.Rainbow
        if(plot.style =="FarmCPU")plot.color= col.Rainbow
        if(plot.style =="Rushville")plot.color= col.Rushville
        if(plot.style =="Congress")plot.color= col.Congress
        if(plot.style =="Ocean")plot.color= col.Ocean
        if(plot.style =="PLINK")plot.color= col.PLINK
            if(plot.style =="Beach")plot.color= col.Beach
            if(plot.style =="Oceanic")plot.color= col.Oceanic
            if(plot.style =="cougars")plot.color= col.cougars
        
        #FarmCPU uses filled dots
        mypch=1
        if(plot.style =="FarmCPU")mypch=20
                
        GI.MP <- GI.MP[order(GI.MP[,2]),]
        GI.MP <- GI.MP[order(GI.MP[,1]),]

        ticks=NULL
        lastbase=0
        
        #print("Manhattan data sorted")
        #print(chm.to.analyze)
        
        #change base position to accumulatives (ticks)
        for (i in chm.to.analyze)
        {
            index=(GI.MP[,1]==i)
            ticks <- c(ticks, lastbase+mean(GI.MP[index,2]))
            GI.MP[index,2]=GI.MP[index,2]+lastbase
            lastbase=max(GI.MP[index,2])
        }
        
        #print("Manhattan chr processed")
        #print(length(index))
        #print(length(ticks))
        #print((ticks))
        #print((lastbase))
        
        x0 <- as.numeric(GI.MP[,2])
        y0 <- as.numeric(GI.MP[,3])
        z0 <- as.numeric(GI.MP[,1])
        position=order(y0,decreasing = TRUE)
        index0=GAPIT.Pruning(y0[position],DPP=DPP)
        index=position[index0]
        
        x=x0[index]
        y=y0[index]
        z=z0[index]

        #Extract QTN
        QTN=GI.MP[which(GI.MP[,borrowSlot]==1),]
        #print(QTN)
        #Draw circles with same size and different thikness
        size=1 #1
        ratio=10 #5
        base=1 #1
        themax=ceiling(max(y))
        themin=floor(min(y))
        wd=((y-themin+base)/(themax-themin+base))*size*ratio
        s=size-wd/ratio/2
        
        #print("Manhattan XY created")
       ####xiaolei update on 2016/01/09 
        if(plot.style =="FarmCPU"){
        pdf(paste("FarmCPU.", name.of.trait,".Manhattan.Plot.Genomewise.pdf" ,sep = ""), width = width0,height=height0)
        }else{
        pdf(paste("GAPIT.", name.of.trait,".Manhattan.Plot.Genomewise.pdf" ,sep = ""), width = width0,height=height0)
        }
            par(mar = c(3,6,5,1))
            plot(y~x,xlab="",ylab=expression(-log[10](italic(p))) ,las=1,
            cex.axis=1, cex.lab=1.3 ,col=plot.color[z],axes=FALSE,type = "p",pch=mypch,lwd=wd,cex=s+.3,main = paste(name.of.trait,sep="             "),cex.main=2.5)
        
        #Label QTN positions
        if(is.vector(QTN)){
          if(position.only){abline(v=QTN[2], lty = 2, lwd=1.5, col = "grey")}else{
          points(QTN[2], QTN[3], type="p",pch=21, cex=2,lwd=1.5,col="dimgrey")
          points(QTN[2], QTN[3], type="p",pch=20, cex=1,lwd=1.5,col="dimgrey")
          }
        }else{
          if(position.only){abline(v=QTN[,2], lty = 2, lwd=1.5, col = "grey")}else{
          points(QTN[,2], QTN[,3], type="p",pch=21, cex=2,lwd=1.5,col="dimgrey")
          points(QTN[,2], QTN[,3], type="p",pch=20, cex=1,lwd=1.5,col="dimgrey")
          }
        }
        
        #Add a horizontal line for bonferroniCutOff
        abline(h=bonferroniCutOff,col="forestgreen")
        #Add FDR line
        abline(h=FDRcutoff,col="forestgreen",lty=2)
        #print(bonferroniCutOff)
        #Set axises
        # jiabo creat chor_taxa
        #print(chor_taxa)
        if(length(chor_taxa)!=length(ticks))chor_taxa=NULL
        #print(unique(GI.MP[,1]))
        if(!is.null(chor_taxa))
        {axis(1, at=ticks,cex.axis=1,labels=chor_taxa,tick=T,gap.axis=0.25)
        }else{axis(1, at=ticks,cex.axis=1,labels=chm.to.analyze,tick=F)}
        axis(2, at=1:themax,cex.axis=1,las=1,labels=1:themax,gap.axis=3,tick=F)

        box()
        palette("default")
        dev.off()
        #print("Manhattan done Genomewise")
        
    } #Genomewise plot
    
    print("GAPIT.Manhattan accomplished successfully!zw")
} #end of GAPIT.Manhattan
#=============================================================================================
