`GAPIT.Multiple.Manhattan` <-
function(model_store,DPP=50000,chor_taxa=NULL,cutOff=0.01,band=5,seqQTN=NULL,Y=NULL,GM=NULL,interQTN=NULL,
    plot.style="Oceanic",plot.line=TRUE,allpch=NULL,plot.type=c("h","s","w")){
    #Object: Make a Manhattan Plot
    #Output: pdfs of the Multiple Manhattan Plot
    #Authors: Zhiwu Zhang and Jiabo Wang
    # Last update: Feb 22, 2022
    ##############################################################################################
  Nenviron=length(model_store)*(ncol(Y)-1)
  environ_name=NULL
  new_xz=NULL
  for(i in 1:length(model_store))
  {
    for(j in 1:(ncol(Y)-1))
    {
      environ_name=c(environ_name,paste(model_store[i],".",colnames(Y)[-1][j],sep=""))
    }
  }
sig_pos=NULL
simulation=FALSE
    if(!is.null(seqQTN)){    
        #seqQTN=-seqQTN
        simulation=TRUE    
    }
themax.y0=NULL
for(i in 1:length(environ_name))
{
  print(paste("Reading GWAS result with ",environ_name[i],sep=""))
  environ_result=read.csv(paste("GAPIT.",environ_name[i],".GWAS.Results.csv",sep=""),head=T)
  environ_result=environ_result[order(environ_result[,3]),]
  environ_result=environ_result[order(environ_result[,2]),]
  environ_filter=environ_result[!is.na(environ_result[,4]),]
  themax.y=round(max(-log10(environ_filter[,4])),0)
  themax.y0=round(max(c(themax.y,themax.y0)),0)
  y_filter=environ_filter[environ_filter[,4]<(cutOff/(nrow(environ_filter))),]
  write.table(y_filter,paste("GAPIT.Filter_",environ_name[i],"_GWAS_result.txt",sep=""))

  result=environ_result[,1:4]
  result=result[match(as.character(GM[,1]),as.character(result[,1])),]
  rownames(result)=1:nrow(result)
  #print(i)
  if(i==1){
    result0=result
    colnames(result0)[4]=environ_name[i]
    }
  if(i!=1){
    result0=merge(result0,result[,c(1,4)],by.x=colnames(result0)[1],by.y=colnames(result)[1])
    colnames(result0)[i+3]=environ_name[i]
    }
  rownames(result)=1:nrow(result)
  result[is.na(result[,4]),4]=1
  sig_pos=append(sig_pos,as.numeric(rownames(result[result[!is.na(result[,4]),4]<(cutOff/nrow(result)),])))
}
#if(length(sig_pos)!=0)sig_pos=sig_pos[!duplicated(sig_pos)]
 if(length(sig_pos[!is.na(sig_pos)])!=0)
 {     x_matrix=as.matrix(table(sig_pos))
       x_matrix=cbind(as.data.frame(rownames(x_matrix)),x_matrix)
       #print(x_matrix)
       lastbase=0
       map_store=as.matrix(cbind(as.character(GM[,2]),as.numeric(as.vector(GM[,3]))))
       #print(head(map_store))
       #print(as.numeric(map_store[,3]))
        for (j in unique(map_store[,1]))
        {
            index=map_store[,1]==j
            # print(as.numeric(map_store[index,2]))
            map_store[index,2]=as.numeric(map_store[index,2])+lastbase
            lastbase=max(as.numeric(map_store[index,2]))
            #print(lastbase)
        }
       colnames(x_matrix)=c("pos","times")
       new_xz=cbind(x_matrix,map_store[as.numeric(as.character(x_matrix[,1])),,drop=FALSE])
       colnames(new_xz)=c("pos","times","chro","xlab")
       new_xz=new_xz[!duplicated(new_xz),]
       new_xz[new_xz[,2]>=3,2]=3
       new_xz[,2]=4-new_xz[,2]
       new_xz[new_xz[,2]==3,2]=0

       new_xz=as.matrix(new_xz)
       new_xz=new_xz[new_xz[,2]!="0",]
       new_xz=matrix(new_xz,length(as.vector(new_xz))/4,4)
}
# setup colors
chm.to.analyze <- unique(result[,1])
nchr=length(chm.to.analyze)
size=1 #1
ratio=10 #5
base=1 #1
numCHR=nchr
ncycle=ceiling(nchr/5)
ncolor=band*ncycle
thecolor=seq(1,nchr,by= ncycle)
col.Rainbow=rainbow(ncolor+1)     
col.FarmCPU=rep(c("#CC6600","deepskyblue","orange","forestgreen","indianred3"),ceiling(numCHR/5))
col.Rushville=rep(c("orangered","navyblue"),ceiling(numCHR/2))    
col.Congress=rep(c("deepskyblue3","firebrick"),ceiling(numCHR/2))
col.Ocean=rep(c("steelblue4","cyan3"),ceiling(numCHR/2))    
col.PLINK=rep(c("gray10","gray70"),ceiling(numCHR/2))     
col.Beach=rep(c("turquoise4","indianred3","darkolivegreen3","red","aquamarine3","darkgoldenrod"),ceiling(numCHR/5))
col.Oceanic=rep(c(  '#EC5f67',    '#FAC863',  '#99C794',    '#6699CC',  '#C594C5'),ceiling(numCHR/5))
col.cougars=rep(c(  '#990000',    'dimgray'),ceiling(numCHR/2))  
if(plot.style=="Rainbow")plot.color= col.Rainbow
if(plot.style =="FarmCPU")plot.color= col.Rainbow
if(plot.style =="Rushville")plot.color= col.Rushville
if(plot.style =="Congress")plot.color= col.Congress
if(plot.style =="Ocean")plot.color= col.Ocean
if(plot.style =="PLINK")plot.color= col.PLINK
if(plot.style =="Beach")plot.color= col.Beach
if(plot.style =="Oceanic")plot.color= col.Oceanic
if(plot.style =="cougars")plot.color= col.cougars  

if("h"%in%plot.type)
{
    pdf(paste("GAPIT.Manhattan.Mutiple.Plot.high",".pdf" ,sep = ""), width = 20,height=6*Nenviron)
    par(mfrow=c(Nenviron,1))
    mypch=1
    for(k in 1:Nenviron)
    { 
       if(k==Nenviron)
        {
        par(mar = c(3.5,8,0,8))
         # par(pin=c(10,((8-mtext.h)/Nenviron)+mtext.h))

        }else{
            #par(mfrow=c(Nenviron,1))
        par(mar = c(1.5,8,0.5,8))    
        }
       environ_result=read.csv(paste("GAPIT.",environ_name[k],".GWAS.Results.csv",sep=""),head=T)
       result=environ_result[,1:4]
       result=result[match(as.character(GM[,1]),as.character(result[,1])),]
       rownames(result)=1:nrow(result)
       GI.MP=result[,c(2:4)]
       borrowSlot=4
       GI.MP[,borrowSlot]=0 #Inicial as 0
       GI.MP[,5]=1:(nrow(GI.MP))
       GI.MP[,6]=1:(nrow(GI.MP)) 
       GI.MP <- GI.MP[!is.na(GI.MP[,1]),]
       GI.MP <- GI.MP[!is.na(GI.MP[,2]),]
       GI.MP[is.na(GI.MP[,3]),3]=1
    #Retain SNPs that have P values between 0 and 1 (not na etc)
       GI.MP <- GI.MP[GI.MP[,3]>0,]
       GI.MP <- GI.MP[GI.MP[,3]<=1,]
    #Remove chr 0 and 99
       GI.MP <- GI.MP[GI.MP[,1]!=0,]
       total_chromo=length(unique(GI.MP[,1]))
    # print(dim(GI.MP))
       if(!is.null(seqQTN))GI.MP[seqQTN,borrowSlot]=1
       numMarker=nrow(GI.MP)
       bonferroniCutOff=-log10(cutOff/numMarker)
       GI.MP[,3] <-  -log10(GI.MP[,3])
       GI.MP[,5]=1:numMarker
       y.lim <- ceiling(max(GI.MP[,3]))  
       chm.to.analyze <- unique(GI.MP[,1])
       nchr=length(chm.to.analyze)
       GI.MP[,6]=1:(nrow(GI.MP))
       MP_store=GI.MP
       index_GI=MP_store[,3]>=0
       MP_store <- MP_store[index_GI,]
       ticks=NULL
       lastbase=0
       for(i in chm.to.analyze)
          {
           index=(MP_store[,1]==i)
           ticks <- c(ticks, lastbase+mean(MP_store[index,2]))
           MP_store[index,2]=MP_store[index,2]+lastbase
           lastbase=max(MP_store[index,2])
          }
       x0 <- as.numeric(MP_store[,2])
       y0 <- as.numeric(MP_store[,3])
       z0 <- as.character(MP_store[,1])
       # convert chromosome character to number
       chor_taxa=as.character(unique(MP_store[,1]))
       chor_taxa=chor_taxa[order(as.numeric(as.character(chor_taxa)))]
       chr_letter=grep("[A-Z]|[a-z]",chor_taxa)
       if(!setequal(integer(0),chr_letter))
         {     
           z0=as.character(MP_store[,1])
           for(i in 1:(length(chor_taxa)))
              {
                index=z0==chor_taxa[i]
                z0[index]=i    
              }
          }
       z0=as.numeric(z0)
       max.x=NULL
       for(i in chm.to.analyze)
          {
           index=(MP_store[,1]==i)
           max.x=c(max.x,max(x0[index]))
          }
       max.x=c(min(x0),max.x)
       x1=sort(x0)
       position=order(y0,decreasing = TRUE)
       values=y0[position]
       if(length(values)<=DPP)
         {
         index=position[c(1:length(values))]
         }else{       
         values=sqrt(values)  #This shift the weight a little bit to the low building.
        #Handler of bias plot
         rv=runif(length(values))
         values=values+rv
         values=values[order(values,decreasing = T)]         
         theMin=min(values)
         theMax=max(values)
         range=theMax-theMin
         interval=range/DPP
         ladder=round(values/interval)
         ladder2=c(ladder[-1],0)
         keep=ladder-ladder2
         index=position[which(keep>=0)]
         }        
       x=x0[index]
       y=y0[index]
       z=z0[index]
        #Extract QTN
       QTN=MP_store[which(MP_store[,borrowSlot]==1),]
        #Draw circles with same size and different thikness
       
       themax=ceiling(max(y))
       themax2=ceiling((ceiling(themax/4))*4)

       themin=floor(min(y))
       wd=((y-themin+base)/(themax-themin+base))*size*ratio
       s=size-wd/ratio/2
       plot(y~x,xlab="",ylab="" ,ylim=c(0,themax2),xlim=c(min(x),max(x)),
           cex.axis=4, cex.lab=4, ,col=plot.color[z],axes=FALSE,type = "p",
           pch=mypch,lwd=wd,cex=s+2.5,cex.main=4)
       mtext(side=2,expression(-log[10](italic(p))),line=3.5, cex=2.5)
       if(!simulation)
         {
          abline(v=QTN[2], lty = 2, lwd=1.5, col = "grey")}else{
          points(QTN[,2], QTN[,3], pch=20, cex=2.5,lwd=2.5,col="dimgrey")
         }        
       if(plot.line)
         {
          if(!is.null(nrow(new_xz)))  
            {
             abline(v=as.numeric(new_xz[,4]),col="grey",lty=as.numeric(new_xz[,2]),untf=T,lwd=3)
            }else{
             abline(v=as.numeric(new_xz[1]),col=plot.color[as.numeric(new_xz[3])],lty=as.numeric(new_xz[2]),untf=T,lwd=3)
            }
         }
        #Add a horizontal line for bonferroniCutOff
       abline(h=bonferroniCutOff,lty=1,untf=T,lwd=3,col="forestgreen")
       axis(2, yaxp=c(0,themax2,4),cex.axis=2.3,tick=T,las=1,lwd=2.5)
       if(k==Nenviron)axis(1, at=max.x,cex.axis=2.5,labels=rep("",length(max.x)),tick=T,lwd=2.5)
       if(k==Nenviron)axis(1, at=ticks,cex.axis=2.5,labels=chm.to.analyze,tick=F,line=1)
       mtext(side=4,paste(environ_name[k],sep=""),line=3.2,cex=2)
    }#end of environ_name
       dev.off()
}#end of plot.type

if("w"%in%plot.type)
{
 pdf(paste("GAPIT.Manhattan.Mutiple.Plot.wide",".pdf" ,sep = ""), width = 16,height=8.5)
 par(mfrow=c(Nenviron,1))
 mtext.h=0.5
 for(k in 1:Nenviron)
 { 
  if(k==Nenviron)
        {#par(mfrow=c(Nenviron,1))
          # print(par())
        par(mar = c(2,8,0,8))
         # par(pin=c(10,((8-mtext.h)/Nenviron)+mtext.h))

        }else{
            #par(mfrow=c(Nenviron,1))
        par(mar = c(1.5,8,0.5,8))    
         # par(pin=c(10,(8-mtext.h)/Nenviron))

        }
  environ_result=read.csv(paste("GAPIT.",environ_name[k],".GWAS.Results.csv",sep=""),head=T)
  #print(environ_result[as.numeric(new_xz[,1]),])
  result=environ_result[,1:4]
    result=result[match(as.character(GM[,1]),as.character(result[,1])),]
    rownames(result)=1:nrow(result)
    GI.MP=result[,c(2:4)]
    borrowSlot=4
    GI.MP[,borrowSlot]=0 #Inicial as 0
    GI.MP[,5]=1:(nrow(GI.MP))
    GI.MP[,6]=1:(nrow(GI.MP))
    
    
    GI.MP <- GI.MP[!is.na(GI.MP[,1]),]
    GI.MP <- GI.MP[!is.na(GI.MP[,2]),]
    GI.MP[is.na(GI.MP[,3]),3]=1
    
    #Retain SNPs that have P values between 0 and 1 (not na etc)
    GI.MP <- GI.MP[GI.MP[,3]>0,]
    GI.MP <- GI.MP[GI.MP[,3]<=1,]
    #Remove chr 0 and 99
    GI.MP <- GI.MP[GI.MP[,1]!=0,]
    total_chromo=length(unique(GI.MP[,1]))
    # print(dim(GI.MP))
    if(!is.null(seqQTN))GI.MP[seqQTN,borrowSlot]=1
    numMarker=nrow(GI.MP)
    bonferroniCutOff=-log10(cutOff/numMarker)
    GI.MP[,3] <-  -log10(GI.MP[,3])
    GI.MP[,5]=1:numMarker
    y.lim <- ceiling(max(GI.MP[,3]))
    
    chm.to.analyze <- unique(GI.MP[,1])
    # chm.to.analyze=chm.to.analyze[order(chm.to.analyze)]
    nchr=length(chm.to.analyze)
    GI.MP[,6]=1:(nrow(GI.MP))
    MP_store=GI.MP
        index_GI=MP_store[,3]>=0
        MP_store <- MP_store[index_GI,]
        ticks=NULL
        lastbase=0
        for (i in chm.to.analyze)
        {
            index=(MP_store[,1]==i)
            ticks <- c(ticks, lastbase+mean(MP_store[index,2]))
            MP_store[index,2]=MP_store[index,2]+lastbase
            lastbase=max(MP_store[index,2])
        }
        
        x0 <- as.numeric(MP_store[,2])
        y0 <- as.numeric(MP_store[,3])
        z0 <- as.character(MP_store[,1])
       # convert chromosome character to number
       chor_taxa=as.character(unique(MP_store[,1]))
       chor_taxa=chor_taxa[order(as.numeric(as.character(chor_taxa)))]
       chr_letter=grep("[A-Z]|[a-z]",chor_taxa)
       if(!setequal(integer(0),chr_letter))
         {     
           z0=as.character(MP_store[,1])
           for(i in 1:(length(chor_taxa)))
              {
                index=z0==chor_taxa[i]
                z0[index]=i    
              }
          }
       z0=as.numeric(z0)
        x1=sort(x0)

        position=order(y0,decreasing = TRUE)
        values=y0[position]
        if(length(values)<=DPP)
        {
         index=position[c(1:length(values))]
            }else{
         
        values=sqrt(values)  #This shift the weight a little bit to the low building.
        #Handler of bias plot
        rv=runif(length(values))
        values=values+rv
        values=values[order(values,decreasing = T)]

        theMin=min(values)
        theMax=max(values)
        range=theMax-theMin
        interval=range/DPP

        ladder=round(values/interval)
        ladder2=c(ladder[-1],0)
        keep=ladder-ladder2
        index=position[which(keep>=0)]
        }
        
        
        x=x0[index]
        y=y0[index]
        z=z0[index]
        # print(length(x))

        #Extract QTN
        #if(!is.null(seqQTN))MP_store[seqQTN,borrowSlot]=1
        #if(!is.null(interQTN))MP_store[interQTN,borrowSlot]=2
        QTN=MP_store[which(MP_store[,borrowSlot]==1),]
        #Draw circles with same size and different thikness
        themax=ceiling(max(y))
        themax2=ceiling((ceiling(themax/4))*4)
        themin=floor(min(y))
        # ratio=5
        wd=((y-themin+base)/(themax-themin+base))*size*ratio
        s=size-wd/ratio/2
        mypch=1
        bamboo=4
        # plot(y~x,xlab="",ylab="" ,ylim=c(0,themax),
        #     cex.axis=4, cex.lab=4, ,col=plot.color[z],axes=FALSE,type = "p",pch=mypch,lwd=0.5,cex=0.7,cex.main=2)
        plot(y~x,xlab="",ylab="" ,ylim=c(0,themax2),xlim=c(min(x),max(x)),
           cex.axis=4, cex.lab=4, ,col=plot.color[z],axes=FALSE,type = "p",
           pch=mypch,lwd=wd,cex=s+2.5,cex.main=4)
        mtext(side=2,expression(-log[10](italic(p))),line=3, cex=1)
        if(plot.line)
        {
          if(!is.null(nrow(new_xz)))  {abline(v=as.numeric(new_xz[,4]),col="grey",lty=as.numeric(new_xz[,2]),untf=T,lwd=2)
             }else{abline(v=as.numeric(new_xz[1]),col=plot.color[as.numeric(new_xz[3])],lty=as.numeric(new_xz[2]),untf=T,lwd=2)
             }
        }
        if(!simulation){abline(v=QTN[2], lty = 2, lwd=1.5, col = "grey")}else{
          points(QTN[,2], QTN[,3], type="p",pch=21, cex=1.5,lwd=1.5,col="dimgrey")
          points(QTN[,2], QTN[,3], type="p",pch=20, cex=1.5,lwd=1.5,col="dimgrey")
          }
        #Add a horizontal line for bonferroniCutOff
        abline(h=bonferroniCutOff,lty=1,untf=T,lwd=1,col="forestgreen")
        axis(2, yaxp=c(0,themax2,bamboo),cex.axis=1.5,las=1,tick=F)
        if(k==Nenviron)axis(1, at=ticks,cex.axis=1.5,labels=chm.to.analyze,tick=F)
        mtext(side=4,paste(environ_name[k],sep=""),line=2,cex=1,base_family="Arial")
 box()
 }#end of environ_name
 dev.off()
}#end of plot.type

if("s"%in%plot.type)
{
    # wd=((y-themin+base)/(themax-themin+base))*size*ratio
 wd=2
 if(is.null(allpch)) allpch=c(0,1,2,5,6,22,21,24,23,25)
  
 pdf(paste("GAPIT.Manhattan.Mutiple.Plot.symphysic",".pdf" ,sep = ""), width = 30,height=18)
 par(mfrow=c(1,1))
 par(mar = c(5,8,5,1))
 themax.y02=ceiling((ceiling(themax.y0/4))*4)

 plot(1~1,col="white",xlab="",ylab="" ,ylim=c(0,themax.y02),xlim=c(min(x),max(x)),yaxp=c(0,themax.y02,4),
    cex.axis=4, cex.lab=4,axes=FALSE,
    pch=mypch,lwd=wd,cex=s+1.3,cex.main=4)
    
        #Add a horizontal line for bonferroniCutOff
 axis(1, at=max.x,cex.axis=2,labels=rep("",length(max.x)),tick=T,lwd=2.5)
 axis(1, at=ticks,cex.axis=2,labels=chm.to.analyze,tick=F,line=1)
 axis(2, yaxp=c(0,themax.y02,4),cex.axis=2,tick=T,las=1,lwd=2.5)
 abline(h=bonferroniCutOff,lty=1,untf=T,lwd=3,col="forestgreen")
 if(plot.line)
    {
    if(!is.null(nrow(new_xz)))  
        {
        abline(v=as.numeric(new_xz[,4]),col="grey",lty=as.numeric(new_xz[,2]),untf=T,lwd=3)
        }else{
        abline(v=as.numeric(new_xz[1]),col=plot.color[as.numeric(new_xz[3])],lty=as.numeric(new_xz[2]),untf=T,lwd=3)
        }
    }
 mtext(side=2,expression(-log[10](italic(p))),line=4, cex=2.5)
 legend("top",legend=paste(environ_name,sep=""),ncol=length(environ_name),
       col="black",pch=allpch[1:Nenviron],lty=0,lwd=1,cex=2,
       bty = "o", bg = "white",box.col="white")

 for(k in 1:Nenviron)
  { 

    environ_result=read.csv(paste("GAPIT.",environ_name[k],".GWAS.Results.csv",sep=""),head=T)
    result=environ_result[,1:4]
    result=result[match(as.character(GM[,1]),as.character(result[,1])),]
    rownames(result)=1:nrow(result)
    GI.MP=result[,c(2:4)]
    borrowSlot=4
    GI.MP[,borrowSlot]=0 #Inicial as 0
    GI.MP[,5]=1:(nrow(GI.MP))
    GI.MP[,6]=1:(nrow(GI.MP)) 
    GI.MP <- GI.MP[!is.na(GI.MP[,1]),]
    GI.MP <- GI.MP[!is.na(GI.MP[,2]),]
    GI.MP[is.na(GI.MP[,3]),3]=1
    
    #Retain SNPs that have P values between 0 and 1 (not na etc)
    GI.MP <- GI.MP[GI.MP[,3]>0,]
    GI.MP <- GI.MP[GI.MP[,3]<=1,]
    #Remove chr 0 and 99
    GI.MP <- GI.MP[GI.MP[,1]!=0,]
    total_chromo=length(unique(GI.MP[,1]))
    # print(dim(GI.MP))
    if(!is.null(seqQTN))GI.MP[seqQTN,borrowSlot]=1
    numMarker=nrow(GI.MP)
    bonferroniCutOff=-log10(cutOff/numMarker)
    GI.MP[,3] <-  -log10(GI.MP[,3])
    GI.MP[,5]=1:numMarker
    y.lim <- ceiling(max(GI.MP[,3]))
    
    chm.to.analyze <- unique(GI.MP[,1])
    # print(chm.to.analyze)
    # chm.to.analyze=chm.to.analyze[order(chm.to.analyze)]
    nchr=length(chm.to.analyze)
    # print(chm.to.analyze)
    GI.MP[,6]=1:(nrow(GI.MP))
    MP_store=GI.MP
    index_GI=MP_store[,3]>=0
    MP_store <- MP_store[index_GI,]
    ticks=NULL
    lastbase=0
    for (i in chm.to.analyze)
        {
            index=(MP_store[,1]==i)
            ticks <- c(ticks, lastbase+mean(MP_store[index,2]))
            MP_store[index,2]=MP_store[index,2]+lastbase
            lastbase=max(MP_store[index,2])
        }
        
    x0 <- as.numeric(MP_store[,2])
    y0 <- as.numeric(MP_store[,3])
    z0 <- as.character(MP_store[,1])
       # convert chromosome character to number
       chor_taxa=as.character(unique(MP_store[,1]))
       chor_taxa=chor_taxa[order(as.numeric(as.character(chor_taxa)))]
       chr_letter=grep("[A-Z]|[a-z]",chor_taxa)
       if(!setequal(integer(0),chr_letter))
         {     
           z0=as.character(MP_store[,1])
           for(i in 1:(length(chor_taxa)))
              {
                index=z0==chor_taxa[i]
                z0[index]=i    
              }
          }
       z0=as.numeric(z0)
       max.x=NULL
    for (i in chm.to.analyze)
        {
            index=(MP_store[,1]==i)
            max.x=c(max.x,max(x0[index]))
        }
    max.x=c(min(x0),max.x)
    x1=sort(x0)

    position=order(y0,decreasing = TRUE)
    values=y0[position]
    if(length(values)<=DPP)
        {
         index=position[c(1:length(values))]
        }else{       
         values=sqrt(values)  #This shift the weight a little bit to the low building.
        #Handler of bias plot
         rv=runif(length(values))
         values=values+rv
         values=values[order(values,decreasing = T)]         
         theMin=min(values)
         theMax=max(values)
         range=theMax-theMin
         interval=range/DPP

         ladder=round(values/interval)
         ladder2=c(ladder[-1],0)
         keep=ladder-ladder2
         index=position[which(keep>=0)]
        }        
    x=x0[index]
    y=y0[index]
    z=z0[index]
        # print(length(x))

        #Extract QTN
        #if(!is.null(seqQTN))MP_store[seqQTN,borrowSlot]=1
        #if(!is.null(interQTN))MP_store[interQTN,borrowSlot]=2
    QTN=MP_store[which(MP_store[,borrowSlot]==1),]
        #Draw circles with same size and different thikness
    themax=ceiling(max(y))
    # themax.y02=ceiling((ceiling(themax.y0/4)+1)*4)
    # print(themax.y02)
    themin=floor(min(y))
    mypch=allpch[k]
    
   if(k!=1) par(new=T)
    
    par(new=T)
    plot(y~x,xlab="",ylab="" ,ylim=c(0,themax.y02),xlim=c(min(x),max(x)),yaxp=c(0,themax.y02,4),
    cex.axis=4, cex.lab=4, ,col=plot.color[z],axes=FALSE,
    pch=mypch,lwd=wd,cex=s+1.5,cex.main=4)
    
    if(!simulation)
       {
        abline(v=QTN[2], lty = 2, lwd=1.5, col = "grey")
        }else{
        points(QTN[,2], QTN[,3], pch=20, cex=2,lwd=2.5,col="dimgrey")
       }        
    
 }#end of environ_name
 dev.off()
}#end of plot.type

print("GAPIT.Manhattan.Mutiple.Plot has done !!!")
return(list(multip_mapP=result0,xz=new_xz))
} #end of GAPIT.Manhattan
#=============================================================================================
