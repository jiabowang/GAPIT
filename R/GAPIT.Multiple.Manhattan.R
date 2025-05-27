`GAPIT.Multiple.Manhattan` <-
function(model_store,DPP=50000,chor_taxa=NULL,cutOff=0.01,band=5,seqQTN=NULL,byTraits=FALSE,
    Y.names=NULL,GM=NULL,interQTN=NULL,WS=10e5,outpch=NULL,inpch=NULL,
    plot.style="Oceanic",plot.line=TRUE,plot.type=c("h","s","w")){
    #Object: Make a Manhattan Plot
    #Output: pdfs of the Multiple Manhattan Plot
    #Authors: Zhiwu Zhang and Jiabo Wang
    # Last update: AUG 24, 2022
    ##############################################################################################
  Nenviron=length(model_store)*length(Y.names)
  environ_name=NULL
  new_xz=NULL
  if(byTraits)
  {
    for(i in 1:length(Y.names))
    {
       for(j in 1:length(model_store))
       {
      # environ_name=c(environ_name,paste(model_store[i],".",Y.names[j],sep=""))
          environ_name=c(environ_name,paste(model_store[j],".",Y.names[i],"(NYC)", sep=""))
       }
    }
  }else{
    for(i in 1:length(model_store))
    {
       for(j in 1:length(Y.names))
       {
      # environ_name=c(environ_name,paste(model_store[i],".",Y.names[j],sep=""))
          environ_name=c(environ_name,paste(model_store[i],".",Y.names[j],"(NYC)", sep=""))
       }
    }
  }
sig_pos=NULL
simulation=FALSE
    if(!is.null(seqQTN)){    
        #seqQTN=-seqQTN
        simulation=TRUE    
    }
themax.y0=NULL
store.x=NULL
y_filter0=NULL
for(i in 1:length(environ_name))
{
  print(paste("Reading GWAS result with ",environ_name[i],sep=""))
  environ_result=read.csv(paste("GAPIT.Association.GWAS_Results.",environ_name[i],".csv",sep=""),head=T)
  num.markers=nrow(environ_result)
  environ_result=environ_result[order(environ_result[,3]),]
  environ_result=environ_result[order(environ_result[,2]),]
  environ_filter=environ_result[!is.na(environ_result[,4]),]
  themax.y=round(max(-log10(environ_filter[,4])),0)
  themax.y0=round(max(c(themax.y,themax.y0)),0)
  chm.to.analyze <- unique(environ_result[,2])
  nchr=length(chm.to.analyze)

  y_filter=environ_filter[environ_filter[,4]<(cutOff/(num.markers)),,drop=FALSE]
  traits=environ_name[i]
  # print(head(y_filter))
  # print(traits)
  if(nrow(y_filter)>0)y_filter=cbind(as.matrix(y_filter[,1:5]),traits)
  y_filter0=rbind(y_filter0,y_filter)
  # write.table(y_filter,paste("GAPIT.Filter_",environ_name[i],"_GWAS_result.txt",sep=""))

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
  # map_store=max.x
  sig_pos=append(sig_pos,as.numeric(rownames(result[result[!is.na(result[,4]),4]<(cutOff/nrow(result)),,drop=FALSE])))
}
  write.csv(y_filter0,paste("GAPIT.Association.Filter_GWAS_results.csv",sep=""),quote=FALSE)

# print(sig_pos)
#if(length(sig_pos)!=0)sig_pos=sig_pos[!duplicated(sig_pos)]
 if(length(sig_pos[!is.na(sig_pos)])>1)
 {
 # {     x_matrix=as.matrix(table(sig_pos))
 #       x_matrix=cbind(as.data.frame(rownames(x_matrix)),x_matrix)
       #print(x_matrix)
        lastbase=0
        map_store=cbind(as.data.frame(GM[,2]),as.numeric(GM[,3]))
        ticks=NULL
        # print(head(map_store))
        max.x=NULL
        for (j in unique(map_store[,1]))
        {
            index=map_store[,1]==j
            ticks <- c(ticks, lastbase+mean(map_store[index,2]))
            map_store[index,2]=as.numeric(map_store[index,2])+lastbase
            lastbase=max(as.numeric(map_store[index,2]))
            max.x=c(max.x,max(as.numeric(map_store[index,2])))
        }
        # print(ticks)
       max.x=c(min(as.numeric(map_store[,2])),max.x)
       store.x=c(store.x,as.numeric(map_store[,2]))
       # colnames(x_matrix)=c("pos","times")
       new_xz0=cbind(sig_pos,map_store[as.numeric(as.character(sig_pos)),,drop=FALSE])
       common=as.numeric(new_xz0[,3])
       scom=sort(common)
       de.sc=scom[-1]-scom[-length(scom)]
       dayu1.index=duplicated(scom)|c(abs(de.sc)<WS,FALSE)

       # print(table(dayu1.index))
       if(sum(dayu1.index)>0)
       {
       scom2=scom[dayu1.index]
       # scom2=scom2[!duplicated(scom2)]
       # print(new_xz0)
       # print(scom2)
       sc.index=as.character(new_xz0[,3])%in%scom2
       # print(table(sc.index))
       new_xz=new_xz0[sc.index,,drop=FALSE]
       # print(new_xz)
       new_xz=cbind(new_xz[,1],2,new_xz[,-1])
       new_xz[duplicated(new_xz[,4]),2]=1
       colnames(new_xz)=c("pos","times","chro","xlab")
       new_xz=new_xz[!duplicated(new_xz),]
       new_xz=as.matrix(new_xz)
       new_xz=new_xz[new_xz[,2]!="0",]
       new_xz=matrix(as.numeric(new_xz),length(as.vector(new_xz))/4,4)
       }
       # print(head(new_xz))
}else{
        lastbase=0
        map_store=cbind(as.data.frame(GM[,2]),as.numeric(GM[,3]))
        ticks=NULL
        max.x=NULL
        # print(head(map_store))
        for (j in unique(map_store[,1]))
        {
            index=map_store[,1]==j
            ticks <- c(ticks, lastbase+mean(map_store[index,2]))
            map_store[index,2]=as.numeric(map_store[index,2])+lastbase
            lastbase=max(as.numeric(map_store[index,2]))
            max.x=c(max.x,max(as.numeric(map_store[index,2])))
        }
       max.x=c(min(as.numeric(map_store[,2])),max.x)
       store.x=c(store.x,as.numeric(map_store[,2]))
}
print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2")
print(new_xz)
# setup colors
# print(head(result))
# chm.to.analyze <- unique(result[,2])
nchr=length(chm.to.analyze)
size=1 #1
ratio=10 #5
base=1 #1
numCHR=nchr
numMarker=nrow(GM)
bonferroniCutOff=-log10(cutOff/numMarker)

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
    Max.high=6*Nenviron
    if(Max.high>8)Max.high=40
    pdf(paste("GAPIT.Association.Manhattans_High",".pdf" ,sep = ""), width = 20,height=6*Nenviron)
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
       environ_result=read.csv(paste("GAPIT.Association.GWAS_Results.",environ_name[k],".csv",sep=""),head=T)
       result=environ_result[,1:4]
       result=result[order(result[,3]),]
       result=result[order(result[,2]),]
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
       # print(head(MP_store))
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
       # print(ticks)
       x1=sort(x0)
       position=order(y0,decreasing = TRUE)
       values=y0[position]
       if(length(values)<=DPP)
         {
         index=position[c(1:length(values))]
         }else{       
          # values=sqrt(values)  #This shift the weight a little bit to the low building.
        #Handler of bias plot
        cut0=ceiling(-log10(cutOff/length(values))/2)
        rv=runif(length(values))
        values=values+rv*(values+cut0)
        index=position[which(values>cut0)]
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
           pch=mypch,lwd=wd,cex=s+2.5,cex.main=2)
       mtext(side=2,expression(-log[10](italic(p))),line=3.5, cex=2.5)
       if(!simulation)
         {
          abline(v=QTN[2], lty = 2, lwd=1.5, col = "grey")}else{
          points(QTN[,2], QTN[,3], pch=21, cex=2.8,lwd=1.5,col="dimgrey")
          points(QTN[,2], QTN[,3], pch=20, cex=1.8,lwd=1.5,col="dimgrey")
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
       # print(ticks)
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
 pdf(paste("GAPIT.Association.Manhattans_Wide",".pdf" ,sep = ""), width = 16,height=8.5)
 par(mfrow=c(Nenviron,1))
 mtext.h=0.5
 size=2
 ratio=5
 for(k in 1:Nenviron)
 { 
  if(k==Nenviron)
        {#par(mfrow=c(Nenviron,1))
          # print(par())
        par(mar = c(2.5,8,0,8))
         # par(pin=c(10,((8-mtext.h)/Nenviron)+mtext.h))

        }else{
            #par(mfrow=c(Nenviron,1))
        par(mar = c(2,8,0.5,8))    
         # par(pin=c(10,(8-mtext.h)/Nenviron))

        }
  environ_result=read.csv(paste("GAPIT.Association.GWAS_Results.",environ_name[k],".csv",sep=""),head=T)
  #print(environ_result[as.numeric(new_xz[,1]),])
  result=environ_result[,1:4]
    result=result[order(result[,3]),]
    result=result[order(result[,2]),]
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
        # values=sqrt(values)  #This shift the weight a little bit to the low building.
        #Handler of bias plot
        cut0=ceiling(-log10(cutOff/length(values))/2)
        rv=runif(length(values))
        values=values+rv*(values+cut0)
      
        index=position[which(values>cut0)]
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
        # size=5
        wd=((y-themin+base)/(themax-themin+base))*size*ratio
        # wd=0.5
        s=size-wd/ratio/2
        mypch=1
        bamboo=4
        # plot(y~x,xlab="",ylab="" ,ylim=c(0,themax),
        #     cex.axis=4, cex.lab=4, ,col=plot.color[z],axes=FALSE,type = "p",pch=mypch,lwd=0.5,cex=0.7,cex.main=2)
        plot(y~x,xlab="",ylab="" ,ylim=c(0,themax2),xlim=c(min(x),max(x)),
           cex.axis=4, cex.lab=4, ,col=plot.color[z],axes=FALSE,type = "p",
           pch=mypch,lwd=wd,cex=s,cex.main=2)
        mtext(side=2,expression(-log[10](italic(p))),line=3, cex=1)
        if(plot.line)
        {
          if(!is.null(nrow(new_xz)))  {abline(v=as.numeric(new_xz[,4]),col="grey",lty=as.numeric(new_xz[,2]),untf=T,lwd=2)
             }else{abline(v=as.numeric(new_xz[1]),col=plot.color[as.numeric(new_xz[3])],lty=as.numeric(new_xz[2]),untf=T,lwd=2)
             }
        }
        if(!simulation){abline(v=QTN[2], lty = 2, lwd=1.5, col = "grey")}else{
          # print("$$$")
          points(QTN[,2], QTN[,3], type="p",pch=21, cex=2.8,lwd=1.5,col="dimgrey")
          points(QTN[,2], QTN[,3], type="p",pch=20, cex=1.5,lwd=1.5,col="dimgrey")
          }
        #Add a horizontal line for bonferroniCutOff
        abline(h=bonferroniCutOff,lty=1,untf=T,lwd=1,col="forestgreen")
        axis(2, yaxp=c(0,themax2,bamboo),cex.axis=1.5,las=1,tick=F)
        if(k==Nenviron)axis(1, at=ticks,cex.axis=1.5,line=0.001,labels=chm.to.analyze,tick=F)
        mtext(side=4,paste(environ_name[k],sep=""),line=2,cex=1,base_family="Arial")
 box()
 }#end of environ_name
 dev.off()
}#end of plot.type

if("s"%in%plot.type)
{
    # wd=((y-themin+base)/(themax-themin+base))*size*ratio
 wd=2
 if(is.null(outpch))
 {
   allpch0=c(0,1,2,5,6)
 }else{
   allpch0=outpch
 }
 if(is.null(inpch))
 {
   add.pch=c("+","*","-","#","<",">","^","$","=","|","?",as.character(1:9),letters[1:26],LETTERS[1:26]) 
 }else{
   add.pch=inpch
 }
 n.vals=ceiling(Nenviron/length(allpch0))-1
 s=size-wd/ratio/2
 DPP=500
 


 pdf(paste("GAPIT.Association.Manhattans_Symphysic",".pdf" ,sep = ""), width = 30,height=18)
 par(mfrow=c(1,1))
 par(mar = c(5,8,5,1))
 themax.y02=ceiling((ceiling(themax.y0/4))*4)

 plot(1~1,col="white",xlab="",ylab="" ,ylim=c(0,themax.y02),xlim=c(min(store.x,na.rm=TRUE),max(store.x,na.rm=TRUE)),yaxp=c(0,themax.y02,4),
    cex.axis=4, cex.lab=4,axes=FALSE,
    pch=1,cex.main=4)
 # print(ticks)   
        #Add a horizontal line for bonferroniCutOff
 axis(1, at=max.x,cex.axis=2,labels=rep("",length(max.x)),tick=T,lwd=2.5)
 axis(1, at=ticks,cex.axis=2,labels=chm.to.analyze,tick=F,line=1)
 axis(2, yaxp=c(0,themax.y02,4),cex.axis=2,tick=T,las=1,lwd=2.5)
 if(!is.null(cutOff))abline(h=bonferroniCutOff,lty=1,untf=T,lwd=3,col="forestgreen")
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
 # legend("top",legend=paste(environ_name,sep=""),ncol=length(environ_name),
 #       col="black",pch=allpch[1:Nenviron],lty=0,lwd=1,cex=2,
 #       bty = "o", bg = "white",box.col="white")
 # step.vals=0
 # if(1>2){
 for(k in 1:Nenviron)
  { 
    step.vals=ceiling(k/length(allpch0))-1

    environ_result=read.csv(paste("GAPIT.Association.GWAS_Results.",environ_name[k],".csv",sep=""),head=T)
    result=environ_result[,1:4]
    result=result[order(result[,3]),]
    result=result[order(result[,2]),]
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
          # values=sqrt(values)  #This shift the weight a little bit to the low building.
        #Handler of bias plot
        cut0=ceiling(-log10(cutOff/length(values))/2)
        rv=runif(length(values))
        values=values+rv*(values+cut0)

        index=position[which(values>cut0)]
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
    mypch=allpch0[k-step.vals*length(allpch0)]
    
   # if(k!=1) par(new=T)
    
    par(new=T)
    plot(y~x,xlab="",ylab="" ,ylim=c(0,themax.y02),xlim=c(min(x),max(x)),yaxp=c(0,themax.y02,4),
    cex.axis=4, cex.lab=4,col=plot.color[z],axes=FALSE,
    pch=mypch,lwd=1,cex=s+2.5,cex.main=4)
    if(step.vals!=0)
    {
      points(y~x,pch=add.pch[step.vals],col=plot.color[z],cex=s+0.5,cex.main=4)
    }
    if(!simulation)
       {
        abline(v=QTN[2], lty = 2, lwd=1.5, col = "grey")
        }else{
        points(QTN[,2], QTN[,3], pch=20, cex=2,lwd=2.5,col="dimgrey")
       }        
    
 }#end of environ_name
 dev.off()
 # }
 ## Plot legend
 nchar.traits=1.5
 # environ_name=paste(environ_name,"1234",sep="")
 nchar0=max(nchar(environ_name))
 # print(nchar0)
 if(Nenviron>5)
 {
  yourpch=c(rep(allpch0,n.vals),allpch0[1:(Nenviron-length(allpch0)*n.vals)])
 }else{
  yourpch=allpch0[1:Nenviron]
 }
  yourpch2=NULL
  for(pp in 1:n.vals)
  {
    yourpch2=c(yourpch2,rep(add.pch[pp],length(allpch0)))
  }
  # yourpch2=c(rep(allpch0,n.vals),allpch0[Nenviron-length(allpch0)*n.vals])
  if(Nenviron>5){
  yourpch2=yourpch2[1:(Nenviron-length(allpch0))]
  yourpch2=c(rep(NA,length(allpch0)),yourpch2)
  }
 
 max.row=25
 max.pch=ifelse(Nenviron<max.row,Nenviron,max.row)
 n.col.pch=ceiling(Nenviron/max.row)
 ratio.cex=ceiling(Nenviron/5)
  if(ratio.cex>5)ratio.cex=5
  cex.Ne=ratio.cex*(0.05*ratio.cex+0.35) #the size of cex
  if(ratio.cex<3)cex.Ne=1
  c.t.d=c(0.5,1,1,1,1.5)[ratio.cex] # the different between size of cex and text
  cex.di=0.3*ratio.cex # the different size between cex and signal
  text.di=c(.02,0.01,0.01,0.02,0.02)[ratio.cex] #the different distance between cex and text

  high.Ne=2*ratio.cex # the total highth of figure
  cex.betw=c(.38,.5,.7,0.9,1)[ratio.cex] # distance between cexes
  x.di=c(1.12,1.09,0.5,0.8,1.3)[ratio.cex]/2 # distance between markers in x axis
  # print(nchar0)
  # x.di0=c(0)
  if(n.col.pch>1){
    text.di=(0.01*nchar0)/3+0.02
    # x.di=0.52*n.col.pch
    x.di=(0.1*(ceiling(nchar0/5)-1)+(n.col.pch-1)*0.1)*ceiling(nchar0/5)#*n.col.pch
  }
  # print(x.di)
 # if(Nenviron>5){
 #  cex.Ne=3
 #  cex.di=1.5
 #  text.di=.02
 #  high.Ne=10
 #  cex.betw=0.9
 #  }else{
 #  cex.Ne=1
 #  cex.di=0.3
 #  high.Ne=Nenviron/2 
 #  text.di=.02
 #  cex.betw=0.9
 #  }


 write.csv(environ_name,"GAPIT.Association.Manhattans_Symphysic_Traitsnames.csv",quote=FALSE)
 pdf(paste("GAPIT.Association.Manhattans_Symphysic_Legend",".pdf" ,sep = ""), width = 4+(x.di*(n.col.pch+1)),height=high.Ne)
 par(mfrow=c(1,1))
 par(mar = c(cex.Ne+1,2,cex.Ne+1,2))
 # print(length(yourpch))
 # print(length(yourpch2))
 plot(0,0,xlab="",ylab="" ,axes=FALSE,
  xlim=c(0,x.di*(n.col.pch)),ylim=c(0,max.pch),col="white")
 for(kk in 1:n.col.pch)
 {
 par(new=T)

 if(kk==n.col.pch)
 {
  if(n.col.pch==1)
  {
  # print(kk)
  max.pch2=Nenviron-(n.col.pch-1)*max.row
  
  plot(rep(0,max.pch2),(max.pch:(max.pch-max.pch2+1))*cex.betw,xlab="",ylab="" ,axes=FALSE,col="black",
  xlim=c(0,x.di*(n.col.pch)),ylim=c(0,max.pch),lwd=1,cex=cex.Ne,
  pch=yourpch[((kk-1)*max.row+1):Nenviron])
  if(Nenviron>5) points(rep(0,max.pch2),((max.pch):(max.pch-max.pch2+1))*cex.betw,
  xlim=c(0,x.di*(n.col.pch)),ylim=c(0,max.pch),lwd=1,cex=cex.Ne-cex.di,
  pch=yourpch2[((kk-1)*max.row+1):Nenviron])
  text(rep((0+text.di),max.pch2),(max.pch:(max.pch-max.pch2+1))*cex.betw,labels=environ_name[((kk-1)*max.row+1):Nenviron],pos=4,cex=cex.Ne-c.t.d)
  }else{
  # print(kk)
  max.pch2=Nenviron-(n.col.pch-1)*max.row
  
  plot(rep((kk-1)*x.di,max.pch2),(max.pch:(max.pch-max.pch2+1))*cex.betw,xlab="",ylab="" ,axes=FALSE,col="black",
  xlim=c(0,x.di*(n.col.pch)),ylim=c(0,max.pch),lwd=1,cex=cex.Ne,
  pch=yourpch[((kk-1)*max.row+1):Nenviron])
  if(Nenviron>5) points(rep((kk-1)*x.di,max.pch2),((max.pch):(max.pch-max.pch2+1))*cex.betw,
  xlim=c(0,x.di*(n.col.pch)),ylim=c(0,max.pch),lwd=1,cex=cex.Ne-cex.di,
  pch=yourpch2[((kk-1)*max.row+1):Nenviron])
  text(rep(((kk-1)*x.di+text.di),max.pch2),(max.pch:(max.pch-max.pch2+1))*cex.betw,labels=environ_name[((kk-1)*max.row+1):Nenviron],pos=4,cex=cex.Ne-c.t.d)
    
  }
 }else{
  if(kk==1)
  {
  print(kk)
  plot(rep(0,max.pch),(max.pch:1)*cex.betw,xlab="",ylab="" ,axes=FALSE,
  xlim=c(0,x.di*(n.col.pch)),ylim=c(0,max.pch),lwd=1,cex=cex.Ne,
  pch=yourpch[((kk-1)*max.row+1):((kk-1)*max.row+max.row)])
  if(Nenviron>5) points(rep(0,max.pch),((max.pch):1)*cex.betw,
  xlim=c(1,x.di*(n.col.pch)),ylim=c(0,max.pch),lwd=1,cex=cex.Ne-cex.di,
  pch=yourpch2[((kk-1)*max.row+1):((kk-1)*max.row+max.row)])
  text(rep((0+text.di),max.pch),(max.pch:1)*cex.betw,labels=environ_name[((kk-1)*max.row+1):((kk-1)*max.row+max.row)],pos=4,cex=cex.Ne-c.t.d)
  }else{
  print(kk)
  plot(rep((kk-1)*x.di,max.pch),(max.pch:1)*cex.betw,xlab="",ylab="" ,axes=FALSE,
  xlim=c(0,x.di*(n.col.pch)),ylim=c(0,max.pch),lwd=1,cex=cex.Ne,
  pch=yourpch[((kk-1)*max.row+1):((kk-1)*max.row+max.row)])
  if(Nenviron>5) points(rep((kk-1)*x.di,max.pch),((max.pch):1)*cex.betw,
  xlim=c(0,x.di*(n.col.pch)),ylim=c(0,max.pch),lwd=1,cex=cex.Ne-cex.di,
  pch=yourpch2[((kk-1)*max.row+1):((kk-1)*max.row+max.row)])
  text(rep(((kk-1)*x.di+text.di),max.pch),(max.pch:1)*cex.betw,labels=environ_name[((kk-1)*max.row+1):((kk-1)*max.row+max.row)],pos=4,cex=cex.Ne-c.t.d)
   
  }

 }
}#end of plot.type
 dev.off()

}
print("GAPIT.Association.Manhattans has done !!!")
return(list(multip_mapP=result0,xz=new_xz))
} #end of GAPIT.Manhattan
#=============================================================================================
