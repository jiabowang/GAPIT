`GAPIT.Multiple.Manhattan` <-
function(model_store,DPP=50000,cutOff=0.01,band=5,seqQTN=NULL,Y=NULL,GM=NULL,interQTN=NULL,plot.style="Oceanic",plot.line=FALSE){
    #Object: Make a Manhattan Plot
    #Options for plot.type = "Separate_Graph_for_Each_Chromosome" and "Same_Graph_for_Each_Chromosome"
    #Output: A pdf of the Manhattan Plot
    #Authors: Alex Lipka, Zhiwu Zhang, Meng Li and Jiabo Wang
    # Last update: Oct 10, 2016
	#Add r2 between candidata SNP and other markers in on choromosome
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
for(i in 1:length(environ_name))
{
  print(paste("Reading GWAS result with ",environ_name[i],sep=""))
  environ_result=read.csv(paste("GAPIT.",environ_name[i],".GWAS.Results.csv",sep=""),head=T)
  environ_filter=environ_result[!is.na(environ_result[,4]),]
  y_filter=environ_filter[environ_filter[,4]<(cutOff/(nrow(environ_filter))),]
  write.table(y_filter,paste("Filter_",environ_name[i],"_GWAS_result.txt",sep=""))

  result=environ_result[,1:4]

  result=result[match(as.character(GM[,1]),as.character(result[,1])),]
  # result=result[order(result[,2]),]
  # result=result[order(result[,1]),]
  #print(head(result))
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
       map_store=as.matrix(cbind(as.numeric(GM[,2]),as.numeric(as.vector(GM[,3]))))
       #print(head(map_store))
       #print(as.numeric(map_store[,3]))
        for (j in unique(map_store[,1]))
        {
            index=map_store[,1]==j
            #print(table(index))
            map_store[index,2]=as.numeric(map_store[index,2])+lastbase
            lastbase=max(as.numeric(map_store[index,2]))
            #print(lastbase)
        }
       
       colnames(x_matrix)=c("pos","times")
       #colnames(xz)=c("pos","col")
       new_xz=cbind(x_matrix,map_store[as.numeric(as.character(x_matrix[,1])),])
       #new_xz[,4]=0
       colnames(new_xz)=c("pos","times","chro","xlab")
       
       new_xz=new_xz[!duplicated(new_xz),]
       new_xz[new_xz[,2]>=3,2]=3
       new_xz[,2]=4-new_xz[,2]
       new_xz[new_xz[,2]==3,2]=0

       new_xz=as.matrix(new_xz)
       new_xz=new_xz[new_xz[,2]!="0",]
       new_xz=matrix(new_xz,length(as.vector(new_xz))/4,4)
       #print(new_xz)
       plot.line=TRUE
       #print(new_xz)
}
#print(as.numeric(new_xz[,4]))
#print(head(result0))
# print(as.numeric(new_xz[,1]))
pdf(paste("GAPIT.Manhattan.Mutiple.Plot",colnames(result0)[-c(1:3)],".pdf" ,sep = ""), width = 20,height=6*Nenviron)
par(mfrow=c(Nenviron,1))
for(k in 1:Nenviron)
{ if(k==Nenviron){#par(mfrow=c(Nenviron,1))
        par(mar = c(3,8,1,8))
        }else{
            #par(mfrow=c(Nenviron,1))
        par(mar = c(0,8,1,8))
        
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
    total_chromo=max(GI.MP[,1])
    # print(dim(GI.MP))
    if(!is.null(seqQTN))GI.MP[seqQTN,borrowSlot]=1
    numMarker=nrow(GI.MP)
    bonferroniCutOff=-log10(cutOff/numMarker)
    GI.MP[,3] <-  -log10(GI.MP[,3])
    GI.MP[,5]=1:numMarker
    y.lim <- ceiling(max(GI.MP[,3]))
    
    chm.to.analyze <- unique(GI.MP[,1])
    chm.to.analyze=chm.to.analyze[order(chm.to.analyze)]
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
        z0 <- as.numeric(MP_store[,1])
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
        size=1 #1
        ratio=10 #5
        base=1 #1
        numCHR=nchr
        themax=ceiling(max(y))
        themin=floor(min(y))
        wd=((y-themin+base)/(themax-themin+base))*size*ratio
        s=size-wd/ratio/2
        ncycle=ceiling(nchr/5)
        ncolor=5*ncycle
        ncolor=band*ncycle

        thecolor=seq(1,nchr,by= ncycle)
        mypch=1
        #plot.color= rainbow(ncolor+1)
        col.Rainbow=rainbow(ncolor+1)     
        col.FarmCPU=rep(c("#CC6600","deepskyblue","orange","forestgreen","indianred3"),ceiling(numCHR/5))
        col.Rushville=rep(c("orangered","navyblue"),ceiling(numCHR/2))    
        col.Congress=rep(c("deepskyblue3","firebrick"),ceiling(numCHR/2))
        col.Ocean=rep(c("steelblue4","cyan3"),ceiling(numCHR/2))    
        col.PLINK=rep(c("gray10","gray70"),ceiling(numCHR/2))     
        col.Beach=rep(c("turquoise4","indianred3","darkolivegreen3","red","aquamarine3","darkgoldenrod"),ceiling(numCHR/5))
        #col.Oceanic=rep(c( '#EC5f67',  '#F99157',  '#FAC863',  '#99C794',  '#5FB3B3',  '#6699CC',  '#C594C5',  '#AB7967'),ceiling(numCHR/8))
        #col.Oceanic=rep(c( '#EC5f67',    '#FAC863',  '#99C794',    '#6699CC',  '#C594C5',  '#AB7967'),ceiling(numCHR/6))
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
    
        #plot.color=rep(c( '#EC5f67',    '#FAC863',  '#99C794',    '#6699CC',  '#C594C5'),ceiling(ncolor/5))

            plot(y~x,xlab="",ylab="" ,ylim=c(0,themax),
            cex.axis=4, cex.lab=4, ,col=plot.color[z],axes=FALSE,type = "p",pch=mypch,lwd=wd,cex=s+2.5,cex.main=4)
            mtext(side=2,expression(-log[10](italic(p))),line=3, cex=2.5)
        #Label QTN positions
        #print(head(QTN))
        #print(head(interQTN))
          if(!simulation){abline(v=QTN[2], lty = 2, lwd=1.5, col = "grey")}else{
            #print("$$$$$$")
          points(QTN[,2], QTN[,3], pch=20, cex=2.5,lwd=2.5,col="black")
          #points(interQTN[,2], interQTN[,3], type="p",pch=8, cex=1,lwd=1.5,col="dimgrey")
          }
        
        #}
        if(plot.line){
          #print(x)
          #print(as.numeric(new_xz[,2]))
          # if(!is.null(nrow(new_xz)))  {abline(v=as.numeric(new_xz[,4]),col=plot.color[as.numeric(new_xz[,3])],lty=as.numeric(new_xz[,2]),untf=T,lwd=3)
          if(!is.null(nrow(new_xz)))  {abline(v=as.numeric(new_xz[,4]),col="grey",lty=as.numeric(new_xz[,2]),untf=T,lwd=3)
             }else{abline(v=as.numeric(new_xz[1]),col=plot.color[as.numeric(new_xz[3])],lty=as.numeric(new_xz[2]),untf=T,lwd=3)
             }
        }
        #Add a horizontal line for bonferroniCutOff
        abline(h=bonferroniCutOff,lty=2,untf=T,lwd=3,col="red")
        axis(2, xaxp=c(1,themax,5),cex.axis=2.5,tick=F)
        if(k==Nenviron)axis(1, at=ticks,cex.axis=2.7,labels=chm.to.analyze,tick=F)
        mtext(side=4,paste(environ_name[k],sep=""),line=3,cex=2.5)
box()



}#end of environ_name

dev.off()
print("GAPIT.Manhattan.Mutiple.Plot has done !!!")
return(list(multip_mapP=result0,xz=new_xz))
} #end of GAPIT.Manhattan
#=============================================================================================
