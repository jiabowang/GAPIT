`GAPIT.Multiple_Synthesis` <-
function(model_store,DPP=500,chor_taxa=NULL,cutOff=0.01,band=5,seqQTN=NULL,Y.names=NULL,GM=NULL,interQTN=NULL,
    plot.style="Oceanic",plot.line=TRUE,allpch=NULL,plot.type=c("s")){
    #Object: Make a Manhattan Plot with mulitple traits
    #Output: pdfs of the Multiple Manhattan Plot
    #Authors: Zhiwu Zhang and Jiabo Wang
    # Last update: MAY 9, 2022
    ##############################################################################################
  if(!require(rgl)) install.packages("rgl")
  if(!require(rglwidget)) install.packages("rglwidget")
  library(rgl)

  Nenviron=length(model_store)*length(Y.names)
  environ_name=NULL
  new_xz=NULL
  for(i in 1:length(model_store))
  {
    for(j in 1:length(Y.names))
    {
      environ_name=c(environ_name,paste(model_store[i],".",Y.names[j],sep=""))
    }
  }
sig_pos=NULL
simulation=FALSE
    if(!is.null(seqQTN)){    
        #seqQTN=-seqQTN
        simulation=TRUE    
    }
taxa0=as.character(GM[,1])
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

if("s"%in%plot.type)
{
  # setup vals
library("plotly")
vals0=c("square","diamond","cross","x","star",
      "triangle-up","triangle-down","triangle-left","triangle-right","triangle-ne",
      "triangle-se","triangle-sw","triangle-nw","pentagon","hexagon",
      "hexagon2","octagon","circl","hexagram","star-triangle-up",
      "star-triangle-down","star-square","star-diamond","diamond-tall","diamond-wide")
n.vals=ceiling(Nenviron/length(vals0))-1
if(n.vals==0) vals=paste(vals0,"-open",sep="")
if(n.vals==1) vals=c(paste(vals0,"-open",sep=""),paste(vals0,"-dot",sep=""))
if(n.vals==2) vals=c(paste(vals0,"-open",sep=""),paste(vals0,"-dot",sep=""),vals0)
if(n.vals>2) vals=c(paste(vals0,"-open",sep=""),paste(vals0,"-dot",sep=""),vals0,paste(vals0,"-open-dot",sep=""))

# vals=vals0[1:Nenviron]
# print(vals)
if(Nenviron<=length(vals))
{
  vals=vals[1:Nenviron]
}else{
  vals=append(rep(vals,floor(Nenviron/length(vals))),vals[1:(Nenviron-length(vals))])
}
print(vals)

 x.all=NULL
 y.all=NULL
 z.all=NULL
 s.all=NULL
 taxa.all=NULL
 for(k in 1:Nenviron)
  { 

    environ_result=read.csv(paste("GAPIT.",environ_name[k],".GWAS.Results.csv",sep=""),head=T)
    result=environ_result[,1:4]
    result=result[order(result[,3]),]
    result=result[order(result[,2]),]
    result=result[match(as.character(GM[,1]),as.character(result[,1])),]
    rownames(result)=1:nrow(result)
    GI.MP=result[,c(2:4)]
    taxa0=as.character(result[,1])
    borrowSlot=4
    GI.MP[,borrowSlot]=0 #Inicial as 0
    GI.MP[,5]=1:(nrow(GI.MP))
    GI.MP[,6]=1:(nrow(GI.MP)) 
    # GI.MP <- GI.MP[!is.na(GI.MP[,1]),]
    # GI.MP <- GI.MP[!is.na(GI.MP[,2]),]
    GI.MP[is.na(GI.MP[,3]),3]=1
    
    #Retain SNPs that have P values between 0 and 1 (not na etc)
    # GI.MP <- GI.MP[GI.MP[,3]>0,]
    # GI.MP <- GI.MP[GI.MP[,3]<=1,]
    #Remove chr 0 and 99
    # GI.MP <- GI.MP[GI.MP[,1]!=0,]
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
        cut0=ceiling(-log10(0.01/length(values))/2)
        rv=runif(length(values))
        values=values+rv*(values+cut0)
        
        index=position[which(values>cut0)]
        }        
    x=x0[index]
    y=y0[index]
    z=z0[index]
    taxa=taxa0[index]
    s=rep(environ_name[k],length(x))
    x.all=append(x.all,x)
    y.all=append(y.all,y)
    z.all=append(z.all,z)
    s.all=append(s.all,s)
    taxa.all=append(taxa.all,taxa)
 }

  S.index=s.all
  for(s in 1:Nenviron)
  {
    S.index[S.index==environ_name[s]]=vals[s]
  
  }
  col.PLINK=rep(c("gray10","gray70"),ceiling(nchr/2))       
  bonferroniCutOff01=-log10(0.01/numMarker)
 c.s=z.all
 c.s[z.all%%2==0]=col.PLINK[1]
 c.s[!z.all%%2==0]=col.PLINK[2]
 
  Position=x.all
  P_value=y.all
  z.all[z.all<10]=paste("0",z.all[z.all<10],sep="")
  zz=paste("Chr_",z.all,sep="")
  #print(zz)
#  if(!require(plotly)) install.packages("plotly")
  #print("!!!!!")
  #print(head(Position))
 library(plotly)
  p <- plotly::plot_ly()%>%
   add_markers(
    # type = 'scatter',
    x = Position,
    y = P_value,
    # colorscale='Viridis',
    # reversescale =T,
    hoverinfo = "text",
    marker=list(
    symbol = S.index,
    # symbol="circl-open",
      size = 10,
      line = list(
        # color = c.s,
        color=c("steelblue"),
        width = 2)),
    # symbol=vals[k],
    text = ~paste("SNP: ", taxa.all, "<br>Chro: ", zz,"<br>Trait: ", s.all)#,
    # color = ~as.character(zz)#,
    # colors=col.PLINK
    )%>%

   # plotly::add_trace(y=bonferroniCutOff01,name = 'CutOff-0.01',color=I("red"),mode="line",width=1.4,text="")%>%
   # plotly::add_trace(y=bonferroniCutOff05,name = 'CutOff-0.05',color=I("red"),mode="line",line=list(width=1.4,dash='dot'),text="")%>%
   layout(title = "Interactive.Multiple_Synthesis.Manhattan.Plot",
                  xaxis = list(title = "Chromsome",zeroline = FALSE,showticklabels = FALSE),
                  yaxis = list (title = "-Log10(p)"))
  
    htmltools::save_html(p, paste("GAPIT.Interactive.Multiple_Synthesis.Manhattan.Plot.html",sep=""))
   
   S.uni=unique(S.index)
   T.uni=unique(s.all)
   
   q <- plotly::plot_ly()%>%
   add_markers(
    # type = 'scatter',
    x = 1,
    y = 1:length(S.uni),
    # colorscale='Viridis',
    # reversescale =T,
    hoverinfo = "text",
    marker=list(
    symbol = S.index,
    # symbol="circl-open",
      size = 10,
      line = list(
        color = c.s,
        # color=c("black","red"),
        width = 2)),
    # symbol=vals[k],
    text = ~paste("SNP: ", taxa.all, "<br>Chro: ", zz,"<br>Trait: ", s.all)#,
    # color = ~as.character(zz)#,
    # colors=col.PLINK
    )%>%

   # plotly::add_trace(y=bonferroniCutOff01,name = 'CutOff-0.01',color=I("red"),mode="line",width=1.4,text="")%>%
   # plotly::add_trace(y=bonferroniCutOff05,name = 'CutOff-0.05',color=I("red"),mode="line",line=list(width=1.4,dash='dot'),text="")%>%
   layout(title = "Interactive.Multiple_Synthesis.Manhattan.Plot",
                  xaxis = list(title = "Chromsome",zeroline = FALSE,showticklabels = FALSE),
                  yaxis = list (title = "-Log10(p)"))
  
    htmltools::save_html(p, paste("GAPIT.Interactive.Multiple_Synthesis.Manhattan.Plot.html",sep=""))
   
}#end of plot.type

print("GAPIT.Manhattan.Mutiple.Plot has done !!!")
# return(list(multip_mapP=result0,xz=new_xz))
} #end of GAPIT.Manhattan
#=============================================================================================

