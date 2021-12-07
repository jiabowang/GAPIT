`GAPIT.LDplot`<-function(G=NULL,LD.chromosome=NULL,LD.location=NULL,LD.range1=1000000,LD.range2=1000000,marker=NULL){
    #Object: Make a LD Plot using significant maker within a window
    #Output: A pdf of the LD Plot
    #Authors:  Jiabo Wang
    # Last update: Mar 12, 2021
    ##############################################################################################
	chromosome=(G[,3]==LD.chromosome[1])
    bp=as.numeric(as.vector(G[,4]))
    if(is.null(LD.range1)&is.null(LD.range2)) stop("Please provide LD.range for 1 or/and 2 !!!")
    if(is.null(LD.range2)&!is.null(LD.range1))LD.range2=LD.range1
    if(is.null(LD.range1)&!is.null(LD.range2))LD.range1=LD.range2
    deviation=bp-as.numeric(as.vector(LD.location[1]))
    location1=deviation< as.numeric(as.vector(LD.range1[1])  )&deviation>=0
    location2=deviation*(-1)< as.numeric(as.vector(LD.range2[1])  )&deviation<=0
    location=location1|location2
    index=chromosome&location
    GLD=G[index,]

 if(nrow(GLD)>1)
    {
    print("Plot LD...")
    hapmapgeno= data.frame(as.matrix(t(GLD[,-c(1:11)])))
    # print(hapmapgeno)
    hapmapgeno[hapmapgeno=="NN"]=NA
    hapmapgeno[hapmapgeno=="XX"]=NA
    hapmapgeno[hapmapgeno=="--"]=NA
    hapmapgeno[hapmapgeno=="++"]=NA
    hapmapgeno[hapmapgeno=="//"]=NA
    print(dim(GLD))
    # print(GLD[,1:4])
    LDdist=as.numeric(as.vector(GLD[,4]))
    LDsnpName=GLD[,1]
    colnames(hapmapgeno)=LDsnpName
    sigsnp=match(LD.location,LDdist)
    print(sigsnp)

#Prune SNM names
#LDsnpName=LDsnpName[GAPIT.Pruning(LDdist,DPP=7)]
    # LDsnpName=as.character(LDsnpName[c(1,sigsnp,length(LDsnpName))]) #keep the first and last snp names only
    LDsnpName=as.character(LDsnpName[c(sigsnp)]) #keep the first and last snp names only
    # LDsnpName=as.character(LDsnpName)
    print(LDsnpName)
    LDsnpName=gsub("SNP-","",LDsnpName)
    color.rgb <- grDevices::colorRampPalette(rev(c("snow","red")),space="rgb")
    
    # print(color.rgb)
    print("Getting genotype object")
    # print(hapmapgeno)
    LDsnp=genetics::makeGenotypes(hapmapgeno,sep="",method=genetics::as.genotype)   #This need to be converted to genotype object
    # print(LDsnp)
    # print(LDdist)
    print("Calling LDheatmap...")
#pdf(paste("GAPIT.LD.pdf",sep=""), width = 12, height = 12)
    grDevices::pdf(paste("GAPIT.LD.chromosom",LD.chromosome,"(",round(max(0,LD.location-LD.range2)/1000000),"_",round((LD.location+LD.range1)/1000000),"Mb)",".pdf",sep=""), width = 12, height = 12)
    graphics::par(mar = c(25,25,25,25))

    # MyHeatmap <- try(LDheatmap(LDsnp, LDdist, LDmeasure="r", add.map=TRUE,flip=TRUE,
    # color=color.rgb(20), name="myLDgrob", add.key=TRUE,geneMapLabelY=0.1) )  
    MyHeatmap <- LDheatmap::LDheatmap(LDsnp, LDdist, flip=TRUE,
    color=color.rgb(20),SNP.name = LDsnpName, name="myLDgrob" ) 
    
    
  #Modify the plot
#      library(grid)
      grid::grid.edit(grid::gPath("myLDgrob", "geneMap","SNPnames"), gp = grid::gpar(cex=0.35,col="blue")) #Edit SNP name
      
      if(!is.null(marker))
      {LDheatmap::LDheatmap.highlight(MyHeatmap, i = marker[1], j=marker[2], col = "blue",lwd=1.5)
      }
      # LDheatmap(MyHeatmap, SNP.name = LDsnpName)
      # grid.edit(gPath("myLDgrob","heatMap","heatmap"),gp=gpar(col="white",lwd=1))
      # grid.edit(gPath("myLDgrob", "Key", "title"), gp=gpar(cex=.5, col="blue"))  #edit key title size and color
      # grid.edit(gPath("myLDgrob", "heatMap", "title"), gp=gpar(just=c("center","bottom"), cex=0.8, col="black")) #Edit gene map title
      
    grDevices::dev.off()
    print("LD heatmap crated")
    }else{ # alternative of if(nrow(GLD)>1)
    print("Warning: There are less than two SNPs on the region you sepcified. No LD plot!")
    } #end of #if(nrow(GLD)>1)
}

