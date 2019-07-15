`GAPIT.GS.Visualization` <-
function(gsBLUP = gsBLUP, BINS=BINS, name.of.trait = name.of.trait){
#Object: To build heat map to show distribution of BLUP and PEV
#Output: pdf
#Authors: Zhiwu Zhang 
# Last update: May 15, 2011 
##############################################################################################
nBin=BINS

BLUP= gsBLUP[,5]
PEV = gsBLUP[,6]

if(BLUP[1]=="NaN"){
  warning ("It was not converged. BLUP was not created!")
}
if(BLUP[1]!="NaN" )
{


BLUP.max=try(max(BLUP))
BLUP.min=try(min(BLUP))
if(inherits(BLUP.max, "try-error"))  return()

  range.BLUP=BLUP.max-BLUP.min
  range.PEV=max(PEV)-min(PEV)
  
  interval.BLUP=range.BLUP/nBin
  interval.PEV=range.PEV/nBin
  
  
  bin.BLUP=floor(BLUP/max(BLUP)*nBin)*max(BLUP)/nBin
  bin.PEV=floor(PEV/max(PEV)*nBin)*max(PEV)/nBin
  
  
  distinct.BLUP=unique(bin.BLUP)
  distinct.PEV=unique(bin.PEV)
  
  if((length(distinct.BLUP)<2)  | (length(distinct.PEV)<2) ) return() #nothing to plot
  
  Position.BLUP=match(bin.BLUP,distinct.BLUP,nomatch = 0)
  Position.PEV=match(bin.PEV,distinct.PEV,nomatch = 0)
  
  value=matrix(1,length(Position.BLUP))
  KG<- (tapply(as.numeric(value), list(Position.BLUP, Position.PEV), sum))
  
  rownames(KG)=round(distinct.BLUP, digits = 4)
  colnames(KG)=round(distinct.PEV, digits = 4)
  
  #Sort the rows and columns in order from smallest to largest
  
  rownames(KG) <- rownames(KG)[order(as.numeric(rownames(KG)))]
  colnames(KG) <- colnames(KG)[order(as.numeric(colnames(KG)))]
  rownames(KG) <- round(as.numeric(rownames(KG)))
  colnames(KG) <- round(as.numeric(colnames(KG)))
  #write.table(KG, "Input_Matrix_for_GS_Heat_Map.txt", quote = FALSE, sep = "\t", row.names = FALSE,col.names = FALSE)

  pdf(paste("GAPIT.", name.of.trait,".GPS.BLUPvsPEV", ".pdf", sep = ""),width = 9)
  #par(mfrow = c(1,1), mar = c(1,1,5,5), lab = c(5,5,7))
  par(mar = c(5,5,6,5))
  
  nba_heatmap <- heatmap.2(KG, Rowv=NA, Colv=NA,  col =  rev(heat.colors(256)), #  scale="column", 
  xlab = "PEV", ylab = "BLUP", main = " ", scale="none", symkey=FALSE, trace="none")

  #nba_heatmap <- heatmap.2(KG,  cexRow =.2, cexCol = 0.2, scale="none", symkey=FALSE, trace="none" )
 
  
  #cexRow =0.9, cexCol = 0.9)
  dev.off() 
}
#print("GAPIT.GS.Visualization accomplished successfully!")

}   #GAPIT.GS.Visualization ends here
#=============================================================================================

