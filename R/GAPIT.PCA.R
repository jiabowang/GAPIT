`GAPIT.PCA` <-
function(X,taxa, PC.number = min(ncol(X),nrow(X)),radius=1,
  file.output=TRUE,PCA.total=0,PCA.col=NULL,
  PCA.3d=FALSE,PCA.legend=NULL){
# Object: Conduct a principal component analysis, and output the prinicpal components into the workspace,
#         a text file of the principal components, and a pdf of the scree plot
# Authors: Alex Lipka and Hyun Min Kang
# Last update: May 31, 2011  
############################################################################################## 
#Conduct the PCA 
print("Calling prcomp...")
PCA.X <- stats::prcomp(X)
eigenvalues <- PCA.X$sdev^2
evp=eigenvalues/sum(eigenvalues)
nout=min(10,length(evp))
xout=1:nout
if(is.null(PCA.col)) PCA.col="red"
# if(!is.null(PCA.legend)) PCA.col0=
print("Joining taxa...")
#Extract number of PCs needed
PCs <- cbind(taxa,as.data.frame(PCA.X$x))
if(file.output) utils::write.table(PCs[,1:(PCA.total+1)], "GAPIT.Genotype.PCA.csv", quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)

# if(file.output) utils::write.table(PCA.X$rotation[,1:PC.number], "GAPIT.Genotype.PCA_loadings.csv", quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)

if(file.output) utils::write.table(eigenvalues, "GAPIT.Genotype.PCA_eigenvalues.csv", quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)

print("Creating PCA graphs...")
#Create a Scree plot 
if(file.output & PC.number>1) {
grDevices::pdf("GAPIT.Genotype.PCA_eigenValue.pdf", width = 12, height = 12)
  graphics::par(mar=c(5,5,4,5)+.1,cex=2)
  #par(mar=c(10,9,9,10)+.1)
  plot(xout,eigenvalues[xout],type="b",col="blue",xlab="Principal components",ylab="Variance")
  graphics::par(new=TRUE)
  plot(xout,evp[xout]*100,type="n",col="red",xaxt="n",yaxt="n",xlab="",ylab="")
  graphics::axis(4)
  graphics::mtext("Percentage (%)",side=4,line=3,cex=2)
grDevices::dev.off()

grDevices::pdf("GAPIT.Genotype.PCA_2D.pdf", width = 8, height = 8)
graphics::par(mar = c(5,5,5,5),xpd=TRUE)
maxPlot=min(as.numeric(PC.number[1]),3)

for(i in 1:(maxPlot-1))
{
   for(j in (i+1):(maxPlot))
   {
      plot(PCA.X$x[,i],PCA.X$x[,j],xlab=paste("PC",i," (evp=",round(evp[i],4)*100,"%)",sep=""),ylab=paste("PC",j," (evp=",round(evp[j],4)*100,"%)",sep=""),pch=19,col=PCA.col,cex.axis=1.3,cex.lab=1.4, cex.axis=1.2, lwd=2,las=1)
      if(!is.null(PCA.legend)) legend(as.character(PCA.legend$pos),legend=PCA.legend$taxa,pch=19,col=PCA.legend$col,
           ncol=PCA.legend$ncol,box.col="white",bty = "n", bg = par("bg"),inset=-0.05)
   }
}
grDevices::dev.off()

#output 3D plot
if(PCA.3d==TRUE)
{   
  
        mycols=cbind(taxa,PCA.col)
        colnames(mycols)=c("taxa","color")
        write.csv(mycols,"color_file.csv",quote=F,row.names=F)
        GAPIT.3D.PCA.python("color_file.csv")
#    if(!require(scatterplot3d)) install.packages("scatterplot3d")
#    library(scatterplot3d)

    grDevices::pdf("GAPIT.Genotype.PCA_3D.pdf", width = 7, height = 7)
    graphics::par(mar = c(5,5,5,5),xpd=TRUE)
    scatterplot3d::scatterplot3d(PCA.X$x[,1],
                  PCA.X$x[,2],
                  PCA.X$x[,3],
                  xlab = paste("PC",1," (evp=",round(evp[1],4)*100,"%)",sep=""),
                  ylab = paste("PC",2," (evp=",round(evp[2],4)*100,"%)",sep=""),
                  zlab = paste("PC",3," (evp=",round(evp[3],4)*100,"%)",sep=""),
                  pch = 20,
                  color = PCA.col,
                  col.axis = "blue",
                  cex.symbols = 1,
                  cex.lab = 1.4,
                  cex.axis = 1.2,
                  lwd = 3,
                  angle = 55,
                  scale.y = 0.7)
    if(!is.null(PCA.legend)) legend(as.character(PCA.legend$pos),legend=PCA.legend$taxa,pch=19,col=PCA.legend$col,
      ncol=PCA.legend$ncol,box.col="white",bty = "n", bg = par("bg"),inset=-0.05)
    grDevices::dev.off()
  }#PCA.3d==TRUE
}#file.output & PC.number>1

#Remove duplicate (This is taken care by QC)
#PCs.unique <- unique(PCs[,1])
#PCs <-PCs[match(PCs.unique, PCs[,1], nomatch = 0), ]



print("Exporting PCs...")
#Write the PCs into a text file

#Return the PCs
return(list(PCs=PCs,EV=PCA.X$sdev^2,nPCs=NULL))
}

