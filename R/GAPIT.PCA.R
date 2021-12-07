`GAPIT.PCA` <-
function(X,taxa, PC.number = min(ncol(X),nrow(X)),file.output=TRUE,PCA.total=0,PCA.col=NULL,PCA.3d=FALSE){
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

print("Creating PCA graphs...")
#Create a Scree plot 
if(file.output & PC.number>1) {
grDevices::pdf("GAPIT.PCA.eigenValue.pdf", width = 12, height = 12)
  graphics::par(mar=c(5,5,4,5)+.1,cex=2)
  #par(mar=c(10,9,9,10)+.1)
  plot(xout,eigenvalues[xout],type="b",col="blue",xlab="Principal components",ylab="Variance")
  graphics::par(new=TRUE)
  plot(xout,evp[xout]*100,type="n",col="red",xaxt="n",yaxt="n",xlab="",ylab="")
  graphics::axis(4)
  graphics::mtext("Percentage (%)",side=4,line=3,cex=2)
grDevices::dev.off()

grDevices::pdf("GAPIT.PCA.2D.pdf", width = 8, height = 8)
graphics::par(mar = c(5,5,5,5))
maxPlot=min(as.numeric(PC.number[1]),3)

for(i in 1:(maxPlot-1)){
for(j in (i+1):(maxPlot)){
plot(PCA.X$x[,i],PCA.X$x[,j],xlab=paste("PC",i,sep=""),ylab=paste("PC",j,sep=""),pch=19,col=PCA.col,cex.axis=1.3,cex.lab=1.4, cex.axis=1.2, lwd=2,las=1)

}
}
grDevices::dev.off()

#output 3D plot
if(PCA.3d==TRUE)
{   
  if(1>2)
  {
#  if(!require(lattice)) install.packages("lattice")
#   library(lattice)
   pca=as.data.frame(PCA.X$x)
   
   grDevices::png(file="example%03d.png", width=500, heigh=500)
    for (i in seq(10, 80 , 1)){
        print(lattice::cloud(PC1~PC2*PC3,data=pca,screen=list(x=i,y=i-40),pch=20,color="red",
        col.axis="blue",cex=1,cex.lab=1.4, cex.axis=1.2,lwd=3))
        }
    grDevices::dev.off()
    system("convert -delay 40 *.png GAPIT.PCA.3D.gif")
    
    # cleaning up
    file.remove(list.files(pattern=".png"))
    }

#    if(!require(rgl)) install.packages("rgl")
#    if(!require(rglwidget)) install.packages("rglwidget")
#    library(rgl)
    
    PCA1 <- PCA.X$x[,1]
    PCA2 <- PCA.X$x[,2]
    PCA3 <- PCA.X$x[,3]
    rgl::plot3d(min(PCA1), min(PCA2), min(PCA3),xlim=c(min(PCA1),max(PCA1)),
     ylim=c(min(PCA2),max(PCA2)),zlim=c(min(PCA3),max(PCA3)),
     xlab="PCA1",ylab="PCA2",zlab="PCA3",
     col = grDevices::rgb(255, 255, 255, 100, maxColorValue=255),radius=0.01)
    num_col=length(unique(PCA.col))
    if(num_col==1)
    { 
      sids1 <- rgl::spheres3d(PCA1, PCA2, PCA3, col = PCA.col,radius=1)
      widgets <- rgl::rglwidget(width = 900, height = 900) %>% rgl::toggleWidget(ids = sids1, label = "PCA")
    }else if(num_col==2)
    {
      index1=PCA.col==unique(PCA.col)[1]
      index2=PCA.col==unique(PCA.col)[2]
      
      sids1 <- rgl::spheres3d(PCA1[index1], PCA2[index1], PCA3[index1], col = PCA.col[index1],radius=1)
      sids2 <- rgl::spheres3d(PCA1[index2], PCA2[index2], PCA3[index2], col = PCA.col[index2],radius=1)
      widgets <- rgl::rglwidget(width = 900, height = 900) %>% rgl::toggleWidget(ids = sids1, label = "Population 1") %>% rgl::toggleWidget(ids = sids2, label = "Population 2")
    }else if(num_col==3)
    {
      index1=PCA.col==unique(PCA.col)[1]
      index2=PCA.col==unique(PCA.col)[2]
      index3=PCA.col==unique(PCA.col)[3]
      
      sids1 <- rgl::spheres3d(PCA1[index1], PCA2[index1], PCA3[index1], col = PCA.col[index1],radius=1)
      sids2 <- rgl::spheres3d(PCA1[index2], PCA2[index2], PCA3[index2], col = PCA.col[index2],radius=1)
      sids3 <- rgl::spheres3d(PCA1[index3], PCA2[index3], PCA3[index3], col = PCA.col[index3],radius=1)
      widgets<-rgl::rglwidget(width = 900, height = 900) %>% rgl::toggleWidget(ids = sids1, label = "Population 1") %>% rgl::toggleWidget(ids = sids2, label = "Population 2") %>% rgl::toggleWidget(ids = sids3, label = "Population 3")
    }else if(num_col==4)
    {
      index1=PCA.col==unique(PCA.col)[1]
      index2=PCA.col==unique(PCA.col)[2]
      index3=PCA.col==unique(PCA.col)[3]
      index4=PCA.col==unique(PCA.col)[4]
      
      sids1 <- rgl::spheres3d(PCA1[index1], PCA2[index1], PCA3[index1], col = PCA.col[index1],radius=1)
      sids2 <- rgl::spheres3d(PCA1[index2], PCA2[index2], PCA3[index2], col = PCA.col[index2],radius=1)
      sids3 <- rgl::spheres3d(PCA1[index3], PCA2[index3], PCA3[index3], col = PCA.col[index3],radius=1)
      sids4 <- rgl::spheres3d(PCA1[index4], PCA2[index4], PCA3[index4], col = PCA.col[index4],radius=1)
      widgets <- rgl::rglwidget(width = 900, height = 900) %>% rgl::toggleWidget(ids = sids1, label = "Population 1") %>% rgl::toggleWidget(ids = sids2, label = "Population 2") %>% rgl::toggleWidget(ids = sids3, label = "Population 3") %>% rgl::toggleWidget(ids = sids4, label = "Population 4")
    }
    if (interactive()) widgets
    htmltools::save_html(widgets, "Interactive.PCA.html")
}
#    if(!require(scatterplot3d)) install.packages("scatterplot3d")
#    library(scatterplot3d)

    grDevices::pdf("GAPIT.PCA.3D.pdf", width = 7, height = 7)
    graphics::par(mar = c(5,5,5,5))
    scatterplot3d::scatterplot3d(PCA.X$x[,1],
                  PCA.X$x[,2],
                  PCA.X$x[,3],
                  xlab = paste("PC",1,sep=""),
                  ylab = paste("PC",2,sep=""),
                  zlab = paste("PC",3,sep=""),
                  pch = 20,
                  color = PCA.col,
                  col.axis = "blue",
                  cex.symbols = 1,
                  cex.lab = 1.4,
                  cex.axis = 1.2,
                  lwd = 3,
                  angle = 55,
                  scale.y = 0.7)
    grDevices::dev.off()
}
print("Joining taxa...")
#Extract number of PCs needed
PCs <- cbind(taxa,as.data.frame(PCA.X$x))

#Remove duplicate (This is taken care by QC)
#PCs.unique <- unique(PCs[,1])
#PCs <-PCs[match(PCs.unique, PCs[,1], nomatch = 0), ]



print("Exporting PCs...")
#Write the PCs into a text file
if(file.output) utils::write.table(PCs[,1:(PCA.total+1)], "GAPIT.PCA.csv", quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)

if(file.output) utils::write.table(PCA.X$rotation[,1:PC.number], "GAPIT.PCA.loadings.csv", quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)

if(file.output) utils::write.table(eigenvalues, "GAPIT.PCA.eigenvalues.csv", quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)

#Return the PCs
return(list(PCs=PCs,EV=PCA.X$sdev^2,nPCs=NULL))
}
#=============================================================================================

