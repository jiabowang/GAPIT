`GAPIT.Phenotype.PCA.View` <-function(PC=NULL,myY=NULL){
# Object: Analysis PCA effection for Phenotype data ,result:a pdf of the scree plot
# myG:Genotype data
# myY:Phenotype data

# Authors: You Tang
# Last update: Sep 7, 2015 
############################################################################################## 
print("GAPIT.Phenotype.PCA.View")
if(is.null(PC)){stop("Validation Invalid. Please input four PC value  !")}
if(is.null(myY)){stop("Validation Invalid. Please select read valid Phenotype flies  !")}

y<-myY[!is.na(myY[,2]),c(1:2)]

traitname=colnames(y)[2]

cv1<-PC[!is.na(match(PC[,1],y[,1])),]
y1<-y[!is.na(match(y[,1],cv1[,1])),]

y2<-y1[order(y1[,1]),]
cv2<-cv1[order(cv1[,1]),]
lcor=round(cor(y2[,-1],cv2[,-1])*100)/100

y.range=max(y2[,2])-min(y2[,2])
y.mean=mean(y2[,2])
n.col=54
y.int=round(abs(y2[,2]-y.mean)/y.range*(.5*n.col-1)*2)+1
mycol=rainbow(n.col)
y.col=mycol[y.int]
y.lab=paste("PC",seq(1:4)," (r=",lcor,")",sep="")

pdf(paste("GAPIT.",traitname,"_vs_PC.pdf",sep=""), width =9, height = 6)
#par(mar = c(5,5,5,5))
par(mar = c(5,5,2,2))
par(mfrow=c(2,2))

plot(y2[,2],cv2[,2],bg="lightgray",xlab="Phenotype",ylab=y.lab[1],main="",cex.lab=1.4,col=y.col)
if(ncol(PC)>2) plot(y2[,2],cv2[,3],bg="lightgray",xlab="Phenotype",ylab=y.lab[2],main="",cex.lab=1.4,col=y.col)
if(ncol(PC)>3) plot(y2[,2],cv2[,4],bg="lightgray",xlab="Phenotype",ylab=y.lab[3],main="",cex.lab=1.4,col=y.col)
if(ncol(PC)>4) plot(y2[,2],cv2[,5],bg="lightgray",xlab="Phenotype",ylab=y.lab[4],main="",cex.lab=1.4,col=y.col)

dev.off()


print(paste("GAPIT.Phenotype.PCA.View ", ".output pdf generate.","successfully!" ,sep = ""))

#GAPIT.Phenotype.View
}
#=============================================================================================
