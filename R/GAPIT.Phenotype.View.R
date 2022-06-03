`GAPIT.Phenotype.View` <-function(myY=NULL,traitname="_",memo="_"){
# Object: Analysis for Phenotype data:Distribution of density,Accumulation,result:a pdf of the scree plot
# myY:Phenotype data

# Authors: You Tang
# Last update: Sep 7, 2015 
############################################################################################## 
print("GAPIT.Phenotype.View in press...")
if(is.null(myY)){stop("Validation Invalid. Please select read valid Phenotype flies  !")}

# y<-myY[!is.na(myY[,2]),2]
# obs<-as.matrix(y)

# traitname=colnames(myY)[2]

grDevices::pdf(paste("GAPIT",memo,traitname,"phenotype_view.pdf",sep ="."), width =10, height = 6)
# graphics::par(mar = c(5,5,5,5))

# graphics::par(mfrow=c(2,2))
# plot(obs,pch=1)
# #hist(obs)
# graphics::hist(obs,xlab="Density",main="",breaks=12, cex.axis=1,col = "gray")
# graphics::boxplot(obs)
# plot(stats::ecdf(obs),col="red",bg="lightgray",xlab="Density",ylab="Accumulation",main="")
layout.matrix <- matrix(c(1,2,1,3,4,5), nrow = 2, ncol = 3)
layout(mat = layout.matrix,
       heights = c(100,100), # Heights of the two rows
       widths = c(2, 2,2)) # Widths of the two columns
y=myY[!is.na(myY[,2]),2]
par(mar = c(5, 5, 2, 1))
plot(y,xlab="Individual",las=1,ylab="Observation", cex=.5,main="a")
par(mar = c(5, 5, 2, 1))
hist(y,xlab="Observation",las=1,ylab="Frequency",cex=.5,main="c")
par(mar = c(5, 4, 2, 1))
plot(density(na.omit(y)),las=1,xlab="Observation",ylab="Density", cex=.5,main="d")
par(mar = c(5, 4, 2, 1))
boxplot(y,horizontal=F,las=1,xlab="",ylab="Observation", cex=.5,main="b")
par(mar = c(5, 4, 2, 1))
plot(ecdf(y),xlab="Observation",las=1,ylab="Accumulative density", cex=.5,main="e",col="gray40")


grDevices::dev.off()


print(paste("GAPIT.Phenotype.View ", ".output pdf generate.","successfully!" ,sep = ""))

#GAPIT.Phenotype.View
}
#=============================================================================================
