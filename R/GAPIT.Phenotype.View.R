`GAPIT.Phenotype.View` <-function(myY=NULL,traitname="_",memo="_"){
# Object: Analysis for Phenotype data:Distribution of density,Accumulation,result:a pdf of the scree plot
# myY:Phenotype data

# Authors: You Tang
# Last update: Sep 7, 2015 
############################################################################################## 
print("GAPIT.Phenotype.View in press...")
if(is.null(myY)){stop("Validation Invalid. Please select read valid Phenotype flies  !")}

y<-myY[!is.na(myY[,2]),2]
obs<-as.matrix(y)

traitname=colnames(myY)[2]

grDevices::pdf(paste("GAPIT",memo,traitname,"phenotype_view.pdf",sep ="."), width =10, height = 6)
graphics::par(mar = c(5,5,5,5))

graphics::par(mfrow=c(2,2))
plot(obs,pch=1)
#hist(obs)
graphics::hist(obs,xlab="Density",main="",breaks=12, cex.axis=1,col = "gray")
graphics::boxplot(obs)
plot(stats::ecdf(obs),col="red",bg="lightgray",xlab="Density",ylab="Accumulation",main="")

grDevices::dev.off()


print(paste("GAPIT.Phenotype.View ", ".output pdf generate.","successfully!" ,sep = ""))

#GAPIT.Phenotype.View
}
#=============================================================================================
