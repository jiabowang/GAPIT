`GAPIT.Validation` <-function(Y=NULL, G=NULL,GD=NULL,GM=NULL,PCA.total=3,KI=NULL,CV=NULL,nfold=NULL,model="gBLUP",file.output=F){
# Object: Genetic Prediction with cross validation 
# nfold:folders number
# Authors: Jiabo Wang and Zhiwu Zhang
# Last update: Mar 15, 2022 
############################################################################################## 
if(is.null(Y)){stop("Validation Invalid. Please input phenotype file !")}
if(ncol(Y)>2) stop("Please just input only one trait to do validation!")
if(is.null(GD)&is.null(G)&is.null(KI))stop ("GAPIT Says:GAPIT need genotype!!!")

all.method=model
File.out=file.output
print("Remove NA individuals in the phenotype file !!!")
Y=Y[!is.na(Y[,2]),]
taxa.y=as.character(Y[,1])
taxa.g=as.character(GD[,1])
Y=Y[taxa.y%in%taxa.g,]

sets=sample(cut(1:nrow(Y),nfold,labels=FALSE),nrow(Y))
colnames(Y)[1]=c("Taxa")
# print(Y)
for(i in 1:length(all.method))
{
	ref_Y_all=NULL
	inf_Y_all=NULL
	for(j in 1:nfold)
    {
        training=Y[,c(1,2)]
        training[sets==j,2]=NA
        training_index=is.na(training[,2])
        # testing=Y[training_index,c(1,2)]
        myBLUP=GAPIT(
	        Y=training,
            GD=GD,
	        GM=GM,
	        KI=KI,
	        PCA.total=PCA.total,
	        model=all.method[i],
	        file.output=FALSE)
        pridiction0=merge(Y,myBLUP$Pred[,c(1,3,5,8)],by.x="Taxa",by.y="Taxa")
        # index=pridiction0[,3]!=2
        ref_Y_all=rbind(ref_Y_all,pridiction0[!training_index,])
        inf_Y_all=rbind(inf_Y_all,pridiction0[training_index,])
        # gblup.r=cor(as.numeric(gapit[index,2]),as.numeric(gapit[index,5]))
    }#end of nfold
    # print(round(min(ref_Y_all[,2]),0))
    # print(round(min(ref_Y_all[,5]),0))
    if(File.out){
    grDevices::pdf(paste("GAPIT.Prediction.", all.method[i],".Ref.pdf", sep = ""), width =6, height = 6)
    graphics::par(mar = c(5,5,5,5))
    plot(ref_Y_all[,2],ref_Y_all[,5],pch=1,
        xlab="Observed(Ref)",ylab="Predicted(Ref)",
        cex.lab=1.3,cex.axis=1.2,lwd=2,main=paste(all.method[i]),
        xlim=c(round(min(ref_Y_all[,2]),0),round(max(ref_Y_all[,2]),0)),
        ylim=c(round(min(ref_Y_all[,5]),0),round(max(ref_Y_all[,5]),0)))   #xlim=c(50,110),ylim=c(50,110),
    kr <- stats::lm(ref_Y_all[,5]~ref_Y_all[,2])
    graphics::abline(a = kr$coefficients[1], b = kr$coefficients[2], col = "red",lwd=4,lty=1)
    graphics::legend("bottomright",paste("R^2=",format(kr$coefficients[2], digits = 4),seq=""), col="white",text.col="blue",lwd=2,cex=1.2,bty="n")
    grDevices::dev.off()

    grDevices::pdf(paste("GAPIT.Prediction.",  all.method[i],".Inf.pdf", sep = ""), width = 6, height = 6)
    graphics::par(mar = c(5,5,5,5))
    plot(inf_Y_all[,2],inf_Y_all[,5],pch=1,
    	xlab="Observed(Inf)",ylab="Predicted(Inf)",
    	cex.lab=1.3,lwd=2,cex.axis=1.2,main=paste(all.method[i]),
    	xlim=c(round(min(inf_Y_all[,2]),0),round(max(inf_Y_all[,2]),0)),
    	ylim=c(round(min(inf_Y_all[,5]),0),round(max(inf_Y_all[,5]),0)))
    ki <- stats::lm(inf_Y_all[,5]~inf_Y_all[,2])
    graphics::abline(a = ki$coefficients[1], b = ki$coefficients[2], col = "red",lwd=3,lty=1)
    graphics::legend("bottomright",paste("R^2=",format(ki$coefficients[2], digits = 4),seq=""), col="white",text.col="blue",lwd=2,cex=1.2,bty="n")
    grDevices::dev.off()

    utils::write.csv(inf_Y_all,paste("GAPIT.Inf.Prediction.",all.method[i],".nfold",nfold,".csv",sep=""),row.names=F)
    utils::write.csv(ref_Y_all,paste("GAPIT.Ref.Prediction.",all.method[i],".nfold",nfold,".csv",sep=""),row.names=F)

    }#end of output
}#end of all.method


# print(inf_Y_all)
# return(list(inf_Y_all,ref_Y_all))
}
#end Prediction one time
#=============================================================================================

