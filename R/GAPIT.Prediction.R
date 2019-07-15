`GAPIT.Prediction` <-function(myK=NULL,y=NULL, num=NULL){
# Object: Genetic Prediction one time by cross validation and cMLM,result:a pdf of the scree plot
# myK:Kinship
# Y: phenotype with columns of taxa,Y1,Y2...
# num:folders number
# Authors: Jiabo Wang and You Tang
# Last update: December 31, 2014 
############################################################################################## 
if(is.null(myK)||is.null(y)){stop("Validation Invalid. Please select read valid flies !")}
if(is.null(num))
  {
	num=5  #not input num value,default folders number is 5
  }

y=y[,1:2]
m=nrow(y)
m.sample=round(m/num)


if(num<2||num>m){stop("Validation Invalid. Please select folders num >1 !")}

vali<-matrix(nr=m.sample,nc=num-1)
cali<-matrix(nr=m-m.sample,nc=num-1)

#vali<-list(v1=unique(as.character(sample(y$Taxa, m.sample))))
#cali<-list(c1=y[!(y$Taxa %in% as.matrix(as.data.frame(vali[1]))), 'Taxa'])

vali[,1]<-unique(as.character(sample(y$Taxa, m.sample)))
cali[,1]<-unique(as.character(y[!(y$Taxa %in% vali[,1]), 'Taxa']))

for(j in 2:num)
{
	if(j!=num)
	{
	 vali[,j]<-unique(as.character(sample(y[!(y$Taxa %in% vali[,1:j-1]), 'Taxa'], m.sample) ))
	}
	if(j==num)
	{
		valilast=unique(as.character(y[!(y$Taxa %in% vali[,1:j-1]), 'Taxa']))
	}

	if(j!=num)
		cali[,j]<-unique(as.character(y[!(y$Taxa %in% vali[,j]), 'Taxa']))
	if(j==num)
		calilast <<- y[!(y$Taxa %in% valilast), 'Taxa']
}

	i=sample(1:num, size = 1)

	if(i!=num){
		lines.vali<-vali[,i]
	  }else{
	 	lines.vali<-valilast
	 }
	 #use only genotypes that were genotyped and phenotyped
	 commonGeno_v <- lines.vali[lines.vali %in% myK[,1]]	               
	 yvali<- y[match(commonGeno_v,y$Taxa),]
    
	 if(i!=num){
		lines.cali<-cali[,i]
	 }else{
		lines.cali<-calilast
	  }
	 #use only genotypes that were genotyped and phenotyped
	 commonGeno_c <- lines.cali[lines.cali %in% myK[,1]]
	 ycali<- y[match(commonGeno_c,y$Taxa),]                
	
	Y.raw=ycali[,c(1,2)]#choos a trait

	myY=Y.raw
	myKI=myK
	max.groups=m
#Run GAPIT
#############################################
	
	blupGAPIT <- GAPIT(
	Y=myY,
	KI=myKI,
	#group.from=max.groups,
	group.from=1,
	group.to=max.groups,
	#group.by=10,
	#PCA.total=3,
	SNP.test=FALSE,
	file.output=FALSE
	)

	blup_prediction=blupGAPIT$GPS
 
	blue<-blupGAPIT$Pred$BLUE
	mean_blue<-mean(blue)

	blup_prediction.ref<-blup_prediction[match(commonGeno_c,blup_prediction$Taxa),]
	blup_prediction.inf<-blup_prediction[match(commonGeno_v,blup_prediction$Taxa),]
	inf_BLUP<-blup_prediction.inf$BLUP
	ref_BLUP<-blup_prediction.ref$BLUP

	inf_pred<-inf_BLUP+mean_blue
	ref_pred<-ref_BLUP+mean_blue


	inf_all<-cbind(blup_prediction.inf,inf_pred)
	ref_all<-cbind(blup_prediction.ref,ref_pred)

	inf_Y_all<-merge(y,inf_all,by.x="Taxa",by.y="Taxa")
	ref_Y_all<-merge(y,ref_all,by.x="Taxa",by.y="Taxa")

	name.of.trait=noquote(names(Y.raw)[2])


pdf(paste("GAPIT.Prediction ", name.of.trait,".Predict reference.pdf", sep = ""), width =6, height = 6)
par(mar = c(5,5,5,5))
plot(ref_Y_all[,2],ref_Y_all[,8],pch=1,xlab="Observed(Ref)",ylab="Predicted(Ref)",cex.lab=1.3,cex.axis=1.2,lwd=2)   #xlim=c(50,110),ylim=c(50,110),
kr<-lm(ref_Y_all[,8]~ref_Y_all[,2])
abline(a = kr$coefficients[1], b = kr$coefficients[2], col = "red",lwd=4,lty=1)
#v1<-max(ref_Y_all[,2]])*10/10
#text(v1,kr$coefficients[1]+kr$coefficients[2]*v1,paste("R^2=",format(kr$coefficients[2], digits = 3),seq=""), col = "blue", adj = c(0, -.1))
legend("bottomright",paste("R^2=",format(kr$coefficients[2], digits = 4),seq=""), col="white",text.col="blue",lwd=2,cex=1.2,bty="n")

dev.off()
pdf(paste("GAPIT.Prediction ", name.of.trait,".Predict inference.pdf", sep = ""), width = 6, height = 6)
par(mar = c(5,5,5,5))
plot(inf_Y_all[,2],inf_Y_all[,8],pch=1,xlab="Observed(Inf)",ylab="Predicted(Inf)",cex.lab=1.5,lwd=2,,cex.axis=1.2)#xlim=c(50,110),ylim=c(45,100),
ki<-lm(inf_Y_all[,8]~inf_Y_all[,2])
abline(a = ki$coefficients[1], b = ki$coefficients[2], col = "red",lwd=3,lty=1)
#v0<-max(inf_Y_all[,2])
#text(v0,ki$coefficients[1]+ki$coefficients[2]*v0,paste("R^2=",format(ki$coefficients[2], digits = 4),seq=""), col = "blue", adj = c(0, -.1))

legend("bottomright",paste("R^2=",format(ki$coefficients[2], digits = 4),seq=""), col="white",text.col="blue",lwd=2,cex=1.2,bty="n")

dev.off()
print(paste("GAPIT.Prediction ", name.of.trait,".Predict phenotype.","successfully!" ,sep = ""))
return(list(inf_Y_all,ref_Y_all))
}
#end Prediction one time
#=============================================================================================

