`GAPIT.cross_validation.compare` <-function(myGD=NULL,y=NULL, rel=NULL,tc=NULL){
# Object: GAPIT.cross validation compare to different folders by replicate Times,result:a pdf of the scree barplot and .cvs
# myGD:numeric SNP
# Y: phenotype with columns of taxa,Y1,Y2...
# rel:replications
# tc:comparation folds number and value
# Authors: You Tang,Jiabo Wang and You Zhou
# Last update: December 31, 2014 
##############################################################################################
if(is.null(myGD)||is.null(y)){stop("Validation Invalid. Please select read valid flies !")}
if(is.null(rel))
  {
	rel=10  #not input rel value,default replications number is 10
  }

if(rel<2){stop("Validation Invalid. Please select replications >1 !")}
#rel<-2 ##replications
#t<-2
y<-y[!is.na(y[,2]),] 
y<-y[,c(1,2)]
y<- stats::na.omit(y)
#############
commonGeno <- unique(as.character(y[,1]))[unique(as.character(y[,1])) %in% myGD[,1]]
cG<-data.frame(commonGeno)
names(cG)<-"Taxa"
colnames(y)<-c("Taxa","pheno")
y2<-merge(y,cG,all.x=FALSE, all.y=TRUE, by = c("Taxa"))
Z1 <- myGD[match(y2$Taxa,myGD[,1]),]
myGD<- Z1
y<-y2
##############
X<-myGD[,-1]
k1<-as.matrix(X)
k2=GAPIT.kinship.VanRaden(snps=k1)
myKI<-as.data.frame(k2)
myKI<-cbind(myGD[,1],myKI)
utils::write.table(y,"Y.txt",quote=F,sep="\t",row.names=F,col.names=T)
utils::write.table(myKI,"K.txt",quote=F,row.names=F,col.names=F,sep="\t")
gc()
myK<- utils::read.table("K.txt",head=F)
y = utils::read.table("Y.txt",head=T)

y <- stats::na.omit(y)
y=y[(y[,1] %in% myK[,1]),]
m=nrow(y)
if(is.null(tc))
	tc<-c(2,5,10,20,50)  ##default compare to folders num
tc1<-as.matrix(tc)
	allstorage.ref=matrix(0,rel,nrow(tc1))
	allstorage.inf=matrix(0,rel,nrow(tc1))
for(w in 1:nrow(tc1)){
num<-tc1[w,]
m.sample=floor(m/num)
	storage.ref=matrix(0,rel,num)
	storage.inf=matrix(0,rel,num)
	#storage.REML=matrix(0,rel,num)
for(k in 1:rel)
{
   #################Rand group method 1############
 sets=sample(cut(1:nrow(y),num,labels=FALSE),nrow(y))
 sets = as.data.frame(sets)
 ynew <- cbind(sets,y)

	#i=sample(1:num, size = 1)
for(i in 1:num){
	
	 #use only genotypes that were genotyped and phenotyped
    ref <- y$Taxa[!ynew$sets==i]
      
     lines.cali<- ref     
   # ycali<- y[match(ref,y$Taxa),]
    #use only genotypes that were genotyped and phenotyped

    test <- y$Taxa[ynew$sets==i]
    lines.vali<-test 
    #yvali<- y[match(test,y$Taxa),]  	
 
 #################end Rand group method############

	 #use only genotypes that were genotyped and phenotyped
	 commonGeno_v <- lines.vali[lines.vali %in% myK[,1]]	               
	 yvali<- y[match(commonGeno_v,y[,1]),]    

	 #use only genotypes that were genotyped and phenotyped
	 commonGeno_c <- lines.cali[lines.cali %in% myK[,1]]
	 ycali<- y[match(commonGeno_c,y[,1]),]               
	
	Y.raw=ycali[,c(1,2)]#choos a trait

	myY=Y.raw
	myKI=myK
	max.groups=m

#Run GAPIT
#############################################
	
	myGAPIT <- GAPIT(
	Y=myY,
	KI=myKI,
	# #group.from=max.groups,
	# group.from=max.groups,
	# group.to=max.groups,
	model="gBLUP",
	#group.by=10,
	# PCA.total=3,
	SNP.test=FALSE,
	file.output=FALSE
	)

prediction=myGAPIT$Pred

prediction.ref<-prediction[match(commonGeno_c,prediction$Taxa),]
prediction.inf<-prediction[match(commonGeno_v,prediction$Taxa),]

YP.ref <- merge(y, prediction.ref, by.x = 1, by.y = "Taxa")
YP.inf <- merge(y, prediction.inf, by.x = 1, by.y = "Taxa")

#Calculate correlation and store them
r.ref=stats::cor(as.numeric(as.vector(YP.ref[,2])),as.numeric(as.vector(YP.ref[,6]) ))
r.inf=stats::cor(as.numeric(as.vector(YP.inf[,2])),as.numeric(as.vector(YP.inf[,6]) ))

if(r.inf<0){
#r.inf=cor(as.numeric(as.vector(YP.inf[,2])),as.numeric(as.vector(YP.inf[,2]+YP.inf[,6])))
combine_output=cbind(as.numeric(as.vector(YP.inf[,2])),as.numeric(as.vector(YP.inf[,6]) ))

utils::write.csv(combine_output, paste("Accuracy_folders",num,k,i,rel,".csv",sep=""))
#stop("...........")
}
storage.ref[k,i]=r.ref
storage.inf[k,i]=r.inf

print(paste(" rel= ", rel, " k= ",k," i= ",i,sep = ""))
}
print(paste("finish  replications k= ",k," folders= ",num,sep = ""))
}
#Find missing position-->0.0
index=is.na(storage.inf)
storage.inf[index]=0


allstorage.inf[,w]=as.matrix(rowMeans(storage.inf))
allstorage.ref[,w]=as.matrix(rowMeans(storage.ref))
#as.matrix(rowMeans(storage.ref))

##output rel times and accuracy for every folders 

combine_output=cbind(storage.inf,allstorage.inf[,w])
combine_output1=cbind(storage.ref,allstorage.ref[,w])
colnames(combine_output)=c(paste("folders",c(1:num),sep=""),"mean")
utils::write.csv(combine_output, paste("Accuracy_folders",num,"by CMLM,rel_",rel,".csv",sep=""))
utils::write.csv(combine_output1, paste("Accuracy_folders  ref",num,"by CMLM,rel_",rel,".csv",sep=""))

}	
sr<-nrow(tc1)
##output means accuracy by rel for every folders 
colnames(allstorage.inf)=c(paste(tc1[c(1:sr),]," folders",sep=""))
utils::write.csv(allstorage.inf, paste("Accuracy_folders",nrow(tc1),"by CMLM,rel_",rel,".compare to means",".csv",sep=""))
utils::write.csv(allstorage.ref, paste("Accuracy_folders  ref",nrow(tc1),"by CMLM,rel_",rel,".compare to means",".csv",sep=""))

	name.of.trait=noquote(names(Y.raw)[2])
#rrel=round(rel/2)
#ppp<-matrix(0,sr,2)
ppp<-matrix(0,sr,2)

#if(rrel!=1){
#	aarm<-colMeans(allstorage.inf[1:rrel,])
#	}else{
#	aarm<-allstorage.inf[1,]	
#	}
#aam<-colMeans(allstorage.inf)
aam<-allstorage.inf
aam<-data.frame(aam)
bbm<-allstorage.ref
bbm<-data.frame(bbm)
for(b in 1:sr){
#ppp[b,]<-as.matrix(c(aarm[b],aam[b]))
ppp[b,1]<-as.matrix(mean(aam[,b]))
#colnames(ppp)<-c(rrel,rel)
}
for(c in 1:sr){
ppp[c,2]<-as.matrix(mean(bbm[,c]))
}
ppp<-as.matrix(cbind(ppp,tc1))
#colnames(ppp)<-c(rel)
sj <- stats::runif(1, 0, 1)
#name.of.trait="qqq"
grDevices::pdf(paste("GAPIT.cross_validation ", name.of.trait,sj,".compare to different folders.", ".pdf", sep = ""),width = 4.5, height = 4,pointsize=9)
graphics::par(mar = c(5,6,5,3))
grDevices::palette(c("blue","red",grDevices::rainbow(2)))
plot(ppp[,3],ppp[,2],xaxt="n",ylim=c(0,1.04),xlim=c(min(tc1)-1,max(tc1)+1),bg="lightgray",xlab="Number of folds",ylab="Correlation",type="o",pch=1,col=1,cex=1.0,cex.lab=1.7, cex.axis=1.3, lwd=3,las=1,lty =2)
	graphics::axis(side=1,at=tc1,labels=tc1,cex.lab=1.7)
        graphics::lines(ppp[,1]~ppp[,3], lwd=3,type="o",pch=19,col=2,lty =1)
	graphics::legend("bottomright",horiz = FALSE,c("Reference","Inference"),pch = c(1,19), lty =c(2,1),col=c(1:2),lwd=2,cex=1.2,bty="n")
grDevices::dev.off()
print(paste("GAPIT.cross validation ", name.of.trait,".compare to different folders.","successfully!" ,sep = ""))
return(list(allstorage.inf))
}#end GAPIT.cross validation compare to different folders by replicate Times
#=============================================================================================
