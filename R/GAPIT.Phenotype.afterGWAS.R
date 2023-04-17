`GAPIT.Phenotype.afterGWAS` <-
function(GWAS, GD, G=NULL,GM,Y,model=NULL,cutOff=0.01,seed=1234){
    #Object: Make phenotype distribution based on the significant markers from GWAS
    #Options for 
    #
    #Output: A pdf of henotype distribution boxplot
    #Authors: Jiabo Wang and Zhiwu Zhang
    # Last update: Mar 13, 2023
    ##############################################################################################
print("GAPIT Phenotype distribution with significant markers in process...")
hapmap=FALSE
if(!is.null(G)) hapmap=TRUE
print(hapmap)
# print(dim(G))
set.seed(seed)
# GWAS=myGAPIT$GWAS
# cutOff=0.01
# Y=mysimulation$Y
colnames(Y)[1]=c("Taxa")
# GD=myGD
# GM=myGM

bonferroniCutOff=cutOff/nrow(GWAS)
sigs=GWAS[GWAS[,4]<bonferroniCutOff,,drop=FALSE]
N.sigs=nrow(sigs)
print(dim(sigs))
letter=letters[1:N.sigs]
if(N.sigs<1|N.sigs>15) return()
if(N.sigs<4)
{
    x.layout=N.sigs
    y.layout=1
}else{
    prime=FALSE
    Pi=0
    for(i in 2:(N.sigs-1))
    {
        if((N.sigs%%i)==0)
        {
            prime=TRUE
            break
        }
    }
    x.layout=ifelse(N.sigs==14,5,ifelse(prime,N.sigs/i,ifelse(N.sigs<6,3,ifelse(N.sigs<12,4,5))))
    y.layout=ifelse(N.sigs==14,3,ifelse(prime,i,ifelse(N.sigs<10,2,3)))
}
sigs=sigs[order(sigs[,4]),]
trait.name=colnames(Y)[2]
# x.layout
# y.layout
grDevices::pdf(paste("GAPIT.Phenotype.Distribution_Significantmarkers.",model,".",trait.name,".pdf",sep=""), width =3*x.layout, height = 3*y.layout)

par(mfrow=c(y.layout,x.layout),mar = c(5,5,2,2))

# y=Y[,2]
y.min=round(min(Y[,2],rm.na=T),1)
y.max=round(max(Y[,2],rm.na=T),1)
if(hapmap)
{
  X=t(G[-1,-c(1:11)])
  taxa=as.character(G[1,-c(1:11)])
  map=G[-1,c(1,3,4)]
}else{
  X=GD[,-1]
  taxa=as.character(GD[,1])
  map=GM
}

# print(letter)
# print(dim(X))
for(i in 1:N.sigs)
{
  marker=sigs[i,,drop=FALSE]
  marker.name=as.character(marker[,1])
  marker.index=map[,1]%in%marker.name
  marker.genotype=cbind(taxa,as.data.frame(as.character(X[,marker.index])))
  type.num=length(unique(marker.genotype[,2]))
  marker.type=as.character(unique(marker.genotype[,2]))
  marker.type=marker.type[order(marker.type)]
  yall=merge(Y,marker.genotype,by.x="Taxa",by.y="taxa")
  colnames(yall)=c("Taxa","Values","Genotype")
  marker.taxa=paste(marker.name,":",marker[,2],":",marker[,3],sep="")
  boxplot(Values~Genotype,data=yall,xlab="",ylab="",
    las=1,ylim=c(y.min,y.max),main=letter[i],
    space=0.2,axes=F,outline=FALSE)
  for(j in 1:type.num)
  {
    yj=yall[yall[,3]==marker.type[j],2]
    points((j+runif(length(yj),min=-0.2,max=0.2) ), yj, cex=0.7,pch = 1,  col="blue")

  }# end of j
  axis(2,col="black",col.ticks="black",col.axis="black",tck=-0.02,las=1,cex.axis=1.5)
  axis(1,at=1:type.num,labels=marker.type,col="black",col.ticks="black",col.axis="black",tck=-0.01,tick=F,cex.axis=1.5)
    # axis(1,at=posi,labels=labels,col="black",col.ticks="black",col.axis="black",tck=-0.01,tick=F)
  # axis(3,at=posi[even],labels=labels[even],col="black",col.ticks="black",col.axis="black",tck=-0.01,tick=F)
  mtext(paste(marker.taxa,sep=""),side=1,line=3.4,cex=1.2)
  mtext(trait.name,side=2,line=3,cex=1.2)
  # legend("top", letters[i], col=c("red","black","blue"),xpd=NA,text.col = "black", pch = c(19,19,19), merge = T, bg = "white",ncol=1, cex = 1.5, lwd=-2, bty='n')

}# end of i

grDevices::dev.off()

} #end of function
#=============================================================================================

