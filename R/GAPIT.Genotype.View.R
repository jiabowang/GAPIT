`GAPIT.Genotype.View` <-function(myGI=NULL,myGD=NULL,chr=NULL, w1_start=NULL,w1_end=NULL,mav1=NULL){
# Object: Analysis for Genotype data:Distribution of SNP density,Accumulation,Moving Average of density,result:a pdf of the scree plot
# myG:Genotype data
# chr: chromosome value
# w1_start:Moving Average windows Start Position
# w1_end:Moving Average windows End Position
# mav1:Moving Average set value length
# Authors: You Tang and Zhiwu Zhang
# Last update: March 11, 2016 
##############################################################################################

#if(nrow(myGI)<1000) return() #Markers are not enough for this analysis
  
if(is.null(myGI)){stop("Validation Invalid. Please select read valid Genotype flies  !")}

if(is.null(myGD)){stop("Validation Invalid. Please select read valid Genotype flies  !")}

if(is.null(w1_start)){w1_start=1}

##if(is.null(w1_end)){w1_end=100}

if(is.null(mav1)){mav1=10}


if(is.null(chr)){chr=1}

#heterozygosity of individuals and SNPs (By Zhiwu Zhang)
  #print("Heterozygosity of individuals and SNPs (By Zhiwu Zhang)")
  X=myGD[,-1]
  H=1-abs(X-1)
  het.ind=apply(H,1,mean)
  het.snp=apply(H,2,mean)
  ylab.ind=paste("Frequency (out of ",length(het.ind)," individuals)",sep="")
  ylab.snp=paste("Frequency (out of ",length(het.snp)," markers)",sep="")
  grDevices::pdf("GAPIT.Heterozygosity.pdf", width =10, height = 6)
  graphics::par(mfrow=c(1,2),mar=c(5,5,1,1)+0.1)
  graphics::hist(het.ind,col="gray", main="",ylab=ylab.ind, xlab="Heterozygosity of individuals")
  graphics::hist(het.snp,col="gray", main="",ylab=ylab.snp, xlab="Heterozygosity of markers")
  grDevices::dev.off()
  rm(X, H, het.ind, het.snp) #Feree memory
  
myFig21<-myGI
myFig21<-myFig21[!is.na(as.numeric(as.matrix(myFig21[,3]))),]

n<-nrow(myFig21)
maxchr<-0
for(i in 1:n){
if(as.numeric(as.matrix(myFig21[i,2]))>maxchr){
maxchr<-as.numeric(as.matrix(myFig21[i,2]))
}
}
n_end<-maxchr
if(maxchr==0){
chr=0
}

#n_end<-as.numeric(as.matrix(myFig21[n,2]))
aaa<-NULL
for(i in 0:n_end){
#myChr<-myFig21[myFig21[,2]==i,]
myChr<-myFig21[as.numeric(as.matrix(myFig21[,2]))==i,]
index<-order(as.numeric(as.matrix(as.data.frame(myChr[,3]))))
aaa<-rbind(aaa,myChr[index,])

}
myFig2<-aaa

if(is.null(w1_end)){
if(nrow(myFig2[as.numeric(as.matrix(myFig2[,2]))==chr,])>100){
w1_end=100
}else{
w1_end=nrow(myFig2[as.numeric(as.matrix(myFig2[,2]))==chr,])
}
}


subResult<-matrix(0,n,1)
for(i in 1 :( n-1))
{
k<-as.numeric(as.matrix(myFig2[i+1,3]))-as.numeric(as.matrix(myFig2[i,3]))
if(k>0){
subResult[i]<-k
}
else{
subResult[i]<-0
}}
results<-cbind(myFig2,subResult)

#####Out  Distribution of SNP density ##########


#####Out Accumulation##########

kk0<-order(as.numeric(as.matrix(results[,4])))

myFig22<-results[kk0,]

m<-nrow(myFig22)

kk1<-matrix(1:m,m,1)
results2<-cbind(myFig22,kk1)
max2<-max(myFig22[,4])


grDevices::pdf("GAPIT.Marker.Density.pdf", width =10, height = 6)
graphics::par(mar=c(5,5,4,5)+0.1)
graphics::hist(as.numeric(as.matrix(results[,4])),
               xlab="Density",
               main="Distribution of SNP",
               breaks=12, cex.axis=0.9,
               col = "dimgray",
               cex.lab=1.3)###,xlim=c(0,25040359))

graphics::par(new=T)
plot(results2[,4],results2[,5]/m,xaxt="n", yaxt="n",bg="lightgray",xlab="",ylab="",type="l",pch=20,col="#990000",cex=1.0,cex.lab=1.3, cex.axis=0.9, lwd=3,las=1,xlim=c(0,max2))
graphics::axis(4,col="#990000",col.ticks="#990000",col.axis="#990000")
graphics::mtext("Accumulation Frequency",side=4,line=3,font=2,font.axis=1.3,col="#990000")
graphics::abline(h=0,col="forestgreen",lty=2)
graphics::abline(h=1,col="forestgreen",lty=2)

grDevices::dev.off()




#####Out Moving Average of density##########
#print(unique(myGI[,2]))

myGD<-myGD[,myGI[,2]==chr]
gc()


myGM0<-myGI[myGI[,2]==chr,]


##remove invalid SNPs
#X<-myGD0[,-1]
X<-myGD
colMax=apply(X,2,max)
colMin=apply(X,2,min)
#mono=as.numeric(colMax)-as.numeric(colMin)
mono=colMax-colMin
index=mono<10E-5
X=X[,!index]


myFig3<-myGM0[!index,]


n3<-nrow(myFig3)


kk3<-order(as.numeric(as.matrix(myFig3[,3])))

myFig23<-myFig3[kk3,]


myGD3<-X[,kk3]

##set windows long 
##w1_start<-30
##w1_end<-230
###get windows numeric snp at the same chr
#print(w1_start)
#print(w1_end)
#print(dim(myFig3))


if(nrow(myFig23)<w1_end)w1_end=nrow(myFig23)

results3_100<-myFig23[w1_start:w1_end,]
myGD3_100<-myGD3[,w1_start:w1_end]

km<-w1_end-w1_start+1
##get number of Density about snp
sum_number_Density <-0

for(j in 1:km)
{
sum_number_Density<-sum_number_Density+(j-1)

}

save_Density_Cor<-matrix(0.0,sum_number_Density,3)
save_Density_Cor_name<-matrix("",sum_number_Density,1)

countSDC<-1
for(j in 1:(km-1))
{
for(k in (j+1):km)
{

save_Density_Cor[countSDC,1]<-abs(as.numeric(as.matrix(results3_100[k,3]))-as.numeric(as.matrix(results3_100[j,3])))
save_Density_Cor[countSDC,2]<-stats::cor(myGD3_100[,j],myGD3_100[,k])
#options(digits=8)
#save_Density_Cor[countSDC,3]<-as.numeric(as.matrix(format(cor(myGD3_100[,j],myGD3_100[,k])%*% cor(myGD3_100[,j],myGD3_100[,k]),digits=8)))
save_Density_Cor[countSDC,3]<-save_Density_Cor[countSDC,2]^2
save_Density_Cor_name[countSDC,1]<-paste(results3_100[j,1],"::::",results3_100[k,1],seq="")
countSDC<-countSDC+1
}
}

#result3_30<-as.data.frame(cbind(save_Density_Cor_name,save_Density_Cor))

k3_3<-order(save_Density_Cor[,1])
result3_3<-save_Density_Cor[k3_3,]

##set moving average value

##mav1<-100

result_mav2<-matrix(0.0,sum_number_Density-mav1,1)

mav1_1<-floor(mav1/2)
mav1_1_end<-sum_number_Density-mav1+mav1_1

result_mav1<-result3_3[(mav1_1+1):mav1_1_end,1]

for(g in 1:(sum_number_Density-mav1)){
sum<-0
for(i in g:(g+mav1-1)){

sum<-sum+result3_3[i,3]

}
#result_mav2[g]<-sum/mav1*5
result_mav2[g]<-sum/mav1
}
result_mav<-cbind(result_mav1,result_mav2)

grDevices::pdf("GAPIT.Marker.LD.pdf", width =10, height = 6)
graphics::par(mar = c(5,5,5,5))

plot(as.matrix(result3_3[,1]),as.matrix(result3_3[,3]),bg="dimgray",xlab="Distance",ylab="R Square",pch=1,cex=0.9,cex.lab=1.2, lwd=0.75,las=1)
#,ylim=c(0,round(max(result3_3[,3]))))

 graphics::lines(result_mav[,2]~result_mav[,1], lwd=6,type="l",pch=20,col="#990000")

grDevices::dev.off()




print(paste("GAPIT.Genotype.View ", ".Two pdf generate.","successfully!" ,sep = ""))

#GAPIT.Genotype.View
}
#=============================================================================================
