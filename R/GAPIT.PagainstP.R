`GAPIT.PagainstP`<-
function(container,Y,testY,model_store,traitname0="",lmpred=NULL,type=c("GEBV"),pch0=NULL,color0=NULL,
  Cross.Vali=TRUE,byTraits=TRUE,
  memo=NULL
  )
#model_store is the store of all model names
#Y is the real phenotype
#type could be set as BLUP, BLUE, or Pred
#Cross.Vali indicate whether show cross validation  (default as TRUE)
#
{ 
# Y=myY
# testY=test
# model_store=c("gBLUP","cBLUP","sBLUP")
# model_store=c("BLINK","FarmCPU")
Y.names=colnames(Y)[2]
print("GAPIT.PagainstP has been beging...")
model_store2=model_store
if(length(model_store)==length(model_store[model_store%in%c("gBLUP","cBLUP","sBLUP")]))
  {
    container=paste(model_store2,".",traitname0,sep="")
  }else{
    if(length(model_store)==length(model_store[model_store%in%c("MLM","GLM","CMLM","SUPER","MLMM","MLMM2","FarmCPU","FarmCPU2","BLINK","BLINK2","BLINKC")]))
      {
        model.n=length(model_store)
        lmpred.n=length(lmpred)
        model2=NULL
        for(i in 1:model.n)
          {
            model2=append(model2,rep(model_store[i],2))
          }
        lmpred2=rep(lmpred,model.n)
        pred.way=rep(".MAS",length(lmpred2))
        pred.way[!lmpred2]=".ABLUP"
        container=paste(model2,".",traitname0,pred.way,sep="")
      }else{
        blup.index=model_store%in%c("gBLUP","cBLUP","sBLUP")
        model.n=length(model_store)
        lmpred.n=length(lmpred)
        model2=NULL
        for(i in 1:model.n)
          {
            dupl=2
            if(model_store[i]%in%c("gBLUP","cBLUP","sBLUP"))dupl=1
            model2=append(model2,rep(model_store[i],dupl))
          }
        if(length(lmpred)==1)
          {
            lmpred2=rep(lmpred,model.n)
            pred.way=rep(".MAS",length(lmpred2))
            pred.way[!lmpred2]=".ABLUP"
            container=paste(model2,".",traitname0,pred.way,sep="")
          }else{
            container=NULL
            cm=1
            pred.way=c(".MAS",".ABLUP")
            for(i in 1:length(model2))
              {
                if(model2[i]%in%c("gBLUP","cBLUP","sBLUP"))
                  {
                    # model2.tem=ifelse(model2[i]=="gBLUP","MLM",ifelse(model2[i]=="cBLUP","CMLM","SUPER"))
                    container=append(container,paste(model2[i],".",traitname0,sep=""))
                  }else{
                    container=append(container,paste(model2[i],".",traitname0,pred.way[1+cm%%2],sep=""))
                    cm=cm+1
                  }
              }
          }# end of length(lmpred)==1 else
      }# end of gwas model
  }# end of all model
        
environ_name=container
n=length(container)
# print(n)
# print(environ_name)
# method_store=NULL
obser=Y
colnames(obser)=c("taxa","observed")
colnames(testY)=c("taxa","observed")
cv.index=c(rep(FALSE,nrow(obser)),rep(TRUE,nrow(testY)))
# print(dim(testY))
# print(table(cv.index))
obser=rbind(obser,testY)
cv.index=cv.index[!is.na(obser[,2])]
obser=obser[!is.na(obser[,2]),]
if(Cross.Vali) 
{
  obser=obser[cv.index,]
  cv.index=cv.index[cv.index]
}
print(environ_name)
gs.index=ifelse(type=="BLUP",5,ifelse(type=="BLUE",7,8))
gs_store=obser
for(i in 1:n)
   {
    gs_result=utils::read.csv(paste("GAPIT.Association.Prediction_results.",environ_name[i],".csv",sep=""),head=T)
    m=nrow(gs_result)
    gs_result0=gs_result[,c(1,gs.index)]
    colnames(gs_result0)=c("Taxa",paste(environ_name[i],sep=""))
    # print(head(gs_result0))
    gs_store=merge(gs_store,gs_result0,by.x=colnames(obser)[1],by.y=colnames(gs_result0)[1])
   }
# colnames(gs_store)[-1]=model_store
print(dim(gs_store))
x.max=ceiling(max(gs_store[,2]))
x.min=floor(min(gs_store[,2]))
y.max=ceiling(max(gs_store[,-c(1,2)]))
y.min=floor(min(gs_store[,-c(1,2)]))
if(is.null(pch0))pch0=c(21:25)[1:n]
# if(is.null(color0))color0=rainbow(7)[1:n]
if(is.null(color0))color0=c("turquoise4","indianred3","darkolivegreen3","red","aquamarine3","darkgoldenrod")[1:n]
# if(is.null(color0))color0=c("lightblue","mistyrose","lavender")[1:n]
grDevices::pdf(paste("GAPIT.Association.Prediction_",type,".pdf" ,sep = ""),width = 8,height=5)

par(mfrow=c(1,1))
par(mar=c(5,7,1,1))
hx=seq(x.min,x.max,abs(x.max-x.min)/5)
hy=seq(y.min,y.max,abs(y.max-y.min)/5)

plot(gs_store[1,2],gs_store[1,3],xlab="",ylab="",
      xlim=c(x.min,x.max+0.2*x.max),ylim=c(y.min,y.max),
      las=1,axes=F,
      pch=1,col="white",cex=1,lwd=1)
     abline(h=hy,col="gray")
     abline(v=hx,col="gray")
r.store=NULL
for(i in 1:n)
   {
     par(new=T)
     color1=rep(color0[i],nrow(gs_store))
     color1[cv.index]="white"
     plot(gs_store[,2],gs_store[,i+2],xlab="",ylab="",
      xlim=c(x.min,x.max+0.2*x.max),ylim=c(y.min,y.max),
      las=1,axes=F,
      pch=pch0[i],col=color0[i],cex=1,lwd=1,bg=color1)
     r.store=append(r.store,cor(gs_store[,2],gs_store[,i+2]))
     # abline(h=hy,col="gray")
     # abline(v=hx,col="gray")
   }
axis(1,col="black",col.ticks="black",col.axis="black",tck=-0.02,xaxp=c(floor(x.min),ceiling(x.max),5),cex.axis=1)
axis(2,col="black",col.ticks="black",tck=-0.01,col.axis="black",yaxp=c(floor(y.min),ceiling(y.max),5),las=1,cex.axis=1)
mtext("Observed Phenotype",side=1,line=2.6,col="black",cex=1)
mtext(paste("Predicted ",type,sep="" ),side=2,line=3.5,col="black",cex=1)

legend("bottomright",legend=paste("R (",colnames(gs_store)[-c(1,2)],")= ",round(r.store,2),sep=""),horiz=F,
          col=color0,pch=pch0,lwd=1,cex=0.7,lty=0,ncol=1,
          bty = "o", bg = "white")
grDevices::dev.off()
print("GAPIT.PagainstP Figures have been done!!!")

}# end of function

