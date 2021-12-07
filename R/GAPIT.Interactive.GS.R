`GAPIT.Interactive.GS`<-
function(model_store=NULL,Y=NULL,type=c("Pred"),testY=NULL
  )
#model_store is the store of all model names
#Y is the real phenotype
#
{ 
# Y=myY
# Y=training
# model_store=c("gBLUP","cBLUP","sBLUP")


n=length(model_store)
method_store=NULL
obser=Y
colnames(obser)=c("taxa","observed")

if("gBLUP"%in%model_store)method_store=append(method_store,"MLM")
if("cBLUP"%in%model_store)method_store=append(method_store,"CMLM")
if("sBLUP"%in%model_store)method_store=append(method_store,"SUPER")
index=c("gBLUP","cBLUP","sBLUP")%in%model_store
no_model=c("gBLUP","cBLUP","sBLUP")[!index]
gs_store=NULL
for(i in 1:n)
   {
    gs_result=utils::read.csv(paste("GAPIT.",method_store[i],".Pred.result.csv",sep=""),head=T)
    m=nrow(gs_result)
    gs_store=cbind(gs_store,gs_result[,8])
   }
colnames(gs_store)=model_store
taxa=as.character(gs_result[,1])
refinf=as.numeric(gs_result[,3])
refinf[refinf>1]=4
pred=cbind(as.data.frame(taxa),gs_store,refinf)
colnames(pred)[1]="taxa"
if(!is.null(testY))
  {
    testY=testY[,c(1,2)]
    colnames(testY)=c("taxa","observed")
    obser2=obser[!is.na(obser[,2]),]
    obser=rbind(obser2,testY)
  }


pred_all=merge(pred,obser,by.x="taxa",by.y="taxa")

taxa=as.character(pred_all[,1])

if(!setequal(no_model,character(0))) 
  {
    nn=length(no_model)
    for(j in 1:nn)
       {
        one=rep(NA,nrow(pred_all))
        pred_all=cbind(pred_all,one)
        colnames(pred_all)[ncol(pred_all)]=no_model[j]
       }
  }
# e=20
# #NQTN=100
# #h2=0.25
# taxa=as.character(myGD[,1])
# myY=Y0[Y0[,1]%in%taxa,c(1,e)]
# myGD=myGD[taxa%in%myY[,1],]
# nfold=5
# repli=1
# sets=sample(cut(1:nrow(myY ),nfold,labels=FALSE),nrow(myY ))


# j=1
# training=myY
# training[sets==j,2]=NA
# training_index=is.na(training[,2])
# testing=myY[training_index,]

# cblup_gapit=GAPIT(Y=training,CV=PC,PCA.total=0,KI=myKI,group.from=200,group.to=2000,group.by=600,SNP.test=F,file.output=F)
# gblup_gapit=GAPIT(Y=training,CV=PC,PCA.total=0,KI=myKI,group.from=2000,group.to=2000,group.by=100,SNP.test=F,file.output=F)
# sblup_gapit=GAPIT(Y=training,CV=PC,PCA.total=0,GD=myGD,GM=myGM,group.from=2000,SUPER_GS=TRUE,sangwich.top="MLM",sangwich.bottom="SUPER",LD=0.1,SNP.test=F,file.output=F,inclosure.from=200,inclosure.to=1000,inclosure.by=200,bin.from=10000,bin.to=100000,bin.by=10000)

# cblup_pred=cblup_gapit$Pred[training_index,]
# gblup_pred=gblup_gapit$Pred[training_index,]
# sblup_pred=sblup_gapit$Pred[training_index,]
# testing_index=!is.na(testing[,2])

# gblup_r_once=cor(testing[testing_index,2],gblup_pred[testing_index,8])
# cblup_r_once=cor(testing[testing_index,2],cblup_pred[testing_index,8])
# sblup_r_once=cor(testing[testing_index,2],sblup_pred[testing_index,8])
# result=cbind(testing[testing_index,],gblup_pred[testing_index,8],cblup_pred[testing_index,8],sblup_pred[testing_index,8])
# colnames(result)=c("taxa","observed","gBLUP","cBLUP","sBLUP")

# gblup_r_once
# cblup_r_once
# sblup_r_once

# write.table(result,paste("gcs_",e,".txt",sep=""))




#myY=read.table(paste("gcs_",e,".txt",sep=""),head=T)
Observed=pred_all$observed[pred_all$refinf==1]
Predicted=pred_all$gBLUP[pred_all$refinf==1]

#if(!require(plotly)) install.packages("plotly")
#  library(plotly)

  # p <- plot_ly(
  #   type = 'scatter',
  #   x = ~Observed,
  #   y = ~Predicted,
  #   data=pred_all,
  #   text = ~paste("Taxa: ",taxa,"<br>Observed: ",round(observed,4) , paste("'<br>",colnames(pred_all)[2],":'",sep=""), round(gBLUP,4)),
  #   #size=2*y/max(y),
  #   color = I("red"),
  #   symbol= I(refinf),
  #   name=colnames(pred_all)[2]
  #   )%>%add_trace(
  #   type = 'scatter',
  #   x = ~observed,
  #   y = ~cBLUP,
  #   #data=myY,
  #   text = ~paste("Taxa: ",taxa,"<br>Observed: ",round(observed,4)  , '<br>cBLUP:', round(cBLUP,4)),
  #   #size=2*y/max(y),
  #   color = I("blue"),
  #   symbol= I(refinf),
  #   name=c("cBLUP")
  #   )%>%add_trace(
  #   type = 'scatter',
  #   x = ~observed,
  #   y = ~sBLUP,
  #   #data=myY,
  #   text = ~paste("Taxa: ",taxa,"<br>Observed: ",round(observed,4)  , '<br>sBLUP:', round(sBLUP,4)),
  #   #size=2*y/max(y),
  #   color = I("green"),
  #   symbol= I(refinf),
  #   name=c("sBLUP")
  #   )
  #   htmltools::save_html(p, "Interactive.GS.html")

#####



 p <- plotly::plot_ly(
    type = 'scatter',
    x = ~Observed,
    y = ~Predicted,
    data=pred_all,
    text = ~paste("Taxa: ",taxa[pred_all$refinf==1],"<br>Observed: ",round(Observed,4) , '<br>gBLUP:', round(Predicted,4)),
    #size=2*y/max(y),
    color = I("red"),
    symbol= I(1),
    name=c("gBLUP with Ref")
    ) %>% plotly::add_trace(
    type = 'scatter',
    x = ~observed[pred_all$refinf==1],
    y = ~cBLUP[pred_all$refinf==1],
    #data=myY,
    text = ~paste("Taxa: ",taxa[pred_all$refinf==1],"<br>Observed: ",round(observed[pred_all$refinf==1],4)  , '<br>cBLUP:', round(cBLUP[pred_all$refinf==1],4)),
    #size=2*y/max(y),
    color = I("blue"),
    symbol= I(1),
    name=c("cBLUP with Ref")
    ) %>% plotly::add_trace(
    type = 'scatter',
    x = ~observed[pred_all$refinf==1],
    y = ~sBLUP[pred_all$refinf==1],
    #data=myY,
    text = ~paste("Taxa: ",taxa[pred_all$refinf==1],"<br>Observed: ",round(observed[pred_all$refinf==1],4)  , '<br>sBLUP:', round(sBLUP[pred_all$refinf==1],4)),
    #size=2*y/max(y),
    color = I("green"),
    symbol= I(1),
    name=c("sBLUP with Ref")
    ) %>% plotly::add_trace(
    type = 'scatter',
    x = ~observed[pred_all$refinf>1],
    y = ~cBLUP[pred_all$refinf>1],
    #data=myY,
    text = ~paste("Taxa: ",taxa[pred_all$refinf>1],"<br>Observed: ",round(observed[pred_all$refinf>1],4)  , '<br>cBLUP:', round(cBLUP[pred_all$refinf>1],4)),
    #size=2*y/max(y),
    color = I("blue"),
    symbol= I(4),
    name=c("cBLUP with Inf")
    ) %>% plotly::add_trace(
    type = 'scatter',
    x = ~observed[pred_all$refinf>1],
    y = ~sBLUP[pred_all$refinf>1],
    #data=myY,
    text = ~paste("Taxa: ",taxa[pred_all$refinf>1],"<br>Observed: ",round(observed[pred_all$refinf>1],4)  , '<br>sBLUP:', round(sBLUP[pred_all$refinf>1],4)),
    #size=2*y/max(y),
    color = I("green"),
    symbol= I(4),
    name=c("sBLUP with Inf")
    ) %>% plotly::add_trace(
    type = 'scatter',
    x = ~observed[pred_all$refinf>1],
    y = ~gBLUP[pred_all$refinf>1],
    #data=myY,
    text = ~paste("Taxa: ",taxa[pred_all$refinf>1],"<br>Observed: ",round(observed[pred_all$refinf>1],4)  , '<br>gBLUP:', round(gBLUP[pred_all$refinf>1],4)),
    #size=2*y/max(y),
    color = I("red"),
    symbol= I(4),
    name=c("gBLUP with Inf")
    )
    htmltools::save_html(p, "Interactive.GS.html")

}



