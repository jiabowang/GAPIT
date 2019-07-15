`GAPIT.Interactive.GS`<-
function(model_store=NULL,Y=NULL,myGD=NULL,myGM=NULL,myKI=NULL,myY=NULL,myCV=NULL,rel=NULL,h2=NULL,NQTN=NULL
  )
#model_store is the store of all model names
#Y is the real phenotype
#
{ 

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




myY=read.table(paste("gcs_",e,".txt",sep=""),head=T)
Observed=myY$observed
Predicted=myY$gBLUP
if(!require(plotly)) install.packages("plotly")
  library(plotly)

  p <- plot_ly(
    type = 'scatter',
    x = ~Observed,
    y = ~Predicted,
    data=myY,
    text = ~paste("Taxa: ",taxa,"<br>Observed: ",round(observed,4) , '<br>gBLUP:', round(gBLUP,4)),
    #size=2*y/max(y),
    color = I("red"),
    name=c("gBLUP")
    )%>%add_trace(
    type = 'scatter',
    x = ~observed,
    y = ~cBLUP,
    #data=myY,
    text = ~paste("Taxa: ",taxa,"<br>Observed: ",round(observed,4)  , '<br>cBLUP:', round(cBLUP,4)),
    #size=2*y/max(y),
    color = I("blue"),
    name=c("cBLUP")
    )%>%add_trace(
    type = 'scatter',
    x = ~observed,
    y = ~sBLUP,
    #data=myY,
    text = ~paste("Taxa: ",taxa,"<br>Observed: ",round(observed,4)  , '<br>sBLUP:', round(sBLUP,4)),
    #size=2*y/max(y),
    color = I("green"),
    name=c("sBLUP")
    )

    htmltools::save_html(p, "Interactive.GS.html")


}



