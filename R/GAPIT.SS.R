`GAPIT.SS` <-
function(DP=NULL,IC=NULL,buspred=FALSE,lmpred=TRUE){
#Object: To Sufficient Statistics (SS) for GWAS and GS
#Designed by Zhiwu Zhang
#Writen by Jiabo Wang
#Last update: Novenber 3, 2016
##############################################################################################
print("GAPIT.SS in process...")
#Define the funcitno here
Timmer=GAPIT.Timmer(Infor="GAPIT.SS")
Memory=GAPIT.Memory(Infor="GAPIT.SS")
GR=list()
GR$GVs=NULL

if(DP$SNP.test)
{
  print("GAPIT will be into GWAS approach...")
  ic_GD=IC$GD
  ic_GM=IC$GM
  ic_Y=IC$Y
  ic_KI=IC$K
  ic_PCA=IC$PCA
  Z=DP$Z

  taxa_Y=as.character(ic_Y[,1])
  Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="GAPIT.QC")
  Memory=GAPIT.Memory(Memory=Memory,Infor="GAPIT.QC")

  if(DP$kinship.algorithm!="None" & DP$kinship.algorithm!="SUPER" & is.null(Z))
  {
    Z=as.data.frame(diag(1,nrow(ic_Y)))
    Z=rbind(taxa_Y,Z)
    taxa=c('Taxa',as.character(taxa_Y))
    Z=cbind(taxa,Z)
  }
# print(head(ic_PCA))
# print(dim(DP$CV))
# print(head(DP$PC))
  if(max(ic_PCA[,2])==min(ic_PCA[,2]))ic_PCA=NULL
#print(head(ic_PCA))
# print("@@@@@")
  print(DP$kinship.algorithm)
  if(DP$kinship.algorithm%in%c("FarmCPU","BLINK","MLMM","BLINKC"))
  {
     print("The GAPIT would go into Bus...")
     Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="GAPIT.FarmCPU")
     Memory=GAPIT.Memory(Memory=Memory,Infor="GAPIT.FarmCPU")

     myBus=GAPIT.Bus(Y=ic_Y,CV=ic_PCA,Z=NULL,GK=NULL,KI=ic_KI,GD=ic_GD,GM=ic_GM,GT=IC$GT,name.of.trait=DP$name.of.trait,DP=DP,
                method=DP$kinship.algorithm,GTindex=DP$GTindex,LD=DP$LD,opt=DP$opt,N.sig=DP$N.sig,
                bin.size=DP$bin.size,bin.selection=DP$bin.selection,alpha=DP$alpha,WS=DP$WS,
                cutOff=DP$cutOff,p.threshold=DP$p.threshold,QTN.threshold=DP$QTN.threshold,FDRcut=DP$FDRcut,
                maf.threshold=DP$maf.threshold,method.GLM=DP$method.GLM,method.sub=DP$method.sub,seq.cutoff=DP$seq.cutoff,
                method.sub.final=DP$method.sub.final,method.bin=DP$method.bin,Random.model=DP$Random.model,
				        DPP=DP$DPP,file.output=DP$file.output,Multi_iter=DP$Multi_iter,num_regwas=DP$num_regwas,bin.regwas= DP$bin.regwas)
     GWAS=myBus$GWAS
     Pred=myBus$Pred
     p=GWAS[,4]
     # print(head(GWAS))
     # print(head(ic_GM))
     myBus$seqQTN=match(as.character(GWAS[p<(DP$cutOff/length(p)),1]),as.character(ic_GM[,1]))
     myBus$seqQTN=myBus$seqQTN[!is.na(myBus$seqQTN)]
 # BUS Prediction with gBLUP
# lmpred=TRUE

     print(myBus$seqQTN)
     if(buspred)
     {  
        X=IC$myallGD[,-1]
        for(l in lmpred)
        {
          memo=ifelse(l,"MAS","ABLUP")
          if(length(myBus$seqQTN)>DP$seq.num)
              {
               myBus$seqQTN=order(p)[1:DP$seq.num]
              }
        if(l)
        {
          print("MAS to Predict phenotype !!")
    # colnames(busCV)[1]=c("Taxa")
    # print(length(IC$GT))
          # index=as.character(IC$myallGD[,1])%in%as.character(ic_Y[,1])
    # print(cbind(ic_Y,IC$PCA))
          if(!is.null(myBus$seqQTN))
          {
          # busCV=cbind(IC$myallCV,X[,myBus$seqQTN])
            
            GD1 = as.matrix(IC$GD[,c(1,myBus$seqQTN+1)])
            GD2 = as.matrix(IC$myallGD[,c(1,myBus$seqQTN+1)])
          }else{
            numMarker=nrow(GWAS)
            sp=sort(GWAS$P.value)
            spd=abs(DP$cutOff-sp)
            index_fdr=grep(min(spd),spd)[1]
            FDRcutoff=DP$cutOff*index_fdr/numMarker
            seqQTN=as.numeric(rownames(GWAS[GWAS$P.value<FDRcutoff,]))
          # busCV=cbind(IC$myallCV,X[,seqQTN])
            GD1 = as.matrix(IC$GD[,c(1,seqQTN+1)])
            GD2 = as.matrix(IC$myallGD[,c(1,seqQTN+1)])
          }
          if(!is.null(IC$myallCV)) 
          {
            CV1=as.matrix(IC$PCA[,-1])
            CV2=as.matrix(IC$myallCV[,-1])
            colnames(GD1)[1]="Taxa"
            colnames(GD2)[1]="Taxa"
            # BLUE=rep(NA,nrow(IC$myallGD))
            print("The dimension of CV in lm model :")
            print(dim(CV1))
            # print(dim(GD1))
    # print(colnames(IC$PCA))
    # print(colnames(GD1))
            # XCV1=cbind(CV1, GD1)
            XCV1=merge(IC$PCA,GD1,by.x=colnames(IC$PCA)[1],by.y=colnames(GD1)[1])
            XCV1=matrix(as.numeric(as.matrix(XCV1[,-1])),nrow(XCV1))
            cv.licols=GAPIT.Licols(X=XCV1)
            XCV1=cv.licols$Xsub
            # print(dim(XCV1))
            # print(dim(ic_Y))
            mylm = stats::lm(ic_Y[,2] ~XCV1)
            lm.coeff=mylm$coefficients
            if(stats::var(CV2[,1])==0)lm.coeff=lm.coeff[-2]
            print(lm.coeff)
            XCV2=merge(IC$myallCV, GD2,by.x=colnames(IC$myallCV)[1],by.y=colnames(GD2)[1])
            XCV=matrix(as.numeric(as.matrix(XCV2[,-1])),nrow(XCV2))
            XCV=XCV[,cv.licols$idx]
            XCV=cbind(1,as.matrix(XCV))
            # print(dim(XCV2))
#CV.Extragenetic specified
            QTN.gs=ncol(GD2)-1
            CV.Extragenetic=DP$CV.Extragenetic
            if(ncol(XCV)>(1+CV.Extragenetic))XCVI=XCV[,c((2+CV.Extragenetic):(ncol(XCV)-QTN.gs))]
            XCVN=XCV[,c(1:(1+CV.Extragenetic)),drop=FALSE]
            if(QTN.gs!=0)XCVqtn=XCV[,c((ncol(XCV)-QTN.gs):ncol(XCV))]
            if(ncol(XCV)>(1+CV.Extragenetic))beta.I=lm.coeff[c((2+CV.Extragenetic):(ncol(XCV)-QTN.gs))]
            beta.N=lm.coeff[c(1:(1+CV.Extragenetic))]
            if(QTN.gs!=0)beta.QTN=lm.coeff[c((ncol(XCV)-QTN.gs):ncol(XCV))]
            BLUE.N=XCVN%*%beta.N
            BLUE.QTN=rep(0,length(BLUE.N))    
            if(QTN.gs!=0)BLUE.QTN=XCVqtn%*%beta.QTN
            BLUE.I=rep(0,length(BLUE.N))
            if(ncol(XCV)>(1+CV.Extragenetic))BLUE.I=XCVI%*%beta.I
            BLUE=cbind(BLUE.N,BLUE.I,BLUE.QTN)
            # print(dim(BLUE))
            BLUE=data.frame(cbind(data.frame(XCV2[,1]),data.frame(BLUE)))
            colnames(BLUE)=c("Taxa","BLUE.N","BLUE.I","QTNs")
            Group=1:nrow(BLUE)
            RefInf=rep(2,nrow(BLUE))
            # print(table(index))
            index=as.character(XCV2[,1])%in%as.character(ic_Y[,1])
            RefInf[index]=1
            ID=1:nrow(BLUE)
            BLUP=rep(0,nrow(BLUE))
            PEV=rep(0,nrow(BLUE))
            BB= cbind(BLUE,Group,RefInf,ID,BLUP,PEV)
            gBreedingValue=BB[,3]+BB[,4]+BB[,8]
            Prediction=BB[,2]+BB[,3]+BB[,4]+BB[,8]
            Pred=cbind(BB,gBreedingValue,Prediction)
            colnames(Pred)=c("Taxa","BLUE.N","BLUE.I","QTNs","Group","RefInf","ID","BLUP","PEV","gBreedingValue","Prediction")

          }else{
            busCV=cbind(as.data.frame(DP$GD[,1]),X[,myBus$seqQTN])
            CV1=NULL
            Group=1:nrow(IC$myallCV)
            RefInf=rep(2,nrow(IC$myallCV))
            RefInf[index]=1
            ID=1:nrow(IC$myallCV)
            BLUP=NA
            PEV=NA
            BLUE=NA
            print("The dimension of CV in lm model :")
            print(dim(CV1))
            # print(dim(GD1))
            mylm = stats::lm(ic_Y[!is.na(ic_Y[,2]),2] ~GD1[,-1])
            # print(stats::predict(mylm,as.data.frame(cbind(IC$myallCV[,-1],GD2))))
            Pred = cbind(as.character(DP$GD[,1]),Group,RefInf,ID,BLUP,PEV,BLUE,stats::predict(mylm,as.data.frame(cbind(IC$myallCV[,-1],GD2))))
            colnames(Pred)=c("Taxa","Group","RefInf","ID","BLUP","PEV","BLUE","Prediction")   
          }   
    # print(dim(CV1))
    # print(table(index))
          print("Linear Regression to Predict phenotype Done !!")  
        }else{
          print("GAGBLUP to Predict phenotype !!")
          print(is.null(IC$myallCV))
          if(!is.null(IC$myallCV)) 
          {
            com.taxa=intersect(as.character(IC$myallCV[,1]),as.character(DP$GD[,1]))
            taxa_CV=as.character(IC$myallCV[,1])
            taxa_g=as.character(DP$GD[,1])
     # print(length(taxa_comall))
            CV1=IC$myallCV[taxa_CV%in%com.taxa,]
            CV1 <- CV1[match(com.taxa,as.character(CV1[,1])),]
            # CV1 = IC$myallCV
            ablup.GD=DP$GD[taxa_g%in%com.taxa,]
            ablup.GD=ablup.GD[match(com.taxa,as.character(ablup.GD[,1])),]           
            ablup.X=ablup.GD[,-1]
            # CV1=as.matrix(CV1[match(com.taxa,as.character(CV1[,1])),])
            if(!is.null(myBus$seqQTN))
            {
               busCV=cbind(CV1,ablup.X[,myBus$seqQTN])
               # busCV=ablup.GD[,c(1,myBus$seqQTN)]

            }else{
               busCV=CV1
            }
          }else{
            ablup.GD=DP$GD
            ablup.X=ablup.GD[,-1]
            busCV=cbind(as.data.frame(ablup.GD[,1]),ablup.X[,myBus$seqQTN])
          }
          # print(dim(GWAS))
          pv=GWAS$effect
          # pv=GWAS$P.value

          threshold <- quantile(abs(pv), 0.3)
          # noneff=as.numeric(rownames(GWAS[GWAS$P.value>DP$cutOff,]))
          
          licols=TRUE
          if(licols)
            {
              # pv.index=pv<threshold
              pv.index=abs(pv)>threshold
              # pv.index=pv!=0

              print(table(pv.index))
              # qtn.index=1:length(pv)%in%myBus$seqQTN
              # print(length(pv.index))
              # print(length(qtn.index))
              ablup.X=ablup.X[,pv.index]
            }
            # print(dim(ablup.X))
          if(is.null(DP$KI))
          {
            KI= GAPIT.kinship.VanRaden(snps=as.matrix(ablup.X))
            colnames(KI)=as.character(ablup.GD[,1])
            busKI=cbind(as.data.frame(ablup.GD[,1]),KI)
            colnames(busKI)[1]=c("Taxa")
          }else{
            busKI=DP$KI
          }
          busCV=cbind(as.data.frame(busCV[,1]),matrix(as.numeric(as.matrix(busCV[,-1])),nrow(busCV),ncol(busCV)-1))

          print("The dimension of CV and phenotype in ABLUP model :")
          print(dim(busCV))
   # print(head(busCV))
   # print(apply(busCV[,-1],2,sum))
          # print(dim(ic_Y))
          busGAPIT=GAPIT(
                  Y=ic_Y,
                  KI=busKI,
                  CV=busCV,
                  CV.Extragenetic=DP$CV.Extragenetic,
                  QTN.gs=ncol(busCV)-ncol(CV1),
                  model="gBLUP",
                  file.output=F)
          Pred=busGAPIT$Pred
          print("ABLUP Predict phenotype Done!!")
        }#if lmpred
        if(DP$file.output) 
        {
          utils::write.csv(Pred,paste("GAPIT.Association.Prediction_results.",DP$name.of.trait,".",memo,".csv",sep=""), row.names = FALSE)
        }
        }#lmpred0
     }#buspred
     va=myBus$vg
     ve=myBus$ve
     h2=va/(va+ve)
     mc=NULL
     bc=NULL
     mp=NULL
     TV=NULL
     Compression=NULL
     GVs=myBus$GVs
  }else{
    print("The GAPIT would go into Main...")
    Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="GAPIT.Main")
    Memory=GAPIT.Memory(Memory=Memory,Infor="GAPIT.Main")

    GT=as.matrix(ic_GD[,1])
    if(DP$PCA.total==0) ic_PCA=NULL
    gapitMain <- GAPIT.Main(Y=ic_Y,
                         GD=IC$GD[,-1],
                         allGD=IC$myallGD[,-1],
                         allCV=IC$myallCV,
                         GM=DP$GM,
                         KI=ic_KI,
                         CV=IC$PCA,
                         CV.Extragenetic=DP$CV.Extragenetic,
                         GP=DP$GP,
                         GK=DP$GK,
                         SNP.P3D=DP$SNP.P3D,
                         kinship.algorithm=DP$kinship.algorithm,
						             bin.from=DP$bin.from,
						             bin.to=DP$bin.to,
						             bin.by=DP$bin.by,
						             inclosure.from=DP$inclosure.from,
						             inclosure.to=DP$inclosure.to,
						             inclosure.by=DP$inclosure.by,
				                 group.from=DP$group.from,
						             group.to=DP$group.to,
						             group.by=DP$group.by,
						             kinship.cluster=DP$kinship.cluster,
						             kinship.group=DP$kinship.group,
						             name.of.trait=DP$name.of.trait,
                         file.path=DP$file.path,
						             file.from=DP$file.from, 
						             file.to=DP$file.to, 
						             file.total=DP$file.total, 
						             file.fragment = DP$file.fragment, 
						             file.G=DP$file.G,
						             file.Ext.G=DP$file.Ext.G,
						             file.GD=DP$file.GD, 
						             file.GM=DP$file.GM, 
						             file.Ext.GD=DP$file.Ext.GD,
						             file.Ext.GM=DP$file.Ext.GM, 
                         SNP.MAF= DP$SNP.MAF,
						             FDR.Rate = DP$FDR.Rate,
						             SNP.FDR=DP$SNP.FDR,
						             SNP.effect=DP$SNP.effect,
						             SNP.impute=DP$SNP.impute,
						             PCA.total=DP$PCA.total,
						             #GAPIT.Version=GAPIT.Version,
                         GT=IC$GT, 
						             SNP.fraction = DP$SNP.fraction, 
						             seed =DP$seed, 
						             BINS = DP$BINS,
						             SNP.test=DP$SNP.test,DPP=DP$DPP, 
						             SNP.permutation=DP$SNP.permutation,
                         LD.chromosome=DP$LD.chromosome,
#						             LD.location=LD.location,
#						             LD.range=LD.range,
#						             SNP.CV=SNP.CV,
						             SNP.robust=DP$SNP.robust,
						             model=DP$model,
                         genoFormat="EMMA",
						             hasGenotype=TRUE,
						             byFile=FALSE,
						             fullGD=TRUE,
						             PC=DP$PC,
						             GI=ic_GM,
						             Timmer = DP$Timmer, 
						             Memory = DP$Memory,
                         sangwich.top=DP$sangwich.top,
						             sangwich.bottom=DP$sangwich.bottom,
						             QC=DP$QC,
						             GTindex=DP$GTindex,
						             LD=DP$LD,
						             file.output=FALSE,
						             cutOff=DP$cutOff, 
						             GAPIT3.output=DP$file.output,
                         Model.selection = DP$Model.selection, 
						             Create.indicator = DP$Create.indicator,
						             # QTN=DP$QTN, 
						             # QTN.round=DP$QTN.round,
						             # QTN.limit=DP$QTN.limit,
						             #QTN.update=QTN.update, 
						             # QTN.method=DP$QTN.method,
						             Major.allele.zero=DP$Major.allele.zero,
						             NJtree.group=DP$NJtree.group,
						             NJtree.type=DP$NJtree.type,
						             plot.bin=DP$plot.bin, 
                         QTN.position=DP$QTN.position,
						             plot.style=DP$plot.style,
						             SUPER_GS=DP$SUPER_GS)  
    GWAS=gapitMain$GWAS
    if(DP$Random.model&DP$file.output)GR=GAPIT.RandomModel(Y=ic_Y,X=IC$GD[,-1],GWAS=GWAS,CV=IC$PCA,cutOff=DP$cutOff,name.of.trait=DP$name.of.trait,N.sig=DP$N.sig,GT=IC$GT)
    Pred=gapitMain$Pred
    va=NA#gapitMain$vg
    ve=NA#gapitMain$ve
    h2=gapitMain$h2
    mc=gapitMain$effect.snp
    bc=gapitMain$effect.cv
    mp=gapitMain$P
    TV=gapitMain$TV
    Compression=gapitMain$Compression
    GVs=GR$GVs
  }#!DP$kinship.algorithm%in%c("FarmCPU","MLMM","BLINK","BLINKC")
myPower=NULL
if(!is.null(GWAS))myPower=GAPIT.Power(WS=DP$WS, alpha=DP$alpha, maxOut=DP$maxOut,seqQTN=DP$QTN.position,GM=DP$GM,GWAS=GWAS)
  return (list(GWAS=GWAS,Pred=Pred,FDR=myPower$FDR,Power=myPower$Power,
  Power.Alpha=myPower$Power.Alpha,alpha=myPower$alpha,h2=h2,va=va,ve=ve,
  mc=mc,bc=bc,mp=mp,TV=TV,Compression=Compression,
  Timmer=Timmer,Memory=Memory,GVs=GVs))
}else{
# Here is Genomic Prediction function
  print("GAPIT will be into GS approach...")
gapitMain <- GAPIT.Main(Y=IC$Y,
                        GD=IC$GD[,-1],
                        allGD=IC$myallGD[,-1],
                        GM=DP$GM,
                        KI=IC$KI,
                        Z=DP$Z,
                        CV=IC$PCA,
                        allCV=IC$myallCV,
                        CV.Extragenetic=DP$CV.Extragenetic,
                        GP=DP$GP,
                        GK=DP$GK,
                        SNP.P3D=DP$SNP.P3D,
                        kinship.algorithm=DP$kinship.algorithm,
                        bin.from=DP$bin.from,
                        bin.to=DP$bin.to,
                        bin.by=DP$bin.by,
                        inclosure.from=DP$inclosure.from,
                        inclosure.to=DP$inclosure.to,
                        inclosure.by=DP$inclosure.by,
                        group.from=DP$group.from,
                        group.to=DP$group.to,
                        group.by=DP$group.by,
                        kinship.cluster=DP$kinship.cluster,
                        kinship.group=DP$kinship.group,
                        name.of.trait=DP$name.of.trait,
                        file.path=DP$file.path,
                        file.from=DP$file.from,
                        file.to=DP$file.to,
                        file.total=DP$file.total,
                        file.fragment = DP$file.fragment,
                        file.G=DP$file.G,
                        file.Ext.G=DP$file.Ext.G,
                        file.GD=DP$file.GD,
                        file.GM=DP$file.GM, 
                        file.Ext.GD=DP$file.Ext.GD,
                        file.Ext.GM=DP$file.Ext.GM, 
                        SNP.MAF= DP$SNP.MAF,
                        FDR.Rate = DP$FDR.Rate,
                        SNP.FDR=DP$SNP.FDR,
                        SNP.effect=DP$SNP.effect,
                        SNP.impute=DP$SNP.impute,
                        PCA.total=DP$PCA.total,
                        #GAPIT.Version=GAPIT.Version,
                        GT=DP$GT, 
                        SNP.fraction = DP$SNP.fraction,
                        seed =DP$ seed,
                        BINS = DP$BINS,
                        SNP.test=DP$SNP.test,
                        DPP=DP$DPP,
                        SNP.permutation=DP$SNP.permutation,
                        LD.chromosome=DP$LD.chromosome,
                        QTN.gs=DP$QTN.gs,
                        #LD.location=LD.location,
                        #LD.range=LD.range,
                        #SNP.CV=SNP.CV,
                        SNP.robust=DP$SNP.robust,
                        model=DP$model,
                        genoFormat="EMMA",
                        hasGenotype=TRUE,
                        byFile=FALSE,
                        fullGD=TRUE,
                        PC=DP$PC,
                        GI=DP$GI,
                        Timmer = DP$Timmer, 
                        Memory = DP$Memory,
                        GAPIT3.output=DP$file.output,
                        sangwich.top=DP$sangwich.top,
                        sangwich.bottom=DP$sangwich.bottom,
                        QC=DP$QC,GTindex=DP$GTindex,
                        LD=DP$LD,file.output=FALSE,
                        cutOff=DP$cutOff, 
                        Model.selection = DP$Model.selection, 
                        Create.indicator = DP$Create.indicator,
                        # QTN=DP$QTN,
                        # QTN.round=DP$QTN.round,
                        # QTN.limit=DP$QTN.limit, 
                        #QTN.update=QTN.update, 
                        # QTN.method=DP$QTN.method, 
                        Major.allele.zero=DP$Major.allele.zero,
                        NJtree.group=DP$NJtree.group,
                        NJtree.type=DP$NJtree.type,
                        plot.bin=DP$plot.bin, 
                        QTN.position=DP$QTN.position,
                        plot.style=DP$plot.style,
                        SUPER_GS=TRUE
                        )  
#print(str(gapitMain))
GWAS=gapitMain$GWAS
Pred=gapitMain$Pred
#print(head(Pred))
va=NA#gapitMain$vg
ve=NA#gapitMain$ve
h2=gapitMain$h2
mc=gapitMain$effect.snp
bc=gapitMain$effect.cv
mp=gapitMain$P
Compression=gapitMain$Compression
GAPIT.Compression.Visualization(Compression = Compression, name.of.trait = DP$name.of.trait,file.output=DP$file.output)
# # print(list(GWAS=GWAS,Pred=Pred,FDR=NULL,Power=NULL,
#   Power.Alpha=NULL,alpha=NULL,h2=h2,va=va,ve=ve,Compression=Compression,
#   mc=mc,bc=bc,mp=mp,TV=gapitMain$TV,
#   Timmer=Timmer,Memory=Memory))
return (list(GWAS=GWAS,Pred=Pred,FDR=NULL,Power=NULL,
  Power.Alpha=NULL,alpha=NULL,h2=h2,va=va,ve=ve,Compression=Compression,
  mc=mc,bc=bc,mp=mp,TV=gapitMain$TV,
  Timmer=Timmer,Memory=Memory))
}#end of SNP.TEST

}  #end of GAPIT.SS function
#=============================================================================================

