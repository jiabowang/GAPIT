`GAPIT.FDR.TypeI` <-
function(WS=c(1e0,1e3,1e4,1e5), GM=NULL,seqQTN=NULL,GWAS=NULL,maxOut=100,MaxBP=1e10){
    #Object: To evaluate power and FDR for the top (maxOut) positive interval defined by WS
    #Input: WS- window size
    #Input: GM - m by 3  matrix for SNP name, chromosome and BP
    #Input: seqQTN - s by 1 vecter for index of QTN on GM (+1 for GDP column wise)
    #Input: GWAS - SNP,CHR,BP,P,MAF
    #maxOut: maximum number of rows to report
    #Requirement: None
    #Output: Table and Plots
    #Authors: Xiaolei Liu & Zhiwu Zhang
    # Date  start: April 2, 2013
    # Last update: Mar 16, 2016
    ##############################################################################################
    #print("GAPIT.Power Started")
    if(is.null(seqQTN) | is.null(GM)) return(list(Power=NULL,FDR= NULL,TypeI= NULL,False= NULL,AUC.FDR= NULL,AUC.T1= NULL))
    
    #store number fdr and t1 records
    NQTN=length(seqQTN)
    table=array(NA,dim=c(NQTN,2*length(WS)))
    fdrtable=array(NA,dim=c(NQTN,2*length(WS)))
    t1table=array(NA,dim=c(NQTN,2*length(WS)))
    cutoff=array(NA,dim=c(length(WS),NQTN))
    cut=array(NA,dim=c(1,NQTN))
    #-----------------FDR and Power analysis-------------------------
    #Information needed: GWAS,myGM and QTN(r)
    GWAS=GWAS[order(GWAS[,2],GWAS[,3]),]
    GWAS[is.na(GWAS[,4]),4]=1
    QTN.list=sort(GWAS[seqQTN,4])
    powerlist=seq(1/length(QTN.list),1,length.out=length(QTN.list))
    #calculate number of false positives in each WS
    total.index=1:nrow(GM)
    
    theWS=1
    for (theWS in 1:length(WS)){
        wsws=WS[theWS]
        qtn.pool=ceiling((as.numeric(GWAS[seqQTN,2])*MaxBP+as.numeric(GWAS[seqQTN,3]))/(2*wsws))
        bonf.pool=ceiling((GWAS[total.index,2]*MaxBP+GWAS[total.index,3])/(2*wsws))
        false.number=length(levels(factor(bonf.pool[!(bonf.pool%in%qtn.pool)])))
        for(j in 1:length(qtn.pool)){
            pbin=min(GWAS[bonf.pool==qtn.pool[j],4])
            cut[,j]=pbin
        }
        if(theWS==1){
            totalfalse=false.number
        }else{
            totalfalse=c(totalfalse,false.number)
        }
        cutoff[theWS,]=sort(cut)
    }
    #Calculate FDR and T1
    for(j in 1:ncol(cutoff)){
        theWS=1
        for (theWS in 1:length(WS)){
            p.index=which(GWAS[,4]<=cutoff[theWS,j])
            wsws=WS[theWS]
            qtn.pool=ceiling((GWAS[seqQTN,2]*MaxBP+GWAS[seqQTN,3])/(2*wsws))
            bonf.pool=ceiling((GWAS[p.index,2]*MaxBP+GWAS[p.index,3])/(2*wsws))
            qtn.number=length(levels(factor(bonf.pool[bonf.pool%in%qtn.pool])))
            false.number=length(levels(factor(bonf.pool[!(bonf.pool%in%qtn.pool)])))
            if(theWS==1){
                final=false.number
                final.fdr=false.number/(qtn.number+false.number)
                final.t1=false.number/totalfalse[theWS]
            }else{
                record=false.number
                record.fdr=false.number/(qtn.number+false.number)
                record.t1=false.number/totalfalse[theWS]
                final=c(final,record)
                final.fdr=c(final.fdr,record.fdr)
                final.t1=c(final.t1,record.t1)
            }
        }
        #record FDR and T1
        if(j==1){
            number.record=final
            fdr.record=final.fdr
            t1.record=final.t1
        }else{
            number.record=rbind(number.record,final)
            fdr.record=rbind(fdr.record,final.fdr)
            t1.record=rbind(t1.record,final.t1)
        }
        
    }
    
    table=number.record
    fdrtable=fdr.record
    t1table=t1.record
    #AUC
    auc.final.fdr=NULL
    auc.final.t1=NULL
    for (theWS in 1:length(WS)){
        auc.fdr=GAPIT.AUC(beta=powerlist,alpha=fdrtable[,theWS])
        auc.t1=GAPIT.AUC(beta=powerlist,alpha=t1table[,theWS])
        auc.final.fdr=c(auc.final.fdr,auc.fdr)
        auc.final.t1=c(auc.final.t1,auc.t1)
    }
    return(list(P=cutoff,Power=powerlist,FDR=fdrtable,TypeI=t1table,False=table,AUC.FDR=auc.final.fdr,AUC.T1=auc.final.t1))
    
}#end of `GAPIT.FDR.TypeI`
#=============================================================================================


