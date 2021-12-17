
# ' GAPIT.DP
# ' @description 
# ' To Data and Parameters
# ' @author Zhiwu Zhang and Jiabo Wang
# ' @export
`GAPIT.DP` <-
function(G=NULL,GD=NULL,GM=NULL,KI=NULL,Z=NULL,CV=NULL,CV.Inheritance=NULL,GP=NULL,GK=NULL,
                group.from=30 ,group.to=1000000,group.by=10,DPP=100000, 
                kinship.cluster="average", kinship.group='Mean',kinship.algorithm="VanRaden",                                                    
                bin.from=10000,bin.to=10000,bin.by=10000,inclosure.from=10,inclosure.to=10,inclosure.by=10,
                SNP.P3D=TRUE,SNP.effect="Add",SNP.impute="Middle",PCA.total=0, SNP.fraction = 1, seed = 123, BINS = 20,SNP.test=TRUE, 
                    SNP.MAF=0,FDR.Rate = 1, SNP.FDR=1,SNP.permutation=FALSE,SNP.CV=NULL,SNP.robust="GLM", NJtree.group=NULL,NJtree.type=c("fan","unrooted"),plot.bin=10^6,PCA.col=NULL,PCA.3d=FALSE,
                file.from=1, file.to=1, file.total=NULL, file.fragment = 99999,file.path=NULL,Inter.Plot=FALSE,Inter.type=c("m","q"),
                file.G=NULL, file.Ext.G=NULL,file.GD=NULL, file.GM=NULL, file.Ext.GD=NULL,file.Ext.GM=NULL, 
                ngrid = 100, llim = -10, ulim = 10, esp = 1e-10, Multi_iter=FALSE,num_regwas=10,FDRcut=FALSE,
                LD.chromosome=NULL,LD.location=NULL,LD.range=NULL, p.threshold=NA,QTN.threshold=0.01,maf.threshold=0.03,
                sangwich.top=NULL,sangwich.bottom=NULL,QC=TRUE,GTindex=NULL,LD=0.1,opt="extBIC",
                file.output=FALSE,cutOff=0.01, Model.selection = FALSE,output.numerical = FALSE,Random.model=FALSE,
                output.hapmap = FALSE, Create.indicator = FALSE,QTN=NULL, QTN.round=1,QTN.limit=0, QTN.update=TRUE, QTN.method="Penalty", Major.allele.zero = FALSE,
        method.GLM="fast.lm",method.sub="reward",method.sub.final="reward",method.bin="static",bin.size=c(1000000),bin.selection=c(10,20,50,100,200,500,1000),
        memo="",Prior=NULL,ncpus=1,maxLoop=3,threshold.output=.01, WS=c(1e0,1e3,1e4,1e5,1e6,1e7),alpha=c(.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),maxOut=100,QTN.position=NULL,
        converge=1,iteration.output=FALSE,acceleration=0,iteration.method="accum",PCA.View.output=TRUE,Geno.View.output=TRUE,plot.style="Oceanic",SUPER_GD=NULL,SUPER_GS=FALSE,CG=NULL,model="MLM"){
#Object: To Data and Parameters  
#Designed by Zhiwu Zhang
#Writen by Jiabo Wang
#Last update: Novenber 3, 2016
##############################################################################################
print("GAPIT.DP in process...")
#Judge phenotype  genotype and GAPIT logical
#print(file.from)
#print(kinship.algorithm)
#print(NJtree.group)
myGenotype<-GAPIT.Genotype(G=G,GD=GD,GM=GM,KI=KI,PCA.total=PCA.total,kinship.algorithm=kinship.algorithm,SNP.fraction=SNP.fraction,SNP.test=FALSE,
                file.path=file.path,file.from=file.from, file.to=file.to, file.total=file.total, file.fragment = file.fragment, file.G=file.G, 
                file.Ext.G=file.Ext.G,file.GD=file.GD, file.GM=file.GM, file.Ext.GD=file.Ext.GD,file.Ext.GM=file.Ext.GM,
                SNP.MAF=SNP.MAF,FDR.Rate = FDR.Rate,SNP.FDR=SNP.FDR,SNP.effect=SNP.effect,SNP.impute=SNP.impute,NJtree.group=NJtree.group,NJtree.type=NJtree.type,
                LD.chromosome=LD.chromosome,LD.location=LD.location,LD.range=LD.range,
                GP=GP,GK=GK,bin.size=NULL,inclosure.size=NULL, 
                sangwich.top=sangwich.top,sangwich.bottom=sangwich.bottom,GTindex=NULL,file.output=file.output, Create.indicator = Create.indicator, Major.allele.zero = Major.allele.zero,Geno.View.output=Geno.View.output,PCA.col=PCA.col,PCA.3d=PCA.3d)

# }

KI=myGenotype$KI
PC=myGenotype$PC
print(dim(PC))

genoFormat=myGenotype$genoFormat
hasGenotype=myGenotype$hasGenotype
byFile=myGenotype$byFile
fullGD=myGenotype$fullGD
GD=myGenotype$GD
GI=myGenotype$GI

GT=myGenotype$GT
G=myGenotype$G
chor_taxa=myGenotype$chor_taxa

#if G exist turn to GD and GM

if(output.numerical){
    utils::write.table(GD,  "GAPIT.Genotype.Numerical.txt",
    quote = FALSE, sep = "\t", row.names = TRUE,col.names = NA)
}

if(output.hapmap){
    utils::write.table(myGenotype$G,  "GAPIT.Genotype.hmp.txt",
                       quote = FALSE, sep = "\t", row.names = FALSE,
                       col.names = FALSE)
}

rownames(GD)=GT
colnames(GD)=GI[,1]
GD=cbind(as.data.frame(GT),GD)

  print("GAPIT.DP accomplished successfully for multiple traits. Results are saved")
  return (list(Y=NULL,G=G,GD=GD,GM=GI,KI=KI,Z=Z,CV=CV,CV.Inheritance= CV.Inheritance,GP=GP,GK=GK,PC=PC,GI=GI,
                group.from= group.from ,group.to= group.to,group.by= group.by,DPP= DPP, name.of.trait=NULL,
                kinship.cluster= kinship.cluster, kinship.group= kinship.group,kinship.algorithm= kinship.algorithm,NJtree.group=NJtree.group,NJtree.type=NJtree.type,PCA.col=PCA.col,PCA.3d=PCA.3d,                                              
                bin.from= bin.from,bin.to= bin.to,bin.by= bin.by,inclosure.from= inclosure.from,inclosure.to= inclosure.to,inclosure.by= inclosure.by,opt=opt,
                SNP.P3D= SNP.P3D,SNP.effect= SNP.effect,SNP.impute= SNP.impute,PCA.total= PCA.total, SNP.fraction = SNP.fraction, seed = seed, 
                BINS = BINS,SNP.test=SNP.test, SNP.MAF= SNP.MAF,FDR.Rate = FDR.Rate, SNP.FDR= SNP.FDR,SNP.permutation= SNP.permutation,
                SNP.CV= SNP.CV,SNP.robust= SNP.robust, file.from= file.from, file.to=file.to, file.total= file.total, file.fragment = file.fragment,file.path= file.path, 
                file.G= file.G, file.Ext.G= file.Ext.G,file.GD= file.GD, file.GM= file.GM, file.Ext.GD= file.Ext.GD,file.Ext.GM= file.Ext.GM, 
                ngrid = ngrid, llim = llim, ulim = ulim, esp = esp,Inter.Plot=Inter.Plot,Inter.type=Inter.type,
                LD.chromosome= LD.chromosome,LD.location= LD.location,LD.range= LD.range,Multi_iter=Multi_iter,
                sangwich.top= sangwich.top,sangwich.bottom= sangwich.bottom,QC= QC,GTindex= GTindex,LD= LD,GT=GT,
                file.output= file.output,cutOff=cutOff, Model.selection = Model.selection,output.numerical = output.numerical,
                output.hapmap = output.hapmap, Create.indicator = Create.indicator,Random.model=Random.model,
				 QTN= QTN, QTN.round= QTN.round,QTN.limit= QTN.limit, QTN.update= QTN.update, QTN.method= QTN.method, Major.allele.zero = Major.allele.zero,
               method.GLM= method.GLM,method.sub= method.sub,method.sub.final= method.sub.final,
               method.bin= method.bin,bin.size= bin.size,bin.selection= bin.selection,FDRcut=FDRcut,
        		memo= memo,Prior= Prior,ncpus=1,maxLoop= maxLoop,threshold.output= threshold.output,
        		WS= WS,alpha= alpha,maxOut= maxOut,QTN.position= QTN.position, converge=1,iteration.output= iteration.output,acceleration=0,
        		iteration.method= iteration.method,PCA.View.output= PCA.View.output, 
                p.threshold=p.threshold,QTN.threshold=QTN.threshold,
                maf.threshold=maf.threshold,chor_taxa=chor_taxa,num_regwas=num_regwas,
        		Geno.View.output= Geno.View.output,plot.style= plot.style,SUPER_GD= SUPER_GD,SUPER_GS= SUPER_GS,CG=CG,plot.bin=plot.bin))
}  #end of GAPIT DP function
#=============================================================================================

