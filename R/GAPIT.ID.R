#' GAPIT.ID
#'
#' @description 
#' GAPIT.ID
#'
#' @param DP param a list (118 elements?)
#' @param IC param a list (9 elements?)
#' @param SS param a list (17 elements?)
#' @param RS param
#' @param cutOff param
#' @param DPP param
#' @param Create.indicator param
#' @param FDR.Rate param
#' @param QTN.position param
#' @param plot.style param
#' @param file.output param
#' @param SNP.MAF param
#' @param CG param
#' @param plot.bin param
#'
#'
#' @return 
#' An invisible NULL.
#'
#' @author Zhiwu Zhang and Jiabo Wang
#'
#' @export
`GAPIT.ID` <- function(
  DP=NULL,
  IC=NULL,
  SS=NULL,
  RS=NULL,
  cutOff=0.01,
  DPP=100000,
  Create.indicator=FALSE,
  FDR.Rate = 1,
  QTN.position=NULL,
  plot.style="Oceanic",
  file.output=TRUE,
  SNP.MAF=0,
  CG=NULL,
  plot.bin=10^9 ){
#Object: To Interpretation and Diagnoses 
#Designed by Zhiwu Zhang
#Writen by Jiabo Wang
#Last update: Novenber 3, 2016
##############################################################################################
print("GAPIT.ID in process...")
#Define the funcitno here

if(is.null(DP)&is.null(IC))#inputdata is other method result
{

GWAS=RS
  GI=RS[,1:3]
  GI=GI[order(GI[,2]),]
  GI=GI[order(GI[,1]),]
  #print(QTN.position)
  ps=RS[,4]
  nobs=nrow(RS)
  if(ncol(RS)>4)
   {maf=RS[,5]
    maf_pass=TRUE
   }
  if(ncol(RS)<5)
   {maf_pass=FALSE
   maf=0.5
   }
  rsquare_base=rep(NA,length(ps))
  rsquare=rep(NA,length(ps))
  df=rep(NA,length(nobs))
  tvalue=rep(NA,length(nobs))
  stderr=rep(NA,length(nobs))
  effect.est=rep(NA,length(nobs))

  if(is.na(maf[1]))  maf=matrix(.5,nrow(GWAS),1)
  print("Filtering SNPs with MAF..." )
	index=maf>=SNP.MAF	     
	PWI.Filtered=cbind(GI,ps,maf,nobs,rsquare_base,rsquare)#[index,]
	colnames(PWI.Filtered)=c("SNP","Chromosome","Position ","P.value", "maf", "nobs", "Rsquare.of.Model.without.SNP","Rsquare.of.Model.with.SNP")
  if(!is.null(PWI.Filtered))
  {
  print("Calculating FDR..." )
  PWIP <- GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure(PWI = PWI.Filtered, FDR.Rate = FDR.Rate, FDR.Procedure = "BH")
    print("QQ plot..." )
  if(file.output) GAPIT.QQ(P.values = ps, name.of.trait = name.of.trait,DPP=DPP)
   print("Manhattan plot (Genomewise)..." )
 if(file.output) GAPIT.Manhattan(GI.MP = cbind(GI[,-1],ps), name.of.trait = name.of.trait, DPP=DPP, plot.type = "Genomewise",cutOff=DP$cutOff,seqQTN=QTN.position,plot.style=plot.style,plot.bin=plot.bin)

 print("Manhattan plot (Chromosomewise)..." )

  #if(file.output) GAPIT.Manhattan(GI.MP = PWIP$PWIP[,2:4], name.of.trait = name.of.trait, DPP=DPP, plot.type = "Chromosomewise",cutOff=cutOff)
 if(file.output) GAPIT.Manhattan(GI.MP = cbind(GI[,-1],ps), name.of.trait = name.of.trait, DPP=DPP, plot.type = "Chromosomewise",cutOff=DP$cutOff,plot.bin=plot.bin)

  #Association Table
  print("Association table..." )
  #print(head(cbind(GI[,-1],ps)))
  #print(head)
  #GAPIT.Table(final.table = PWIP$PWIP, name.of.trait = name.of.trait,SNP.FDR=SNP.FDR)
 # GWAS=PWIP$PWIP[PWIP$PWIP[,9]<=DP$SNP.FDR,]
 # print(head(GWAS))
  print("Joining tvalue and stderr" )
  
        DTS=cbind(GI,df,tvalue,stderr,effect.est)
        colnames(DTS)=c("SNP","Chromosome","Position","DF","t Value","std Error","effect")	

  print("Creating ROC table and plot" )

if(file.output){
  if( !is.null(DP) ){
    ys <- nrow(DP$G)
  }
  myROC=GAPIT.ROC(t=tvalue,se=stderr,Vp=stats::var(ys),trait=name.of.trait)
}
  print("ROC table and plot created" )
  print("MAF plot..." )
if(file.output&maf_pass) myMAF1=GAPIT.MAF(MAF=maf,P=ps,E=NULL,trait=name.of.trait)
  if(file.output){
   utils::write.table(GWAS, paste("GAPIT.", name.of.trait, ".GWAS.Results.csv", sep = ""), quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
   utils::write.table(DTS, paste("GAPIT.", name.of.trait, ".Df.tValue.StdErr.csv", sep = ""), quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
   #if(!byPass) write.table(GWAS.2, paste("GAPIT.", name.of.trait, ".Allelic_Effect_Estimates.csv", sep = ""), quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
  }#end file.output
  }#end DP

}else{ #inputdata is GAPIT3 result
  cutOff=DP$cutOff
  DPP=DP$DPP
  Create.indicator=DP$Create.indicator
  FDR.Rate = DP$FDR.Rate
  QTN.position=DP$QTN.position
  plot.style=DP$plot.style
  file.output=DP$file.output
  SNP.MAF=DP$SNP.MAP
  CG=DP$CG
  plot.bin=DP$CG
  name.of.trait=DP$memo
	
GWAS=SS$GWAS
#print(head(GWAS))
Pred=SS$Pred
  GI=GWAS
  
  GI=GI[order(GI[,3]),]
  GI=GI[order(GI[,2]),]
  
  byPass=TRUE
  if(DP$kinship.algorithm%in%c("FarmCPU","MLMM","Blink","BlinkC"))byPass=FALSE
  if(byPass) 
{
 # print(head(SS$GWAS))
      ps=SS$TV$ps
      nobs=SS$TV$nobs
      maf=GWAS$maf
  #maf=SS$TV$maf
      rsquare_base=SS$TV$rsquare_base
      rsquare=SS$TV$rsquare
      df=SS$TV$df
      tvalue=SS$TV$tvalue
      stderr=SS$TV$stderr
      effect.est=SS$mc
      effect=SS$mc
      #GI=cbind(GI,effect)
     
   if(DP$file.output & !is.null(SS$Compression) & !is.na(SS$Compression[1,6])){
     GAPIT.Compression.Visualization(Compression = SS$Compression, 
                                     name.of.trait = DP$name.of.trait,
                                     file.output = file.output)
   }
  
}else{
  maf=GI$maf
  ps=GI$P.value
  nobs=GI$nobs
  rsquare_base=rep(NA,length(ps))
  rsquare=rep(NA,length(ps))
  df=rep(NA,length(ps))
  tvalue=rep(NA,length(ps))
  stderr=rep(NA,length(ps))
  effect.est=GI$effect
  

  }
  if(is.na(maf[1]))  maf=matrix(.5,nrow(GI),1)
if(!is.null(IC$GD)&DP$SNP.test)
{ 
  
  print("Filtering SNPs with MAF..." )
	#index=maf>=DP$SNP.MAF	
	#PWI.Filtered=cbind(GI[,-5],rsquare_base,rsquare)
  PWI.Filtered=cbind(GWAS[,1:6],rsquare_base,rsquare)
	colnames(PWI.Filtered)=c("SNP","Chromosome","Position ","P.value", "maf", "nobs", "Rsquare.of.Model.without.SNP","Rsquare.of.Model.with.SNP")
  
  # if(!is.null(PWI.Filtered))
  # {
  #Run the BH multiple correction procedure of the results
  #Create PWIP, which is a table of SNP Names, Chromosome, bp Position, Raw P-values, FDR Adjusted P-values
  print("Calculating FDR..." )
  PWIP <- GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure(PWI = PWI.Filtered, FDR.Rate = FDR.Rate, FDR.Procedure = "BH")
  #print(str(PWIP))  
  print("QQ plot..." )
   if(DP$file.output) GAPIT.QQ(P.values = GI$P.value, name.of.trait = DP$name.of.trait,DPP=DP$DPP)
   print("Manhattan plot (Genomewise)..." )
   if(DP$file.output) GAPIT.Manhattan(GI.MP = GWAS[,2:4], name.of.trait = DP$name.of.trait, DPP=DP$DPP, plot.type = "Genomewise",cutOff=DP$cutOff,seqQTN=DP$QTN.position,plot.style=DP$plot.style,plot.bin=DP$plot.bin,chor_taxa=DP$chor_taxa)
   #print("@@@@@@@@@@@@@@@@@@@@@@@@")
   print("Manhattan plot (Chromosomewise)..." )

 if(DP$file.output) GAPIT.Manhattan(GI.MP = GWAS[,2:4],GD=IC$GD[,-1], CG=DP$CG,name.of.trait = DP$name.of.trait, DPP=DP$DPP, plot.type = "Chromosomewise",cutOff=DP$cutOff,plot.bin=DP$plot.bin)

  #Association Table
  print("Association table..." )
  
  print("Joining tvalue and stderr" )
  # print(head(DP$chor_taxa))
  # print(head(as.numeric(as.matrix(GWAS[,2]))))
  # print(all.equal(as.character(DP$chor_taxa),as.character(unique(as.numeric(as.matrix(GWAS[,2]))))))
  # print(head(GWAS))
   if(all.equal(as.character(DP$chor_taxa),as.character(unique(sort(as.numeric(as.matrix(GWAS[,2]))))))!=TRUE)
   { 
     chro=as.numeric(as.matrix(GWAS[,2]))
     chor_char=unique(DP$chor_taxa)
     # print(chro)
     # print(chor_char)
     for(i in 1:length(unique(chro)))
     {
      chro[chro==i]=chor_char[i]
     }
     GWAS[,2]=chro
   }
   #print(head(GWAS))
   #print(head(DP$GM))
  # print(length(tvalue))
  # print(length(stderr))
  # print(length(effect.est))
        DTS=cbind(GWAS[,1:3],df,tvalue,stderr,effect.est)
        colnames(DTS)=c("SNP","Chromosome","Position","DF","t Value","std Error","effect")	

  print("Creating ROC table and plot" )
  # print(head(tvalue))
  # print(head(stderr))
  # print(head(var(as.matrix(IC$Y[,2]))))
if(DP$file.output) myROC=GAPIT.ROC(t=as.numeric(tvalue),se=as.numeric(stderr),Vp=stats::var(as.matrix(IC$Y[,2])),trait=DP$name.of.trait)
  print("ROC table and plot created" )

  print("MAF plot..." )
if(DP$file.output) myMAF1=GAPIT.MAF(MAF=maf,P=ps,E=NULL,trait=DP$name.of.trait)

print("GAPIT.Interactive.Manhattan")
print(DP$Inter.type)
#GI=GI[order(GI[,4]),]
#print(head(GI))
if(DP$file.output){

   if(ncol(GI)>1)
     {
      new_GI=merge(PWIP$PWIP,GI[,c("SNP","effect")],by.x="SNP",by.y="SNP")
     }else{
      new_GI=GI
     }
   new_GI=new_GI[order(new_GI[,4]),]
# }
# print("@@@@@@")
if(DP$Inter.Plot) 
  {
  print(DP$Inter.type)
  GAPIT.Interactive.Manhattan(GWAS=new_GI,X_fre=maf,plot.type=DP$Inter.type,name.of.trait = DP$name.of.trait)
   } 
# if(DP$file.output){
   
   # print(as.character(DP$chor_taxa))
   # print(as.character(unique(as.numeric(as.matrix(new_GI[,2])))))
   if(all.equal(as.character(DP$chor_taxa),as.character(sort(unique(as.numeric(as.matrix(new_GI[,2]))))))!=TRUE)
   {
     # print("@@@")
     chro=as.numeric(as.matrix(new_GI[,2]))
     chor_char=unique(DP$chor_taxa)

     for(i in 1:length(unique(chro)))
     {
      chro[chro==i]=chor_char[i]
     }
     new_GI[,2]=chro
   }

   utils::write.table(new_GI, paste("GAPIT.", DP$name.of.trait, ".GWAS.Results.csv", sep = ""), quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
   utils::write.table(DTS, paste("GAPIT.", DP$name.of.trait, ".Df.tValue.StdErr.csv", sep = ""), quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
   #print(head(GWAS.2))
   #if(byPass) write.table(GWAS.2[,1:4], paste("GAPIT.", DP$name.of.trait, ".Allelic_Effect_Estimates.csv", sep = ""), quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
     }#end file.output
  # }#PWI.Filtered
}#end IC$GD)
  print("GAPIT.ID accomplished successfully for multiple traits. Results are saved")
  return(invisible(NULL))
}#is.null(DP)&is.null(IC)

}  #end of GAPIT.ID function
#=============================================================================================

