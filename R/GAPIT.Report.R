`GAPIT.Report` <-
function(name.of.trait=NULL,GWAS=NULL,pred=pred,ypred=NULL,tvalue=NULL,stderr=NULL,Vp=1,
DPP=100000,cutOff=.01,threshold.output=.01,MAF=NULL,seqQTN=NULL,MAF.calculate=FALSE,plot.style="rainbow"){
#Object: Out put plots and tables
#Input: GWAS,name.of.trait, DPP 
#Requirement: None
#Output: Graphs and tables
#Output: return ycor if ypred is not null
#Authors: Zhiwu Zhang
# Date  start: April 2, 2013
# Last update: April 2, 2013
##############################################################################################
#print("GAPIT.Report Started")
#print(seqQTN)
#Manhattan Plots
#print("Manhattan plot (Genomewise)..." )
if(plot.style=="FarmCPU"){
    GAPIT.Manhattan(GI.MP = GWAS[,2:4], name.of.trait = name.of.trait, DPP=DPP, plot.type = "Genomewise",cutOff=cutOff,seqQTN=seqQTN,plot.style=plot.style)
}
if(plot.style=="rainbow"){
GAPIT.Manhattan(GI.MP = GWAS[,2:4], name.of.trait = name.of.trait, DPP=DPP, plot.type = "Genomewise",cutOff=cutOff,seqQTN=seqQTN,plot.style=plot.style)
    #}
#print("Manhattan plot (Chromosomewise)..." )
GAPIT.Manhattan(GI.MP = GWAS[,2:4], name.of.trait = name.of.trait, DPP=DPP, plot.type = "Chromosomewise",cutOff=cutOff,plot.style=plot.style)
}


#QQ plots
#print("QQ plotting..." )
#if(plot.style=="rainbow"){
#    GAPIT.QQ(P.values = GWAS[,4], name.of.trait = name.of.trait,DPP=DPP)
#}
#if(plot.style=="nature"){
GAPIT.QQ(P.values = GWAS[,4], name.of.trait = name.of.trait,DPP=DPP,plot.style=plot.style)
    #}
#Association Table
#print("Create association table..." )
index=1:nrow(GWAS)
if(threshold.output<1)index=which(GWAS[,4]<threshold.output)
if(plot.style=="FarmCPU"){
write.table(GWAS[index,], paste("FarmCPU.", name.of.trait, ".GWAS.Results.csv", sep = ""), quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
}
if(plot.style=="rainbow"){
write.table(GWAS[index,], paste("GAPIT.", name.of.trait, ".GWAS.Results.csv", sep = ""), quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
}

#Prediction
#print("Create prediction table..." )
#if(!is.null(pred)) write.table(pred, paste("GAPIT.", name.of.trait, ".Pred.csv", sep = ""), quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
#print("Create prediction table for unknown phenotype...")
#if(!is.null(ypred)) write.table(ypred, paste("GAPIT.", name.of.trait, ".unknownY.Pred.csv", sep = ""), quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
if(!is.null(pred) || !is.null(ypred)){
myPred=FarmCPU.Pred(pred=pred,ypred=ypred,name.of.trait=name.of.trait)
}

#ROC
#print("Creating ROC table and plot" )
myROC=GAPIT.ROC(t=tvalue,se=stderr,Vp=Vp,trait=name.of.trait,plot.style=plot.style)

#MAF
#print("Creating MAF table and plot" )
if(MAF.calculate){
    myMAF=GAPIT.MAF(MAF=MAF,P=GWAS[,4],E=NULL,trait=name.of.trait,threshold.output=threshold.output,plot.style=plot.style)}
#print("Report accomplished" )
}#The function GAPIT.Report ends here
#=============================================================================================
