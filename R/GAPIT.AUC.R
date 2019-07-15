 #Object: To calculate Area Under (ROC) Curve (AUC)
 #Straitegy: NA
 #Output: P value
 #intput: beta-power and alpha-fdr or type I error
 #Authors: Zhiwu Zhang
 #Last update: December 18, 2015
##############################################################################################
GAPIT.AUC=function(beta=NULL,alpha=NULL){
	n=length(beta)
	#plot(alpha,beta,type="b")
	db=beta[-1]-beta[-n]
	da=1-.5*(alpha[-1]+alpha[-n])
	ab=da*db
	AUC=sum(ab)
	return(AUC)
}
#=============================================================================================
