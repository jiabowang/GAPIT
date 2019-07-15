`GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure` <-
function(PWI = PWI, FDR.Rate = 0.05, FDR.Procedure = "BH"){
#Object: Conduct the Benjamini-Hochberg FDR-Controlling Procedure
#Output: PWIP, number.of.significant.SNPs
#Authors: Alex Lipka and Zhiwu Zhang 
# Last update: May 5, 2011 
##############################################################################################
#Make sure that your compouter has the latest version of Bioconductor (the "Biobase" package) and multtest

if(is.null(PWI))
{
PWIP=NULL
number.of.significant.SNPs = 0
}

if(!is.null(PWI))
{  
 
    #library(multtest)
    
    if(dim(PWI)[1] == 1){
     PWIP <- cbind(PWI, PWI[4])
     colnames(PWIP)[9] <- "FDR_Adjusted_P-values"
    }
   
    if(dim(PWI)[1] > 1){ 
    #mt.rawp2adjp Performs the Simes procedure.  The output should be two columns, Left column: originial p-value
    #Right column: Simes corrected p-value
    res <- mt.rawp2adjp(PWI[,4], FDR.Procedure)

    #This command should order the p-values in the order of the SNPs in the data set
  adjp <- res$adjp[order(res$index), ]

  #round(adjp[1:7,],4)
    #Logical statment: 0, if Ho is not rejected; 1, if  Ho is rejected, by the Simes corrected p-value
#  temp <- mt.reject(adjp[,2], FDR.Rate)

    #Lists all number of SNPs that were rejected by the BY procedure
  #temp$r

    #Attach the FDR adjusted p-values to AS_Results

  PWIP <- cbind(PWI, adjp[,2])

    #Sort these data by lowest to highest FDR adjusted p-value
  PWIP <- PWIP[order(PWIP[,4]),]
  
  colnames(PWIP)[9] <- "FDR_Adjusted_P-values"
#  number.of.significant.SNPs = temp$r
  }
  #print("GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure accomplished successfully!")
}  
  #return(list(PWIP=PWIP, number.of.significant.SNPs = number.of.significant.SNPs))
  return(list(PWIP=PWIP))
}#GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure ends here
#=============================================================================================

