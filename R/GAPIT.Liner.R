`GAPIT.Liner` <-
function(Y,GD,CV){
    #Object: To have Y, GD and CV the same size and order
    #Input: Y,GDP,GM,CV
    #Input: GD - n by m +1 dataframe or n by m big.matrix
    #Input: GDP - n by m matrix. This is Genotype Data Pure (GDP). THERE IS NOT COLUMN FOR TAXA. 
    #Input: orientation-Marker in GDP go colmun or row wise
    #Requirement: Y, GDP and CV have same taxa order. GDP and GM have the same order on SNP
    #Output: GWAS,GPS,Pred
    #Authors: Zhiwu Zhang
    # Last update: Febuary 24, 2013
    ##############################################################################################
    #print("GAPIT.Liner Started")
    #print(date())
    #print("Memory used at begining of BUS")
    #print(memory.size())
    #print("dimension of Y,GD and CV at begining")
    #print(dim(Y))
    #print(dim(GD))
    #print(dim(CV))
    
    if(!is.null(CV))taxa=intersect(intersect(GD[,1],Y[,1]),CV[,1])
    if(is.null(CV))taxa=intersect(GD[,1],Y[,1])

    Y=Y[match(taxa, Y[,1], nomatch = 0),]
    GD=GD[match(taxa, GD[,1], nomatch = 0),]
    
    if(!is.null(CV)) CV=CV[match(taxa, CV[,1], nomatch = 0),]
    Y = Y[order(Y[,1]),]
    GD = GD[order(GD[,1]),]
    if(!is.null(CV)) CV = CV[order(CV[,1]),]

    #print("dimension of Y,GD and CV at end")
    #print(dim(Y))
    #print(dim(GD))
    #print(dim(CV))
    
  print("GAPIT.Liner accomplished successfully")
  return (list(Y=Y,GD=GD,CV=CV))
}#The function GAPIT.Liner ends here
#=============================================================================================

