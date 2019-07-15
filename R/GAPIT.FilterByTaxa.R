`GAPIT.FilterByTaxa` <-
function(taxa,Data){
    #Object: To filter a data (Y, CV or GD) by taxa
    #Input: taxa - vector of taxa
    #Input: data - data frame with first column as taxa
    #Requirement: all taxa must be in data
    #Output: filtered data
    #Authors: Zhiwu Zhang
    # Last update: May 22, 2013
##############################################################################################
   #print("GAPIT.FilterByTaxa Started")

    Data=Data[match(taxa, Data[,1], nomatch = 0),]

  return (Data)

}#The function GAPIT.FilterByTaxa ends here
#=============================================================================================

