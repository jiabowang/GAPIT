`GAPIT.Table` <-
function(final.table = final.table, name.of.trait = name.of.trait,SNP.FDR=1){
#Object: Make and export a table of summary information from GWAS
#Output: A table summarizing GWAS results
#Authors: Alex Lipka and Zhiwu Zhang
# Last update: May 10, 2011 
##############################################################################################

#Filter SNPs by FDR
index=(final.table[,7]<=SNP.FDR)
final.table=final.table[index,]

#Export this summary table as an excel file
write.table(final.table, paste("GAPIT.", name.of.trait, ".GWAS.Results.csv", sep = ""), quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)


#print("GAPIT.Table accomplished successfully!")
  

}   #GAPIT.Table ends here
#=============================================================================================

