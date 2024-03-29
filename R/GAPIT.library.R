if(!require(gplots)) install.packages("gplots")
if(!require(genetics)) install.packages("genetics")
if(!require(ape)) install.packages("ape")
if(!require(compiler)) install.packages("compiler")
if(!require(grid)) install.packages("grid")
if(!require(bigmemory)) install.packages("bigmemory")
if(!require(EMMREML)) install.packages("EMMREML")
if(!require(scatterplot3d)) install.packages("scatterplot3d")
if(!require(lme4)) install.packages("lme4")
# if(!require(rgl)) install.packages("rgl")

if(!'multtest'%in% installed.packages()[,"Package"]){
	if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
   BiocManager::install("multtest")
   BiocManager::install("snpStats")
}
