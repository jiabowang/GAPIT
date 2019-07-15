
if(!require(gplots)) install.packages("gplots")
if(!require(LDheatmap)) install.packages("LDheatmap")
if(!require(genetics)) install.packages("genetics")
if(!require(ape)) install.packages("ape")
if(!require(compiler)) install.packages("compiler")

if(!require(EMMREML)) install.packages("EMMREML")
if(!require(scatterplot3d)) install.packages("scatterplot3d")

if(!'multtest'%in% installed.packages()[,"Package"]){
	source("http://www.bioconductor.org/biocLite.R")
	biocLite("multtest")
	biocLite("snpStats")
}
