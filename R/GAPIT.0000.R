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

`GAPIT.0000` <-
function(){
##############################################################################################
#GAPIT: Genome Association and Prediction Integrated Tool
#Objective 1: State of art methods for high  power, accuracy and speed;
#Objective 2: User friendly by design, help documents, and web forum;
#Objective 3: Comprehensive output to interpret data and results;
#Objective 4: Informative tables and high quality figures for reports and publication;

#Methods implimented: 
# 1. GLM (Structure or Q method for GWAS, Pritchard et. al. Genetics, 2000)
# 2. MLM (Q+K, Yu et. al. Nature Genetics, 2006)
# 3. gBLUP (Marker based kinship, Zhang et. al. Journal of Animal Science, 2007)
# 4. PCA (Zhao et. al. Plos Genetics, 2007)
# 5. EMMA (Kang et. al. Genetics, 2008)
# 6. CMLM (Zhang et. al. Nature Genetics, 2010)
# 7. EMMAx (Kang et. al. Nature Genetics, 2010)
# 8. P3D (Zhang et. al. Nature Genetics, 2010)
# 9. FaST-LMM (Lippert et. al. Nature Methods, 2011)
# 10. ECMLM (Li et. al. BMC Bioogy, 2014)
# 11. SUPER (Wang et. al. PLoS One, 2014)

GAPIT.Version="2025.03.03, GAPIT 4.0"
return(GAPIT.Version)
}
#=============================================================================================

