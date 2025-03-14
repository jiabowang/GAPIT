
GAPIT [![](https://img.shields.io/badge/Issues-0%2B-brightgreen.svg)](https://github.com/jiabowang/GAPIT/issues)
===========

Genome Association and Prediction Integrated Tools Version III

![The GAPIT logo](man/figures/LOGO_WEB.png "Genome Association and Prediction Integrated Tools logo")


Citation
=====

The DOI of GAPIT GitHub repository is:
DOI: 10.5281/zenodo.7931838

If you use GAPIT and publish your analysis, please report the program version and cite the appropriate article:

  
The citation for GAPIT3 is:    
Wang J., Zhang Z., GAPIT Version 3: Boosting Power and Accuracy for Genomic Association and Prediction, Genomics, Proteomics & Bioinformatics (2021), doi: https://doi.org/10.1016/j.gpb.2021.08.005.

The citation for GAPIT2 is:    
Tang Y., Liu X., Wang J., Li M., Wang Q., et al., 2016 GAPIT Version 2: An Enhanced Integrated Tool for Genomic Association and Prediction. Plant J. 9, https://10.3835/plantgenome2015.11.0120.
  
The citation for GAPIT is:    
Lipka A. E., Tian F., Wang Q., Peiffer J., Li M., et al., 2012 GAPIT: genome association and prediction integrated tool. Bioinformatics 28: 2397–2399, https://doi.org/10.1093/bioinformatics/bts444.

The citation for cBLUP and sBLUP is:    
Wang J., Zhou Z., Zhang Z., Li H., Liu D., et al., 2018 Expanding the BLUP alphabet for genomic prediction adaptable to the genetic architectures of complex traits. Heredity https://doi.org/10.1038/s41437-018-0075-0.

The citation for BLINK is:    
Huang M, Liu X, Zhou Y, Summers RM, Zhang Z. BLINK: A package for the next level of genome-wide association studies with both individuals and markers in the millions. Gigascience. https://doi.org/10.1093/gigascience/giy154.

The citation for Farm-CPU is:    
Liu X., Huang M., Fan B., Buckler E. S., Zhang Z., 2016 Iterative Usage of Fixed and Random Effect Models for Powerful and Efficient Genome-Wide Association Studies. PLoS Genet. 12: e1005767. https://doi.org/10.1371/journal.pgen.1005767.

The citation for SUPER method is:    
Wang Q., Tian F., Pan Y., Buckler E. S., Zhang Z., 2014 A SUPER Powerful Method for Genome Wide Association Study (Y Li, Ed.). PLoS One 9: e107684, https://doi.org/10.1371/journal.pone.0107684.

The citation for P3D is:    
Zhang Z., Ersoz E., Lai C. Q., Todhunter R. J., Tiwari H. K., et al., 2010 Mixed linear model approach adapted for genome-wide association studies. Nat. Genet. 42: 355–360. https://doi.org/10.1038/ng.546.


Authors: 
-----

Jiabo Wang and Zhiwu Zhang

Contact:
-----

wangjiaboyifeng@163.com (Jiabo)

Source:
-----

[User manual](https://github.com/jiabowang/GAPIT/raw/refs/heads/master/Documents/gapit_help_document.pdf)

[Demo Data](https://github.com/jiabowang/GAPIT/raw/refs/heads/master/Documents/GAPIT_Tutorial_Data.zip)

[Source code](https://raw.githubusercontent.com/jiabowang/GAPIT/refs/heads/master/gapit_functions.txt)

Contents:
-----

   
   * [Start](#start)
      * [Installing GAPIT from source functions](#installation-from-source-functions)
      * [Installing GAPIT from GitHub](#installation-from-github)
      * [Installing GAPIT from archive](#installation-from-an-archive)
   * [Analysis](#analysis)
      * [GWAS](#gwas)
      * [GS](#gs)
   * [Example](#example)


Start
======


GAPIT is a package that is run in the R software environment.
R can be freely downloaded from [http://www.r-project.org](http://www.r-project.org).
We also recommend the integrated development environment RStudio which is also freely available at [http://www.rstudio.com](http://www.rstudio.com).


GAPIT can currently be installed in several ways.

- From source on the internet
- From GitHub
- From an archive


## Installation from source functions

GAPIT can be loaded with a single funciton. 


```
R> source("http://zzlab.net/GAPIT/gapit_functions.txt")
```

Or from GitHub function.
```
R> source("https://raw.githubusercontent.com/jiabowang/GAPIT/refs/heads/master/gapit_functions.txt", encoding = "UTF-8")
```

## Installation from GitHub

Installation can also be made from GitHub when the R package `devtools` is available.

```
R> install.packages("devtools")
R> devtools::install_github("jiabowang/GAPIT", force=TRUE)
R> library(GAPIT)

or
R> install.packages("remotes")
R> remotes::install_github("jiabowang/GAPIT")
R> library(GAPIT)


```

## Installation from an archive


GAPIT can be installed from an archive such as \*.tar.gz or \*.zip archive.
An archive can be downloaded from the "releases" page.
If you would like the latest version of GAPIT from the GitHub site you may want to clone it and then build it (this may require Rtools on Windows).

```
bash$ git clone git@github.com:jiabowang/GAPIT.git
bash$ R CMD build GAPIT
```

Once an archive has been obtained it can be installed from a shell, similar to as follows.


```
bash$ R CMD INSTALL GAPIT_3.5.0.9000.tar.gz
```

Or similarly from within R.

```
R> install.packages("GAPIT_3.5.0.9000.tar.gz", repos = NULL, type="source")
```

## In some case of the BiocManager can not be installed

Installation of same packages such as multtest and biobase can not be intalled from BiocManager. These packages can be downloaded in the Bioconductor website and be installed from local source files.

The website of Bioconductor is here:

```
https://bioconductor.org/packages/3.19/bioc/
```

Analysis
======

The taxa and order of individual among genotype file and phenotype files could be different. The GAPIT can automaticly filter and order the common taxa among among genotype file and phenotype files.

GWAS
-----
GWAS methods include: General Linear Model(GLM),Mixed Linear Model(MLM),compressed Mixed Linear Model(CMLM),SUPER, FarmCPU, MLMM, and BLINK. Users can use model="GLM" to select method. 

GS
-----
GS methods include: gBLUP, cBLUP, sBLUP, and GAGBLUP. Besides GAGBLUP, all GS methods could be seleted by using model="gBLUP"

   
Example
=====


```
# loading packages for GAPIT and GAPIT functions
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")
# loading data set
myY=read.table(file="https://github.com/jiabowang/GAPIT/raw/refs/heads/master/Documents/mdp_traits.txt", head = TRUE)
myGD=read.table("https://github.com/jiabowang/GAPIT/raw/refs/heads/master/Documents/mdp_numeric.txt",head=T)
myGM=read.table("https://github.com/jiabowang/GAPIT/raw/refs/heads/master/Documents/mdp_SNP_information.txt",head=T)
#myG=read.table(file="https://github.com/jiabowang/GAPIT/raw/refs/heads/master/Documents/mdp_genotype_test.hmp.txt", head = FALSE)
# performing simulation phenotype
set.seed(198521)
mysimulation<-GAPIT(h2=0.7,NQTN=20,GD=myGD,GM=myGM)
myY=mysimulation$Y


myGAPIT <- GAPIT(
  Y=myY[,c(1,2)],
  GD=myGD,
  GM=myGM,
  model=c("GLM","MLM","SUPER","MLMM","FarmCPU","Blink"),# choose model
  #model=c("FarmCPU"),
  PCA.total=3,                                          # set total PCAs
  NJtree.group=4,                                       # set the number of clusting group in Njtree plot
  QTN.position=mysimulation$QTN.position,
  Inter.Plot=TRUE,                                      # perform interactive plot
  Multiple_analysis=TRUE,                               # perform multiple analysis
  PCA.3d=TRUE,                                          # plot 3d interactive PCA
  file.output=T
)
```
More details please check and ref user manual.

