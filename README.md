GAPIT3
===========
Genome Association Predict Integrate Tools

![GAPIT](https://github.com/jiabowang/GAPIT3/blob/master/material/LOGO_WEB.png)

Authors: 
-----

Jiabo Wang and Zhiwu Zhang

Contact:
-----

wangjiaboyifeng@163.com (Jiabo)

Source:
-----

[User manual](http://zzlab.net/GAPIT/gapit_help_document.pdf)

[Demo Data](http://zzlab.net/GAPIT/GAPIT_Tutorial_Data.zip)

[Source code](http://zzlab.net/GAPIT/gapit_functions.txt)

Contents:
-----


Start
======

   GAPIT is a package that is run in the R software environment, which can be freely downloaded from [http://www.r-project.org](http://www.r-project.org) or [http://www.rstudio.com](http://www.rstudio.com).

Function Loading
-------
  Now GAPIT can load library by only one funciton. 
      
    source("www.zzlab.net/GAPIT/GAPIT.library.R")
  
  After loading library, we need to source GAPIT function.
  
    source("www.zzlab.net/GAPIT/gapit_functions.txt")
  

Data Preparing
-------


### Phenotype Data

The user has the option of performing GWAS on multiple phenotypes in GAPIT. This is achieved by including all phenotypes in the text file of phenotypic data. Taxa names should be in the first column of the phenotypic data file and the remaining columns should contain the observed phenotype from each individual. Missing data should be indicated by either “NaN” or “NA”. 

![Phenotype file](https://github.com/jiabowang/GAPIT3/blob/master/material/phenotype.png)

### Genotype Data

#### Hapmap Format
Hapmap is a commonly used format for storing sequence data where SNP information is stored in the rows and taxa information is stored in the columns. This format allows the SNP information (chromosome and position) and genotype of each taxa to be stored in one file.


#### Numeric Format
GAPIT also accepts the numeric format. Homozygotes are denoted by “0” and “2” and heterozygotes are denoted by “1” in the “GD” file.  Any numeric value between “0” and “2” can represent imputed SNP genotypes. The first row is a header file with SNP names, and the first column is the taxa name.
The “GM” file contains the name and location of each SNP. The first column is the SNP id, the second column is the chromosome, and the third column is the base pair position. As seen in the example, the first row is a header file.

Anlysis
======
GWAS
-----

* GLM

The GAPIT use Least Squares to solve the modle. The code of GAPIT running GLM is:
      
      myGAPIT_GLM <- GAPIT(
      Y=myY[,c(1,2)],
      GD=myGD,
      GM=myGM,
      model="GLM",
      PCA.total=5,
      file.output=T
      )


* MLM

EMMA method is used in GAPIT, the code of MLM is:

      myGAPIT_MLM <- GAPIT(
      Y=myY[,c(1,2)],
      GD=myGD,
      GM=myGM,
      model="MLM",
      PCA.total=5,
      file.output=T
      )


* CMLM

Compress Mixed Linear Model is published by Zhang in 2010. The code of CMLM is:

      myGAPIT_CMLM <- GAPIT(
      Y=myY[,c(1,2)],
      GD=myGD,
      GM=myGM,
      model="CMLM",
      PCA.total=5,
      file.output=T
      )

* MLMM

Multiple Loci Mixied linear Model is published by Segura in 2012. The code of MLMM is:

      myGAPIT_MLMM <- GAPIT(
      Y=myY[,c(1,2)],
      GD=myGD,
      GM=myGM,
      model="MLMM",
      PCA.total=5,
      file.output=T
      )

* SUPER

Settlement of MLM Under Progressively Exclusive Relation- ship is published by Qishan in 2014. The code of MLMM is:

      myGAPIT_SUPER <- GAPIT(
      Y=myY[,c(1,2)],
      GD=myGD,
      GM=myGM,
      model="SUPER",
      PCA.total=5,
      file.output=T
      )


* Farm-CPU

Fixed and random model Circulating Probability Unification (FarmCPU) is published by Xiaolei in 2016. The code of MLMM is:

      myGAPIT_FarmCPU <- GAPIT(
      Y=myY[,c(1,2)],
      GD=myGD,
      GM=myGM,
      model="FarmCPU",
      PCA.total=5,
      file.output=T
      )



GS
-----

* gBLUP

gBLUP used marker kinship to replace the pedgree relationship matrix. The code is:

      myGAPIT_gBLUP <- GAPIT(
      Y=myY[,c(1,2)],
      GD=myGD,
      GM=myGM,
      group.from=nrow(myGD),
      group.to=nrow(myGD),
      PCA.total=5,
      file.output=T
      )



* cBLUP

cBLUP used group kinship to replace the individual matrix. The code is:

      myGAPIT_cBLUP <- GAPIT(
      Y=myY[,c(1,2)],
      GD=myGD,
      GM=myGM,
      group.from=10,
      group.to=nrow(myGD),
      group.by=10,
      PCA.total=5,
      file.output=T
      )

* sBLUP

sBLUP used SUPER method to build psedue QTN kinship matrix. The code is:

      myGAPIT_sBLUP <- GAPIT(
      Y=myY[,c(1,2)],
      GD=myGD,
      GM=myGM,
      sangwich.top="GLM",sangwich.bottom="SUPER",
      LD=0.1,
      SUPER_GS=T,
      PCA.total=5,
      file.output=T
      )


Result
=====

Example
=====

Citation
=====






