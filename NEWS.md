
#  3.3.0
Released on CRAN 202X-XX-XX


Addressed R CMD check GAPIT_3.5.0.tar.gz issues 2025-10

Note: break used in wrong context: no loop is visible at GAPIT.R:531
Note: possible error in 'GAPIT(Y = Y, CV = CV, ': unused argument (GTindex = GTindex) at GAPIT.Bus.R:24 
Note: possible error in 'layout(title = "Interactive.Manhattan.Plot", ': unused arguments (title = "Interactive.Manhattan.Plot", xaxis = list(title = "Chromsome", zeroline = FALSE, showticklabels = FALSE), yaxis = list(title = "-Log10(p)")) at GAPIT.Interactive.Manhattan.R:85 
Note: possible error in 'layout(title = "Interactive.Multiple_Synthesis.Manhattan.Plot", ': unused arguments (title = "Interactive.Multiple_Synthesis.Manhattan.Plot", xaxis = list(title = "Chromsome", zeroline = FALSE, showticklabels = FALSE), yaxis = list(title = "-Log10(p)")) at GAPIT.Multiple_Synthesis.R:244 
Note: possible error in 'layout(title = "Interactive.Multiple_Synthesis.Manhattan.Plot", ': unused arguments (title = "Interactive.Multiple_Synthesis.Manhattan.Plot", xaxis = list(title = "Chromsome", zeroline = FALSE, showticklabels = FALSE), yaxis = list(title = "-Log10(p)")) at GAPIT.Multiple_Synthesis.R:277 
Note: possible error in 'GAPIT(Y = myY, GD = myGD, ': unused arguments (threshold.output = 0.001, iteration.output = TRUE) at GAPIT.Power.compare_plink.R:123 
Note: possible error in 'GAPIT(Y = myY, GD = myGD, ': unused arguments (threshold.output = 0.001, iteration.output = TRUE) at GAPIT.Power.compare_plink.R:140 
Note: possible error in 'GAPIT(Y = myY, GD = myGD, ': unused arguments (LD = 0.1, threshold.output = 0.001, iteration.output = TRUE) at GAPIT.Power.compare_plink.R:159 
Note: possible error in 'GAPIT(Y = myY, GD = myGD, ': unused arguments (threshold.output = 0.001, iteration.output = TRUE) at GAPIT.Power.compare_plink.R:177 
Note: possible error in 'GAPIT(Y = myY, G = myG, ': unused arguments (threshold.output = 0.001, iteration.output = TRUE) at GAPIT.Power.compare_plink.R:196 

* checking code files for non-ASCII characters ... WARNING
Found the following files with non-ASCII characters:
  R/GAPIT.3D.PCA.python.R
  R/GAPIT.Main.R
Portable packages must use only ASCII characters in their R code and
NAMESPACE directives, except perhaps in comments.
Use \uxxxx escapes for other characters.
Function ‘tools::showNonASCIIfile’ can help in finding non-ASCII
characters in files.

* checking dependencies in R code ... WARNING
'library' or 'require' calls not declared from:
  ‘lme4’ ‘plotly’ ‘rgl’ ‘rglwidget’
'library' or 'require' calls in package code:
  ‘lme4’ ‘plotly’ ‘rgl’ ‘rglwidget’
  Please use :: or requireNamespace() instead.
  See section 'Suggested packages' in the 'Writing R Extensions' manual.
Namespaces in Imports field not imported from:
  ‘grid’ ‘lattice’ ‘rgl’
  All declared Imports should be used.

* checking dependencies in R code ... NOTE
Namespaces in Imports field not imported from:
  ‘grid’ ‘lattice’ ‘rgl’
  All declared Imports should be used.

* checking R code for possible problems ... NOTE
GAPIT.Bread: warning in GAPIT(Y = Y, CV = CV, Z = Z, KI = KI, GD = GD,
  GM = GM, model = ("GLM"), QC = FALSE, CV.Extragenetic =
  CV.Extragenetic, file.output = file.output): partial argument match
  of 'QC' to 'QC.Y'
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Bread.R:24-36)
GAPIT.Bread: warning in GAPIT(Y = Y, CV = CV, Z = Z, KI = KI, GD = GD,
  GM = GM, model = "MLM", QC = FALSE, CV.Extragenetic =
  CV.Extragenetic, file.output = file.output): partial argument match
  of 'QC' to 'QC.Y'
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Bread.R:49-61)
GAPIT.Bread: warning in GAPIT(Y = Y, CV = CV, Z = Z, KI = KI, GD = GD,
  GM = GM, model = "CMLM", QC = FALSE, CV.Extragenetic =
  CV.Extragenetic, file.output = file.output): partial argument match
  of 'QC' to 'QC.Y'
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Bread.R:73-85)

Added library_name:: prefix to function names.
* checking R code for possible problems ... NOTE
GAPIT.Bus: no visible binding for global variable ‘FDR.Rate’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Bus.R:124)
GAPIT.Bus: no visible binding for global variable ‘FDR.Rate’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Bus.R:435)
GAPIT.Cor.matrix: no visible global function definition for ‘cor’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Cor.matrix.R:6)
GAPIT.EMMAxP3D: no visible global function definition for ‘var’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.EMMAxP3D.R:838-843)
GAPIT.Genotype.View: no visible global function definition for ‘median’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Genotype.View.R:109)
GAPIT.Genotype.View: no visible global function definition for ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Genotype.View.R:121)
GAPIT.Genotype.View: no visible global function definition for ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Genotype.View.R:124)
GAPIT.Genotype.View: no visible global function definition for ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Genotype.View.R:125)
GAPIT.Genotype.View: no visible global function definition for ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Genotype.View.R:126)
GAPIT.Genotype.View: no visible global function definition for ‘head’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Genotype.View.R:130)
GAPIT.Genotype.View: no visible global function definition for ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Genotype.View.R:133)
GAPIT.Genotype.View: no visible global function definition for ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Genotype.View.R:134)
GAPIT.Genotype.View: no visible global function definition for ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Genotype.View.R:135)
GAPIT.Genotype.View: no visible global function definition for ‘hist’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Genotype.View.R:138)
GAPIT.Genotype.View: no visible global function definition for ‘abline’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Genotype.View.R:159)
GAPIT.Genotype.View: no visible global function definition for ‘lines’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Genotype.View.R:204)
GAPIT.Genotype.View: no visible global function definition for
  ‘write.csv’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Genotype.View.R:209)
GAPIT.Genotype.View: no visible global function definition for ‘layout’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Genotype.View.R:235-237)
GAPIT.Genotype.View: no visible global function definition for ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Genotype.View.R:238)
GAPIT.Genotype.View: no visible global function definition for ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Genotype.View.R:243)
GAPIT.Genotype.View: no visible global function definition for ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Genotype.View.R:244)
GAPIT.Genotype.View: no visible global function definition for
  ‘write.csv’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Genotype.View.R:249)
GAPIT.Genotype.View: no visible global function definition for ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Genotype.View.R:252)
GAPIT.Genotype.View: no visible global function definition for ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Genotype.View.R:253)
GAPIT.Genotype.View: no visible global function definition for ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Genotype.View.R:256)
GAPIT.Genotype.View: no visible global function definition for ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Genotype.View.R:257)
GAPIT.Genotype.View: no visible global function definition for ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Genotype.View.R:258)
GAPIT.IC: no visible global function definition for ‘rnorm’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.IC.R:115)
GAPIT.LD.decay: no visible binding for global variable ‘myGM’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.LD.decay.R:43)
GAPIT.LD.decay: no visible global function definition for ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.LD.decay.R:80)
GAPIT.LD.decay: no visible global function definition for ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.LD.decay.R:86)
GAPIT.LD.decay: no visible global function definition for ‘abline’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.LD.decay.R:88)
GAPIT.LD.decay: no visible global function definition for ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.LD.decay.R:89)
GAPIT.LD.decay: no visible global function definition for ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.LD.decay.R:90)
GAPIT.LD.decay: no visible global function definition for ‘title’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.LD.decay.R:93)
GAPIT.LD.decay: no visible global function definition for ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.LD.decay.R:104)
GAPIT.LD.decay: no visible global function definition for ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.LD.decay.R:109)
GAPIT.LD.decay: no visible global function definition for ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.LD.decay.R:110)
GAPIT.LD.decay: no visible global function definition for ‘title’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.LD.decay.R:113)
GAPIT.LD.decay: no visible global function definition for ‘lines’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.LD.decay.R:132)
GAPIT.LD.decay: no visible global function definition for ‘hist’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.LD.decay.R:135)
GAPIT.LD.decay: no visible global function definition for ‘hist’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.LD.decay.R:140)
GAPIT.LD.decay: no visible global function definition for ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.LD.decay.R:147)
GAPIT.LD.decay: no visible global function definition for ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.LD.decay.R:150)
GAPIT.LD.decay: no visible global function definition for ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.LD.decay.R:151)
GAPIT.LD.decay: no visible global function definition for ‘mtext’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.LD.decay.R:153)
GAPIT.LD.decay: no visible global function definition for ‘legend’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.LD.decay.R:156-158)
GAPIT.LD.decay: no visible global function definition for ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.LD.decay.R:156-158)
GAPIT.Manhattan: no visible global function definition for ‘pdf’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:97)
GAPIT.Manhattan: no visible global function definition for ‘layout’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:99)
GAPIT.Manhattan: no visible global function definition for ‘cor’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:137)
GAPIT.Manhattan: no visible global function definition for
  ‘colorRampPalette’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:160)
GAPIT.Manhattan: no visible global function definition for ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:175)
GAPIT.Manhattan: no visible global function definition for ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:176)
GAPIT.Manhattan: no visible global function definition for ‘abline’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:183)
GAPIT.Manhattan: no visible global function definition for ‘abline’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:184)
GAPIT.Manhattan: no visible global function definition for ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:185)
GAPIT.Manhattan: no visible global function definition for ‘barplot’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:187)
GAPIT.Manhattan: no visible global function definition for ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:189)
GAPIT.Manhattan: no visible global function definition for ‘dev.off’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:192)
GAPIT.Manhattan: no visible global function definition for ‘rainbow’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:218)
GAPIT.Manhattan: no visible global function definition for ‘pdf’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:294)
GAPIT.Manhattan: no visible global function definition for ‘pdf’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:296)
GAPIT.Manhattan: no visible global function definition for ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:298)
GAPIT.Manhattan: no visible global function definition for ‘abline’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:304)
GAPIT.Manhattan: no visible global function definition for ‘points’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:305)
GAPIT.Manhattan: no visible global function definition for ‘points’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:306)
GAPIT.Manhattan: no visible global function definition for ‘abline’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:309)
GAPIT.Manhattan: no visible global function definition for ‘points’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:310)
GAPIT.Manhattan: no visible global function definition for ‘points’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:311)
GAPIT.Manhattan: no visible global function definition for ‘abline’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:316)
GAPIT.Manhattan: no visible global function definition for ‘abline’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:318)
GAPIT.Manhattan: no visible global function definition for ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:326)
GAPIT.Manhattan: no visible global function definition for ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:327)
GAPIT.Manhattan: no visible global function definition for ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:328)
GAPIT.Manhattan: no visible global function definition for ‘box’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:330)
GAPIT.Manhattan: no visible global function definition for ‘palette’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:331)
GAPIT.Manhattan: no visible global function definition for ‘dev.off’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Manhattan.R:332)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘read.csv’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:45)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘write.csv’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:80)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘rainbow’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:164)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘pdf’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:187)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:188)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:194)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:199)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘read.csv’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:201)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘runif’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:271)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘mtext’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:291)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘abline’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:294)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘points’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:295)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘points’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:296)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘abline’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:302)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘abline’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:304)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘abline’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:309)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:310)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:311)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:312)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘mtext’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:313)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘dev.off’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:315)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘pdf’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:320)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:321)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:330)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:335)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘read.csv’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:339)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘runif’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:416)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘mtext’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:444)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘abline’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:447)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘abline’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:448)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘abline’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:451)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘points’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:453)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘points’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:454)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘abline’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:457)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:458)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:459)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘mtext’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:461)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘box’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:463)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘dev.off’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:465)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘pdf’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:490)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:491)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:492)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:500)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:501)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘abline’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:502)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘abline’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:507)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘abline’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:509)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘mtext’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:512)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘read.csv’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:522)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘runif’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:605)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:628)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘points’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:634)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘abline’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:638)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘points’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:640)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:642)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘dev.off’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:645)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘write.csv’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:692)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘pdf’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:693)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:694)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:695)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:702)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘points’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:714-716)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘text’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:717)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘points’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:726-728)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘text’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:729)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘points’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:739-741)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘text’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:742)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘points’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:748-750)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘text’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:751)
GAPIT.Multiple.Manhattan: no visible global function definition for
  ‘dev.off’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.Manhattan.R:757)
GAPIT.Multiple.QQ: no visible global function definition for ‘runif’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.QQ.R:19)
GAPIT.Multiple.QQ: no visible global function definition for ‘pdf’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.QQ.R:40)
GAPIT.Multiple.QQ: no visible global function definition for ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.QQ.R:41)
GAPIT.Multiple.QQ: no visible global function definition for ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.QQ.R:42)
GAPIT.Multiple.QQ: no visible global function definition for ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.QQ.R:43)
GAPIT.Multiple.QQ: no visible global function definition for ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple.QQ.R:96)
GAPIT.Multiple_Synthesis: no visible global function definition for
  ‘read.csv’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple_Synthesis.R:34)
GAPIT.Multiple_Synthesis: no visible global function definition for
  ‘read.csv’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple_Synthesis.R:122)
GAPIT.Multiple_Synthesis: no visible global function definition for
  ‘runif’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple_Synthesis.R:206)
GAPIT.Multiple_Synthesis: no visible global function definition for
  ‘add_markers’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple_Synthesis.R:244-270)
GAPIT.Multiple_Synthesis: no visible global function definition for
  ‘add_markers’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Multiple_Synthesis.R:277-303)
GAPIT.PCA: no visible global function definition for ‘legend’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.PCA.R:50-51)
GAPIT.PCA: no visible global function definition for ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.PCA.R:50-51)
GAPIT.PCA: no visible global function definition for ‘write.csv’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.PCA.R:62)
GAPIT.PCA: no visible global function definition for ‘legend’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.PCA.R:84-85)
GAPIT.PCA: no visible global function definition for ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.PCA.R:84-85)
GAPIT.PagainstP: no visible global function definition for ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.PagainstP.R:117)
GAPIT.PagainstP: no visible global function definition for ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.PagainstP.R:118)
GAPIT.PagainstP: no visible global function definition for ‘abline’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.PagainstP.R:126)
GAPIT.PagainstP: no visible global function definition for ‘abline’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.PagainstP.R:127)
GAPIT.PagainstP: no visible global function definition for ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.PagainstP.R:131)
GAPIT.PagainstP: no visible global function definition for ‘cor’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.PagainstP.R:139)
GAPIT.PagainstP: no visible global function definition for ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.PagainstP.R:143)
GAPIT.PagainstP: no visible global function definition for ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.PagainstP.R:144)
GAPIT.PagainstP: no visible global function definition for ‘mtext’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.PagainstP.R:145)
GAPIT.PagainstP: no visible global function definition for ‘mtext’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.PagainstP.R:146)
GAPIT.PagainstP: no visible global function definition for ‘legend’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.PagainstP.R:148-150)
GAPIT.Phenotype.View: no visible global function definition for
  ‘layout’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Phenotype.View.R:27-29)
GAPIT.Phenotype.View: no visible global function definition for ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Phenotype.View.R:32)
GAPIT.Phenotype.View: no visible global function definition for ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Phenotype.View.R:34)
GAPIT.Phenotype.View: no visible global function definition for ‘hist’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Phenotype.View.R:35)
GAPIT.Phenotype.View: no visible global function definition for ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Phenotype.View.R:36)
GAPIT.Phenotype.View: no visible global function definition for
  ‘density’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Phenotype.View.R:37)
GAPIT.Phenotype.View: no visible global function definition for
  ‘na.omit’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Phenotype.View.R:37)
GAPIT.Phenotype.View: no visible global function definition for ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Phenotype.View.R:38)
GAPIT.Phenotype.View: no visible global function definition for
  ‘boxplot’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Phenotype.View.R:39)
GAPIT.Phenotype.View: no visible global function definition for ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Phenotype.View.R:40)
GAPIT.Phenotype.View: no visible global function definition for ‘ecdf’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Phenotype.View.R:41)
GAPIT.Phenotype.afterGWAS: no visible global function definition for
  ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Phenotype.afterGWAS.R:53)
GAPIT.Phenotype.afterGWAS: no visible global function definition for
  ‘boxplot’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Phenotype.afterGWAS.R:85-87)
GAPIT.Phenotype.afterGWAS: no visible global function definition for
  ‘points’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Phenotype.afterGWAS.R:91)
GAPIT.Phenotype.afterGWAS: no visible global function definition for
  ‘runif’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Phenotype.afterGWAS.R:91)
GAPIT.Phenotype.afterGWAS: no visible global function definition for
  ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Phenotype.afterGWAS.R:94)
GAPIT.Phenotype.afterGWAS: no visible global function definition for
  ‘axis’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Phenotype.afterGWAS.R:95)
GAPIT.Phenotype.afterGWAS: no visible global function definition for
  ‘mtext’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Phenotype.afterGWAS.R:98)
GAPIT.Phenotype.afterGWAS: no visible global function definition for
  ‘mtext’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Phenotype.afterGWAS.R:99)
GAPIT.Power.compare: no visible binding for global variable ‘myGM’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Power.compare.R:50)
GAPIT.Power.compare: no visible global function definition for
  ‘rainbow’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Power.compare.R:81)
GAPIT.RandomModel: no visible global function definition for ‘var’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.RandomModel.R:36)
GAPIT.RandomModel: no visible global function definition for ‘layout’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.RandomModel.R:200-202)
GAPIT.RandomModel: no visible global function definition for ‘par’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.RandomModel.R:203)
GAPIT.Remove.outliers: no visible global function definition for
  ‘quantile’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Remove.outliers.R:18)
GAPIT.Remove.outliers: no visible global function definition for ‘IQR’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Remove.outliers.R:20)
GAPIT.SS: no visible global function definition for ‘quantile’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.SS.R:213)
GAPIT.Validation: no visible global function definition for ‘cor’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Validation.R:45)
GAPIT.Validation: no visible global function definition for ‘cor’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Validation.R:46)
GAPIT.Validation: no visible global function definition for ‘cor’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Validation.R:61)
GAPIT.Validation: no visible global function definition for ‘cor’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.Validation.R:80)

GAPIT.cross_validation.compare: no visible binding for global variable
  ‘myGD’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.cross_validation.compare.R:23)
GAPIT.cross_validation.compare: no visible binding for global variable
  ‘rel’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.cross_validation.compare.R:126)
GAPIT.cross_validation.compare: no visible binding for global variable
  ‘rel’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.cross_validation.compare.R:132)
GAPIT.cross_validation.compare: no visible binding for global variable
  ‘rel’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.cross_validation.compare.R:150)
GAPIT.cross_validation.compare: no visible binding for global variable
  ‘rel’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.cross_validation.compare.R:151)
GAPIT.cross_validation.compare: no visible binding for global variable
  ‘rel’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.cross_validation.compare.R:157)
GAPIT.cross_validation.compare: no visible binding for global variable
  ‘rel’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.cross_validation.compare.R:158)
GAPIT.get.LL: no visible binding for global variable ‘var’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.get.LL.R:23-28)
emmreml: no visible global function definition for ‘ginv’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.emmreml.R:14)
emmreml: no visible global function definition for ‘optimize’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.emmreml.R:32-33)
emmreml: no visible global function definition for ‘ginv’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.emmreml.R:39)
emmreml: no visible global function definition for ‘ginv’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.emmreml.R:47)
emmreml: no visible global function definition for ‘ginv’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.emmreml.R:70-71)
emmreml: no visible global function definition for ‘ginv’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.emmreml.R:83-84)
emmreml: no visible global function definition for ‘ginv’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.emmreml.R:94)
emmreml: no visible global function definition for ‘pchisq’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.emmreml.R:99)
emmreml: no visible binding for global variable ‘p.adjust.methods’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.emmreml.R:100)
emmreml : <anonymous>: no visible global function definition for
  ‘p.adjust’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.emmreml.R:101-102)
emmreml: no visible global function definition for ‘pchisq’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.emmreml.R:104-105)
emmreml : <anonymous>: no visible global function definition for
  ‘p.adjust’
  (/home/knausb/gits/GAPIT.Rcheck/00_pkg_src/GAPIT/R/GAPIT.emmreml.R:106-107)


Addressed testthat::test_check() 2025-10
* Moved 3d.PDC.py from R/ to inst/
* Added escape (\%) for percentage signs in GAPIT.R
* Added if( exists("eig.R"") ){ rm(eig.R) } to GAPIT.EMMAxP3D
* Added if(SNP.P3D == TRUE & exists("eig.L"") ){ rm(eig.L) } to GAPIT.EMMAxP3D
* Added Lambdas[ Lambdas <= 0 ] <- 1e-9 to GAPIT.emma.REMLE
* Changed Depends: R (<= 4.4.1) to Depends: R (>= 4.4.1)
* updated GAPIT.Rproj

* Added roxygen to GAPIT.
* Omitted `install.packages()` calls.
* Replaced `library()` and `require()` calls with `package::function()`.
* Added tests for models MLM, GLM, CMLM, MMLM, SUPER, FarmCPU, gBLUP, and cBLUP.
* Created manual page for `GAPIT3::GAPIT()`.
* Created manual page for `GAPIT3::GAPIT.Compression.Visualization()`, added parameter `file.output`.
* Created manual page for `GAPIT3::GAPIT.ID()`, file.output is respected for testing, returns an invisible NULL to facilitate testing.
* Omitted uncalled function `BlinkR.SUB()`.
* Documented example data sets.



