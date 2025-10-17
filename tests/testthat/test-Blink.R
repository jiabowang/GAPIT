


test_that("Blink function works", {
  myPhenoFile <- system.file("extdata", "mdp_traits.txt.gz",
                             package = "GAPIT")
  myGenomeDataFile <- system.file("extdata", 
                                 "mdp_numeric.txt.gz",
                                 package = "GAPIT")  
  myGenomeMapFile <- system.file("extdata", 
                                 "mdp_SNP_information.txt.gz",
                                 package = "GAPIT")
  
  myPhenotypes <- read.table(myPhenoFile, header = TRUE)
  myPhenotypes <- myPhenotypes[, c(1, 3)]
  myGD  <- read.table(myGenomeDataFile, header = TRUE)
  myGM  <- read.table(myGenomeMapFile, header = TRUE)
  
  # SNP name is key to GD and GM
  # all(colnames(myGD)[-1] == myGM$SNP)

  mySamps <- intersect(myGD$taxa, myPhenotypes$Taxa)
  rownames(myPhenotypes) <- myPhenotypes$Taxa
  myPhenotypes <- myPhenotypes[mySamps,]
  rownames(myGD) <- myGD$taxa
  myGD <- myGD[mySamps, ]
  # all(myGD$taxa == myPhenotypes$Taxa)

#  myBlink <- Blink(Y = myPhenotypes, GD = myGD, GM = myGM)

  # myBlink <- GAPIT( Y = myPhenotypes,
  #                   G = myGD,
  #                   PCA.total = 3,
  #                   file.output = FALSE,
  #                   model = "Blink"
  # )
  # 
  # expect_true(inherits(myBlink, "list"))
  # expect_true(length(myBlink) == 4)
  # expect_true(all(names(myBlink) == c("GWAS", "myGLM", "PEV", "seqQTN")))

})


# test_that("Blink.BICselection works", {
#   myPhenoFile <- system.file("extdata", "mdp_traits.txt.gz",
#                              package = "GAPIT3")
#   myPhenotypes <- read.table(myPhenoFile, header = TRUE)
#   #  myPhenotypes <- myPhenotypes[, 1:2]
#   myPhenotypes <- myPhenotypes[, c(1, 3)]
# 
#   myGenomeDataFile <- system.file("extdata", 
#                                   "mdp_numeric.txt.gz",
#                                   package = "GAPIT3")
#   myGD  <- read.table(myGenomeDataFile, header = TRUE)
#     
#   # devtools::load_all()
#   # debug(Blink.BICselection)
#   my_Blink <- Blink.BICselection(Y = myPhenotypes, GD = myGD)
# })
