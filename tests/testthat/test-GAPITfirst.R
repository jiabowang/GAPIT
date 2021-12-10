
# testthat::test_dir("tests/")
# library("testthat")

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


test_that("GAPIT mdp (Y and X) file import works", {
  myPhenoFile <- system.file("extdata", "mdp_traits.txt.gz",
                             package = "GAPIT3")
  myGenoFile <- system.file("extdata", "mdp_genotype_test.hmp.txt.gz",
                            package = "GAPIT3")

  myPhenotypes <- read.table(myPhenoFile, header = TRUE)
  myPhenotypes <- myPhenotypes[, 1:3]
  myGenotypes  <- read.table(myGenoFile, header = FALSE)


  expect_true( inherits(myPhenotypes, what = "data.frame") )
  expect_true( ncol(myPhenotypes) >= 2 )
  
  expect_true( inherits(myGenotypes, what = "data.frame") )
  expect_true( ncol(myGenotypes) >= 12 )
})


# gfiles <- c("GAPIT.Heterozygosity.pdf",
#             "GAPIT.Kin.VanRaden.csv",
#             "GAPIT.Kin.VanRaden.pdf",
#             "GAPIT.Marker.Density.pdf", "GAPIT.Marker.LD.pdf",
#             "GAPIT.MLM.dpoll.Df.tValue.StdErr.csv",
#             "GAPIT.MLM.dpoll.GWAS.Results.csv", "GAPIT.MLM.dpoll.Log.csv",
#             "GAPIT.MLM.dpoll.MAF.pdf",
#             "GAPIT.MLM.dpoll.Manhattan.Plot.Chromosomewise.pdf",
#             "GAPIT.MLM.dpoll.Manhattan.Plot.Genomewise.pdf",
#             "GAPIT.MLM.dpoll.Optimum.pdf",
#             "GAPIT.MLM.dpoll.phenotype_view.pdf",
#             "GAPIT.MLM.dpoll.PRED.csv",
#             "GAPIT.MLM.dpoll.QQ-Plot.pdf",
#             "GAPIT.MLM.dpoll.ROC.csv",
#             "GAPIT.MLM.dpoll.ROC.pdf",
#             "GAPIT.MLM.EarDia.Df.tValue.StdErr.csv",
#             "GAPIT.MLM.EarDia.GWAS.Results.csv",
#             "GAPIT.MLM.EarDia.Log.csv", "GAPIT.MLM.EarDia.MAF.pdf",
#             "GAPIT.MLM.EarDia.Manhattan.Plot.Chromosomewise.pdf",
#             "GAPIT.MLM.EarDia.Manhattan.Plot.Genomewise.pdf",
#             "GAPIT.MLM.EarDia.Optimum.pdf",
#             "GAPIT.MLM.EarDia.phenotype_view.pdf", "GAPIT.MLM.EarDia.PRED.csv",
#             "GAPIT.MLM.EarDia.QQ-Plot.pdf", "GAPIT.MLM.EarDia.ROC.csv",
#             "GAPIT.MLM.EarDia.ROC.pdf",
#             "GAPIT.MLM.EarHT.Df.tValue.StdErr.csv",
#             "GAPIT.MLM.EarHT.GWAS.Results.csv", "GAPIT.MLM.EarHT.Log.csv",
#             "GAPIT.MLM.EarHT.MAF.pdf",
#             "GAPIT.MLM.EarHT.Manhattan.Plot.Chromosomewise.pdf",
#             "GAPIT.MLM.EarHT.Manhattan.Plot.Genomewise.pdf",
#             "GAPIT.MLM.EarHT.Optimum.pdf",
#             "GAPIT.MLM.EarHT.phenotype_view.pdf",
#             "GAPIT.MLM.EarHT.PRED.csv", "GAPIT.MLM.EarHT.QQ-Plot.pdf",
#             "GAPIT.MLM.EarHT.ROC.csv", "GAPIT.MLM.EarHT.ROC.pdf",
#             "GAPIT.PCA.2D.pdf", "GAPIT.PCA.3D.pdf", "GAPIT.PCA.csv",
#             "GAPIT.PCA.eigenValue.pdf",
#             "GAPIT.PCA.eigenvalues.csv", "GAPIT.PCA.loadings.csv",
#             "libloc_166_37b72f10c7c2119.rds",
#             "libloc_213_36e896e939ef6a36.rds")


test_that("GAPIT function works, MLM model", {
  myPhenoFile <- system.file("extdata", "mdp_traits.txt.gz",
                             package = "GAPIT3")
  myGenoFile <- system.file("extdata", "mdp_genotype_test.hmp.txt.gz",
                            package = "GAPIT3")

  myPhenotypes <- read.table(myPhenoFile, header = TRUE)
  myPhenotypes <- myPhenotypes[, 1:3]
  myGenotypes  <- read.table(myGenoFile, header = FALSE)

#  setwd(tempdir())
#  getwd()

  myGAPIT <- GAPIT( Y = myPhenotypes,
                    G = myGenotypes,
                    PCA.total = 3,
                    file.output = FALSE,
                    model = "MLM"
                  )

# list.files()
#  unlink(gfiles)
# list.files()

  expect_true(inherits(myGAPIT, "list"))
  expect_true(length(myGAPIT) == 11)
})


test_that("GAPIT function works, GLM model", {
  myPhenoFile <- system.file("extdata", "mdp_traits.txt.gz",
                             package = "GAPIT3")
  myGenoFile <- system.file("extdata", "mdp_genotype_test.hmp.txt.gz",
                            package = "GAPIT3")
  
  myPhenotypes <- read.table(myPhenoFile, header = TRUE)
  myPhenotypes <- myPhenotypes[, 1:3]
  myGenotypes  <- read.table(myGenoFile, header = FALSE)
  
#  setwd(tempdir())
#  getwd()
  
  myGAPIT <- GAPIT( Y = myPhenotypes,
                    G = myGenotypes,
                    PCA.total = 3,
                    file.output = FALSE,
                    model = "GLM"
  )
  
  # list.files()
#  unlink(gfiles)
  # list.files()
  
  expect_true(inherits(myGAPIT, "list"))
  expect_true(length(myGAPIT) == 11)
})



test_that("GAPIT function works, CMLM model", {
  myPhenoFile <- system.file("extdata", "mdp_traits.txt.gz",
                             package = "GAPIT3")
  myGenoFile <- system.file("extdata", "mdp_genotype_test.hmp.txt.gz",
                            package = "GAPIT3")
  
  myPhenotypes <- read.table(myPhenoFile, header = TRUE)
  myPhenotypes <- myPhenotypes[, 1:3]
  myGenotypes  <- read.table(myGenoFile, header = FALSE)
  
#  setwd(tempdir())
#  getwd()
  
  myGAPIT <- GAPIT( Y = myPhenotypes,
                    G = myGenotypes,
                    PCA.total = 3,
                    file.output = FALSE,
                    model = "CMLM"
  )
  
  # list.files()
#  unlink(gfiles)
  # list.files()
  
  expect_true(inherits(myGAPIT, "list"))
  expect_true(length(myGAPIT) == 11)
  expect_true(all(names(myGAPIT) == c("GWAS", "Pred", "mc", "bc", "mp",
                                      "h2", "PCA", "GD", "GM", 
                                      "KI", "Compression")))
})


test_that("GAPIT function works, MMLM model", {
  myPhenoFile <- system.file("extdata", "mdp_traits.txt.gz",
                             package = "GAPIT3")
  myGenoFile <- system.file("extdata", "mdp_genotype_test.hmp.txt.gz",
                            package = "GAPIT3")
  
  myPhenotypes <- read.table(myPhenoFile, header = TRUE)
  myPhenotypes <- myPhenotypes[, 1:3]
  myGenotypes  <- read.table(myGenoFile, header = FALSE)
  
#  setwd(tempdir())
#  getwd()
  
  myGAPIT <- GAPIT( Y = myPhenotypes,
                    G = myGenotypes,
                    PCA.total = 3,
                    file.output = FALSE,
                    model = "MMLM"
  )
  
  # list.files()
#  unlink(gfiles)
  # list.files()
  
  expect_true(inherits(myGAPIT, "list"))
  expect_true(length(myGAPIT) == 11)
  expect_true(all(names(myGAPIT) == c("GWAS", "Pred", "mc", "bc", "mp",
                                      "h2", "PCA", "GD", "GM", 
                                      "KI", "Compression")))
})


test_that("GAPIT function works, SUPER model", {
  myPhenoFile <- system.file("extdata", "mdp_traits.txt.gz",
                             package = "GAPIT3")
  myGenoFile <- system.file("extdata", "mdp_genotype_test.hmp.txt.gz",
                            package = "GAPIT3")
  
  myPhenotypes <- read.table(myPhenoFile, header = TRUE)
  myPhenotypes <- myPhenotypes[, 1:3]
  myGenotypes  <- read.table(myGenoFile, header = FALSE)
  
#  setwd(tempdir())
#  getwd()
  
  myGAPIT <- GAPIT( Y = myPhenotypes,
                    G = myGenotypes,
                    PCA.total = 3,
                    file.output = FALSE,
                    model = "SUPER"
  )
  
  # list.files()
#  unlink(gfiles)
  # list.files()
  
  expect_true(inherits(myGAPIT, "list"))
  expect_true(length(myGAPIT) == 9)
  expect_true(all(names(myGAPIT) == c("GWAS", "Pred", "mc", "mp",
                                      "PCA", "GD", "GM", 
                                      "KI", "Compression")))
})


# test_that("GAPIT function works, FarmCPU model", {
#   myPhenoFile <- system.file("extdata", "mdp_traits.txt.gz",
#                              package = "GAPIT3")
#   myGenoFile <- system.file("extdata", "mdp_genotype_test.hmp.txt.gz",
#                             package = "GAPIT3")
#   
#   myPhenotypes <- read.table(myPhenoFile, header = TRUE)
#   myPhenotypes <- myPhenotypes[, 1:3]
#   myGenotypes  <- read.table(myGenoFile, header = FALSE)
#   
# #  setwd(tempdir())
# #  getwd()
#   
#   myGAPIT <- GAPIT( Y = myPhenotypes,
#                     G = myGenotypes,
#                     PCA.total = 3,
#                     file.output = FALSE,
#                     model = "FarmCPU"
#   )
#   
#   # list.files()
# #  unlink(gfiles)
#   # list.files()
#   
#   expect_true(inherits(myGAPIT, "list"))
#   expect_true(length(myGAPIT) == 5)
#   expect_true(all(names(myGAPIT) == c("GWAS",
#                                       "h2", "PCA",
#                                       "GD", "GM")))
#   
# })


test_that("GAPIT function works, gBLUP model", {
  myPhenoFile <- system.file("extdata", "mdp_traits.txt.gz",
                             package = "GAPIT3")
  myGenoFile <- system.file("extdata", "mdp_genotype_test.hmp.txt.gz",
                            package = "GAPIT3")
  
  myPhenotypes <- read.table(myPhenoFile, header = TRUE)
  myPhenotypes <- myPhenotypes[, 1:3]
  myGenotypes  <- read.table(myGenoFile, header = FALSE)
  
#  setwd(tempdir())
#  getwd()
  
  myGAPIT <- GAPIT( Y = myPhenotypes,
                    G = myGenotypes,
                    PCA.total = 3,
                    file.output = FALSE,
                    model = "gBLUP"
  )
  
  # list.files()
#  unlink(gfiles)
  # list.files()
  
  expect_true(inherits(myGAPIT, "list"))
  expect_true(length(myGAPIT) == 7)
  expect_true(all(names(myGAPIT) == c("Pred",
                                      "h2", "PCA", "GD", "GM", 
                                      "KI", "Compression")))
})

test_that("GAPIT function works, cBLUP model", {
  myPhenoFile <- system.file("extdata", "mdp_traits.txt.gz",
                             package = "GAPIT3")
  myGenoFile <- system.file("extdata", "mdp_genotype_test.hmp.txt.gz",
                            package = "GAPIT3")
  
  myPhenotypes <- read.table(myPhenoFile, header = TRUE)
  myPhenotypes <- myPhenotypes[, 1:3]
  myGenotypes  <- read.table(myGenoFile, header = FALSE)
  
#  setwd(tempdir())
#  getwd()
  
  myGAPIT <- GAPIT( Y = myPhenotypes,
                    G = myGenotypes,
                    PCA.total = 3,
                    file.output = FALSE,
                    model = "cBLUP"
  )
  
  # list.files()
#  unlink(gfiles)
  # list.files()
  
  expect_true(inherits(myGAPIT, "list"))
  expect_true(length(myGAPIT) == 7)
  expect_true(all(names(myGAPIT) == c("Pred",
                                      "h2", "PCA", "GD", "GM", 
                                      "KI", "Compression")))
})


