



test_that("GAPIT function works, FarmCPU model", {
  myPhenoFile <- system.file("extdata", "mdp_traits.txt.gz",
                             package = "GAPIT3")
  myGenoFile <- system.file("extdata", "mdp_genotype_test.hmp.txt.gz",
                            package = "GAPIT3")
  
  myPhenotypes <- read.table(myPhenoFile, header = TRUE)
#  myPhenotypes <- myPhenotypes[, 1:2]
  myPhenotypes <- myPhenotypes[, c(1, 3)]
  myGenotypes  <- read.table(myGenoFile, header = FALSE)
  
  #  setwd(tempdir())
  #  getwd()
  
  #devtools::load_all()
  #debug(FarmCPU.LM)
  myGAPIT <- GAPIT( Y = myPhenotypes,
                    G = myGenotypes,
                    PCA.total = 3,
                    file.output = FALSE,
                    model = "FarmCPU"
  )
  
  # list.files()
  #  unlink(gfiles)
  # list.files()
  
  expect_true(inherits(myGAPIT, "list"))
  expect_true(length(myGAPIT) == 5)
  expect_true(all(names(myGAPIT) == c("GWAS",
                                      "h2", "PCA",
                                      "GD", "GM")))
  
})


test_that("FarmCPU function works", {
  myPhenoFile <- system.file("extdata", "mdp_traits.txt.gz",
                             package = "GAPIT3")
  myGenoFile <- system.file("extdata", "mdp_numeric.txt.gz",
                            package = "GAPIT3")
  
  myPhenotypes <- read.table(myPhenoFile, header = TRUE)
  myPhenotypes <- myPhenotypes[, 1:2]
  myGenotypes  <- read.table(myGenoFile, header = FALSE)
  
  #  setwd(tempdir())
  #  getwd()
  
  # myFarmCPU <- FarmCPU( Y = myPhenotypes,
  #                   GD = myGenotypes,
  #                   file.output = FALSE
  # )
  
  
})
  