


test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


test_that("GAPIT mdp (Y and X) import works", {
  myPhenoFile <- system.file("extdata", "mdp_traits.txt.gz",
                             package = "GAPIT3documentation")
  myGenoFile <- system.file("extdata", "mdp_genotype_test.hmp.txt.gz",
                            package = "GAPIT3documentation")

  myPhenotypes <- read.table(myPhenoFile, header = TRUE)
  myGenotypes  <- read.table(myGenoFile, header = FALSE)

  expect_true( inherits(myPhenotypes, what = "data.frame") )
  expect_true( ncol(myPhenotypes) >= 2 )
  
  expect_true( inherits(myGenotypes, what = "data.frame") )
  expect_true( ncol(myGenotypes) >= 12 )
})


