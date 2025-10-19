
#' 
#' @name
#' Numeric format data
#' 
#' @title
#' Numeric format genotypic data
#' 
#' 
#' @details 
#' Genetic data in numeric format.
#' Samples in rows.
#' The first column is the sample names (taxa).
#' Subsequent columns are genotypes.
#' Homozygotes are scored as either 0 or 2.
#' Heterozygotes are scored as 1.
#' 
#' 
#' 
#' @examples  
#' \dontrun{
#' myGDFile <- system.file("extdata", "mdp_numeric.txt.gz", package = "GAPIT3")
#' myGD  <- read.table(myGDFile, header = FALSE)
#' myGD[1:4, 1:6]
#' }
#'
NULL


