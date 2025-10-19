#' 
#' @name 
#' Phenotypic data
#'
#' @title 
#' Phenotypic data
#' 
#' 
#' @details 
#' Samples are in rows.
#' The first row is a header.
#' The first column contains the sample names (Taxa).
#' Subsequent columns are phenotypic traits.
#' The example dataset consists of 301 samples scored for three traits.
#' 
#' 
#' @examples  
#' \dontrun{
#' myPhenoFile <- system.file("extdata", "mdp_traits.txt.gz", package = "GAPIT3")
#' myPhenotypes <- read.table(myPhenoFile, header = TRUE)
#' head(myPhenotypes)
#' }
#'
NULL

