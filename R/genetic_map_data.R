
#' 
#' @name
#' Genetic Map data
#' 
#' @title
#' Genetic Map data
#' 
#' @details 
#' Genetic Map to provide genomic coordinates for GD.
#' Variants are in rows.
#' The first column is the variant name (SNP).
#' The second column is the chromosome (Chromosome).
#' The third column is the position (position).
#' 
#' @seealso Numeric format data
#' 
#' 
#' @examples  
#' \dontrun{
#' myGMFile <- system.file("extdata", "mdp_SNP_information.txt.gz", package = "GAPIT3")
#' myGM  <- read.table(myGMFile, header = FALSE)
#' myGM[1:4, ]
#' }
#'
NULL


