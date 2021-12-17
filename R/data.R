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


#' 
#' @name
#' HAPMAP format data
#' 
#' @title
#' HAPMAP format data
#' 
#' 
#' @details 
#' Genotypic data in HAPMAP format.
#' 
#' Variants are in rows.
#' The first 11 columns contain information about each variant.
#' GAPIT uses columns 1 (variant name; rs), 3 (chromosome; chrom), and 4 (position; pos).
#' Subsequent columns are samples.
#' The example dataset consists of 3093 variants scored for 280 samples.
#' 
#' 
#' @examples  
#' \dontrun{
#' myGenoFile <- system.file("extdata", "mdp_genotype_test.hmp.txt.gz", package = "GAPIT3")
#' myGenotypes  <- read.table(myGenoFile, header = FALSE)
#' myGenotypes[1:4, 1:11]
#' myGenotypes[1:4, c(1, 12:20)]
#' }
#'
#'
NULL



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


