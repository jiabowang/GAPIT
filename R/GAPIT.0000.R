.gapit_set_repos <- function(){
  r <- getOption("repos")
  if(is.null(r) || length(r)==0 || is.na(r["CRAN"]) || r["CRAN"]=="" || r["CRAN"]=="@CRAN@"){
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  }
}
.gapit_set_download_options <- function(){
  if(is.null(getOption("timeout")) || is.na(getOption("timeout"))){
    options(timeout = 600)
  }else{
    options(timeout = max(600, getOption("timeout")))
  }
}
.gapit_install_cran <- function(pkg){
  .gapit_set_repos()
  .gapit_set_download_options()
  try(install.packages(pkg, dependencies = TRUE), silent = TRUE)
}
.gapit_bioc_mirrors <- function(){
  c(
    "https://bioconductor.org",
    "https://bioconductor.posit.co",
    "https://mirrors.tuna.tsinghua.edu.cn/bioconductor",
    "https://mirrors.nju.edu.cn/bioconductor",
    "https://mirrors.ustc.edu.cn/bioc"
  )
}
.gapit_bioc_version_for_r <- function(){
  rv <- getRversion()
  if (rv >= "4.4.0") return("3.19")
  if (rv >= "4.3.0") return("3.18")
  if (rv >= "4.2.0") return("3.16")
  if (rv >= "4.1.0") return("3.14")
  return(BiocManager::version())
}
.gapit_install_bioc <- function(pkgs){
  pkgs <- as.character(pkgs)
  pkgs <- pkgs[!is.na(pkgs) & nzchar(pkgs)]
  pkgs <- unique(pkgs)
  if(!requireNamespace("BiocManager", quietly = TRUE)){
    .gapit_set_repos()
    .gapit_set_download_options()
    try(install.packages("BiocManager"), silent = TRUE)
  }
  if(requireNamespace("BiocManager", quietly = TRUE)){
    .gapit_set_download_options()
    ver <- .gapit_bioc_version_for_r()
    try(BiocManager::install(version = ver, ask = FALSE), silent = TRUE)
    mirrors <- .gapit_bioc_mirrors()
    for(m in mirrors){
      options(BioC_mirror = m)
      try(BiocManager::install(pkgs, ask = FALSE, update = FALSE), silent = TRUE)
      still_missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
      if(length(still_missing) == 0) return(invisible(TRUE))
      if(identical(.Platform$OS.type, "windows")){
        old_method <- getOption("download.file.method")
        on.exit(options(download.file.method = old_method), add = TRUE)
        options(download.file.method = "wininet")
        try(BiocManager::install(still_missing, ask = FALSE, update = FALSE), silent = TRUE)
        still_missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
        if(length(still_missing) == 0) return(invisible(TRUE))
      }
    }
  }
  invisible(all(sapply(pkgs, requireNamespace, quietly = TRUE)))
}
.gapit_require_or_install <- function(pkg, bioc = FALSE){
  if(!requireNamespace(pkg, quietly = TRUE)){
    cat("R:", R.version.string, "| Platform:", R.version$platform, "\n")
    if(isFALSE(bioc)) .gapit_install_cran(pkg) else .gapit_install_bioc(pkg)
  }
  if(!requireNamespace(pkg, quietly = TRUE)){
    stop(sprintf("Package '%s' not available after install attempt.", pkg), call. = FALSE)
  }
}
.gapit_require_or_install("gplots")
.gapit_require_or_install("genetics")
.gapit_require_or_install("ape")
.gapit_require_or_install("compiler")
.gapit_require_or_install("grid")
.gapit_require_or_install("bigmemory")
.gapit_require_or_install("EMMREML")
.gapit_require_or_install("scatterplot3d")
.gapit_require_or_install("lme4")
.gapit_require_or_install("multtest", bioc = TRUE)
.gapit_require_or_install("snpStats", bioc = TRUE)

`GAPIT.0000` <-
function(){
##############################################################################################
#GAPIT: Genome Association and Prediction Integrated Tool
#Objective 1: State of art methods for high  power, accuracy and speed;
#Objective 2: User friendly by design, help documents, and web forum;
#Objective 3: Comprehensive output to interpret data and results;
#Objective 4: Informative tables and high quality figures for reports and publication;

#Methods implimented: 
# 1. GLM (Structure or Q method for GWAS, Pritchard et. al. Genetics, 2000)
# 2. MLM (Q+K, Yu et. al. Nature Genetics, 2006)
# 3. gBLUP (Marker based kinship, Zhang et. al. Journal of Animal Science, 2007)
# 4. PCA (Zhao et. al. Plos Genetics, 2007)
# 5. EMMA (Kang et. al. Genetics, 2008)
# 6. CMLM (Zhang et. al. Nature Genetics, 2010)
# 7. EMMAx (Kang et. al. Nature Genetics, 2010)
# 8. P3D (Zhang et. al. Nature Genetics, 2010)
# 9. FaST-LMM (Lippert et. al. Nature Methods, 2011)
# 10. ECMLM (Li et. al. BMC Bioogy, 2014)
# 11. SUPER (Wang et. al. PLoS One, 2014)

GAPIT.Version="2026.3.10, GAPIT 4.1"
return(GAPIT.Version)
}
#=============================================================================================
