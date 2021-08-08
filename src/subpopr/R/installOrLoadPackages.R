# this function will load all required R packages and, if doInstall=TRUE, it will install any that are missing
installOrLoadPackages <- function(doInstall=FALSE){
  biocManPackages <- c("BiocParallel")
  requiredPackageNames <- c("fpc","ape","getopt","optparse","readr","dplyr",
                            "ggplot2","tidyr","gridExtra","DT","cluster","ggrepel",
                            "data.table","rmarkdown","batchtools","futile.logger",
                            "Cairo","reshape2")
  pkgsLoaded <- sapply(requiredPackageNames,installOrLoadPackage,doInstall=doInstall)
}

# this function will load a package and, if doInstall=TRUE, it will install it if not yet installed
installOrLoadPackage <- function(pkgName, doInstall=FALSE){
  if(doInstall){
    if(pkgName %in% biocManPackages){
      if(!requireNamespace("BiocManager", quietly=TRUE)){install.packages("BiocManager", repos="https://cloud.r-project.org/")}
      if(!requireNamespace(pkgName, quietly=TRUE)){BiocManager::install(pkgName)}
    }else{
      if(!requireNamespace(pkgName, quietly=TRUE)){
        install.packages(pkgName, repos="https://cloud.r-project.org/")
      }
    }
  }
  library(pkgName,character.only = TRUE)
}
