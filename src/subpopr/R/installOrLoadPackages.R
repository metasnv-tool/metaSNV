requiredPackageNames_biocManPackages <- c("BiocParallel")
requiredPackageNames <- c("fpc","ape","getopt","optparse","readr","dplyr","forcats","kableExtra",
                          "ggplot2","tidyr","gridExtra","DT","cluster","ggrepel",
                          "data.table","rmarkdown",# for report rendering
                          "batchtools","futile.logger",
                          "Cairo","reshape2","BiocParallel")

# this function will just check whether packages are installed and return a vector of packages names that are missing
# if all packages are installed, will return NULL
getMissingPackages <- function(){
    #logical named vector with true if package if installed, false otherwise; names are package names
    pkgsInstalled <- sapply(requiredPackageNames,requireNamespace,quietly=TRUE)
    if(any(!pkgsInstalled)){
      return(names(pkgsInstalled[!pkgsInstalled]))
    }else{
      return(NULL)
    }
}


# this function will load all required R packages and, if doInstall=TRUE, it will install any that are missing
installOrLoadPackages <- function(doInstall=FALSE, doSuppressPackageStartupMessages=FALSE){

  # check that required variables exist
  if(!exists("requiredPackageNames") ||
    is.null(requiredPackageNames) ||
    length(requiredPackageNames) == 0){
    stop("Expected variable does not exist. Variable called requiredPackageNames, which ",
         "should be a vector of package names required, does not exist. Aborting.")
  }

  if(!exists("requiredPackageNames_biocManPackages") ||
     is.null(requiredPackageNames_biocManPackages) ||
     length(requiredPackageNames_biocManPackages) == 0){
    stop("Expected variable does not exist. Variable called requiredPackageNames_biocManPackages, which ",
         "should be a vector of package names required, does not exist. Aborting.")
  }

  pkgsLoaded <- sapply(requiredPackageNames,installOrLoadPackage,
                       doInstall=doInstall, doSuppressPackageStartupMessages=doSuppressPackageStartupMessages)
}

# this function will load a package and, if doInstall=TRUE, it will install it if not yet installed
installOrLoadPackage <- function(pkgName, doInstall=FALSE, doSuppressPackageStartupMessages=FALSE){
  if(doInstall){
    if(pkgName %in% requiredPackageNames_biocManPackages){
      if(!requireNamespace("BiocManager", quietly=TRUE)){install.packages("BiocManager", repos="https://cloud.r-project.org/")}
      if(!requireNamespace(pkgName, quietly=TRUE)){BiocManager::install(pkgName)}
    }else{
      if(!requireNamespace(pkgName, quietly=TRUE)){
        install.packages(pkgName, repos="https://cloud.r-project.org/")
      }
    }
  }
  if(doSuppressPackageStartupMessages){
    suppressPackageStartupMessages(library(pkgName,character.only = TRUE))
  }else{
    library(pkgName,character.only = TRUE)
  }

}
