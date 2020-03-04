


prepareMetaSNVresults <- function(species, metaSNVDir){

  metaSNVDirDist <- paste0(metaSNVdir,"/distances/")
  metaSNVDirSNVsFiltered <- paste0(metaSNVdir,"/filtered/pop/")
  return(
    prepareMetaSNVresults_motus(species, metaSNVDirDist, metaSNVDirSNVsFiltered)
    )

}


prepareMetaSNVresults_motus <- function(species, metaSNVDirDist, metaSNVDirSNVsFiltered){
  # check for required files
  if(!dir.exists(metaSNVDirDist)){
    stop("Missing directory: ", metaSNVDirDist)
  }
  if(!dir.exists(metaSNVDirSNVsFiltered)){
    stop("Missing directory: ", metaSNVDirSNVsFiltered)
  }
  if( ! file.exists(distanceMatrixFileMann) ){
    stop("Missing file: ", distanceMatrixFileMann)
  }
  if( ! file.exists(distanceMatrixFileAllele) ){
    stop("Missing file: ", distanceMatrixFileAllele)
  }
  if( ! file.exists(freqCompFile) ){
    stop("Missing file: ", freqCompFile)
  }

  metaSNVresultPaths <- list()
  metaSNVresultPaths$distanceMann <- paste0(metaSNVDirDist,"/",species,".filtered.mann.dist")
  metaSNVresultPaths$distanceAllele <- paste0(metaSNVDirDist,"/",species,".filtered.allele.dist")
  metaSNVresultPaths$filteredSNVs <- paste0(metaSNVDirSNVsFiltered,"/",species,".filtered.freq")

  metaSNVresults <- list()
  metaSNVresults$distMann <- read.table(metaSNVresultPaths$distanceMann,header=T,row.names=1,check.names=F)
  metaSNVresults$distAllele <- read.table(metaSNVresultPaths$distanceAllele,header=T,row.names=1,check.names=F)
  metaSNVresults$snvFreqs.filtered <- read.table(metaSNVresultPaths$filteredSNVs,header=T,row.names=1,check.names=F)
}
