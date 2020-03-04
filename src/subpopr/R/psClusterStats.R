
#' collect results from PS files
#' report number of clusters per species at specified cut-off value
#' @param psFilePath
#' @param psCutoff cutoff value for calling clusters (default = 0.8)
#' @return nClus number of clusters (int)
getNumberOfClusters <- function(psFilePath, psCutoff){
  if(file.size(psFilePath) <= 1 ){return(NA)}
  psVals <- read.table(psFilePath,header=T,row.names=1,sep='\t')
  nClus <- as.integer(max(rownames(psVals)[psVals$x > psCutoff]))
  return(nClus)
}


#' compare different PS cutoff values
#' @param distMethod distance method to use (default = "mann")
#' @param cutOffVals vector of PS cut off values to test, e.g. c(0.6,0.7,0.8,0.9)
#' @return data frame with columns: speciesID, number of clusters, PS cutoff
collectPSCutoffClusterResults <- function(resultsDir,
                                          #writeDir=NULL,
                                          distMeth = "mann",
                                          cutOffVals = c(0.6,0.7,0.8,0.9)){
  # if(is.null(writeDir)){
  #   writeDir <- resultsDir
  # }
  allPSfiles <- list.files(recursive = T,
                           path = paste0(resultsDir,"/"),
                           pattern = paste0('.*_',distMeth,'_PS_values.tab'),
                           full.names = T)

  psCutoffClusterResults <- data.frame()

  # for each species and each cluster .pos file
  for (d in allPSfiles) {
    # get e.g. 537011
    species <- strsplit(basename(d) ,'_')[[1]][1]
    nClus <- sapply(cutOffVals, getNumberOfClusters,psFilePath = d)
    data.frame(species=species,cutOffVals=cutOffVals,nClus=nClus)
    psCutoffClusterResults <- rbind.data.frame(psCutoffClusterResults,
                                     data.frame(species=species,
                                                cutOffVals=cutOffVals,
                                                nClus=nClus))

  }
  psCutoffClusterResults$distMethod = distMeth
  # write.table(psCutoffClusterResults,
  #            file=paste0(writeDir,"/psCutoffClusterResults_",distMeth,".tsv"),
  #            sep = "\t",row.names = F,quote = F)

  return(psCutoffClusterResults)
}
