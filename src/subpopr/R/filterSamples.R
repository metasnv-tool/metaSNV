#' remove outlier samples, they can mess up clustering
#'@param distMatrix matrix of dissimilarities / distances or similarities
#'@param maxTimesSd mean dissimilarity for this observation can
#'be at most this many standard deviations away (+/-) from the mean of
#'the per-sample mean dissimilarities
#'@param maxNoutliers if there are more than this many "outliers" then maybe they aren't outliers
#'so don't do anything, just return the original matrix
#'@return list with two named entries: "distMatrix" = dist matrix with outliers removed
#'and "outliersRemoved" = vector of samples names removed
removeOutliersFromDistMatrixMeanDissim <- function(distMatrix, maxTimesSd=3, maxNoutliers=3,
                                         warn = T, species="") {
  dd <- distMatrix
  diag(dd) <- NA
  # based on mean dissimilarities
  meanDissimPerSample <- rowSums(dd,na.rm = T)/ncol(dd)
  upperThreshold <- mean(meanDissimPerSample)+maxTimesSd*sd(meanDissimPerSample)
  lowerThreshold <- mean(meanDissimPerSample)-maxTimesSd*sd(meanDissimPerSample)
  samplesToRm <- names(meanDissimPerSample[meanDissimPerSample>upperThreshold|
                                             meanDissimPerSample<lowerThreshold])

  if(length(samplesToRm) == 0 | length(samplesToRm) > maxNoutliers){
    return(distMatrix)
  }
  warning(paste("Removing",length(samplesToRm),"'outlier' samples from species:",species,
                "(",paste(samplesToRm,collapse = ","),")"))
  distMatrixClean <- distMatrix[!row.names(distMatrix) %in% samplesToRm,
                                        !colnames(distMatrix) %in% samplesToRm]
  return(list(distMatrix=distMatrixClean,
              outliersRemoved=samplesToRm))
}


#' remove outlier samples, they can mess up clustering
#'@param distMatrix matrix of dissimilarities / distances or similarities
#'@param maxTimesSd mean dissimilarity for this observation can
#'be at most this many standard deviations away (+/-) from the mean of
#'the per-sample mean dissimilarities
#'@param maxNoutliers if there are more than this many "outliers" then maybe they aren't outliers
#'so don't do anything, just return the original matrix
#'@return list with two named entries: "distMatrix" = dist matrix with outliers removed
#'and "outliersRemoved" = vector of samples names removed
removeOutliersFromDistMatrixMinDissim <- function(distMatrix, maxTimesSd=3, maxNoutliers=3,
                                         warn = T, species="") {
  dd <- distMatrix
  diag(dd) <- NA
  # based on min dissim -- e.g. if 3 samples only look like each other, they're ok
  # if a sample only looks like itself, then it's not useful
  minDissimPerSample <- apply(dd,1,min,na.rm=T)

  upperThreshold <- mean(minDissimPerSample)+maxTimesSd*sd(minDissimPerSample)
  lowerThreshold <- mean(minDissimPerSample)-maxTimesSd*sd(minDissimPerSample)
  samplesToRm <- names(minDissimPerSample[minDissimPerSample>upperThreshold|
                                               minDissimPerSample<lowerThreshold])

  if( length(samplesToRm) == 0 ){
    distMatrixClean <- distMatrix
    outliersRemoved <- c()
  } else if(length(samplesToRm) > maxNoutliers){
    warning(paste("Too many samples (",length(samplesToRm)," samples) look like 'outliers' for species:",species,
                  "(samples are: ",paste(samplesToRm,collapse = ","),"). Not removing any samples."))
    distMatrixClean <- distMatrix
    outliersRemoved <- c()
  }else{
      warning(paste("Removing",length(samplesToRm),"'outlier' samples from species:",species,
                "(",paste(samplesToRm,collapse = ","),")"))
    distMatrixClean <- distMatrix[!row.names(distMatrix) %in% samplesToRm,
                                !colnames(distMatrix) %in% samplesToRm]
    outliersRemoved <- samplesToRm
  }
  return(list(distMatrix=distMatrixClean,
              outliersRemoved=outliersRemoved))
}
