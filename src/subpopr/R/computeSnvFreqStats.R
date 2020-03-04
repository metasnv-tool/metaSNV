computeSnvFreqStats <- function(snvFreqs) {
  #Compute allele frequency distributional properties
  # for each sample, what is the proportion of SNVs that have an extreme frequency e.g. <10 or >90
  # Samples that didn't have vertical coverage at the SNV position (-1) are ignored in 10/90 calc.
  freq_data_sample_20_80 <- apply(snvFreqs,2,function(x) {
    x <- x[x!=-1]
    x <- x[!is.na(x)]
    sum(x < 20 | x > 80)/length(x)
  })
  freq_data_sample_10_90 <- apply(snvFreqs,2,function(x) {
    x <- x[x!=-1]
    x <- x[!is.na(x)]
    sum(x < 10 | x > 90)/length(x)
  })
  freq_data_sample_5_95 <- apply(snvFreqs,2,function(x) {
    x <- x[x!=-1]
    x <- x[!is.na(x)]
    sum(x < 5 | x > 95)/length(x)
  })

  freq_data_sample <- data.frame(freq_data_sample_20_80,freq_data_sample_10_90,freq_data_sample_5_95)

  return(freq_data_sample)
}

#' Compute allele frequency distribution for each sample, what is the proportion of
#' SNVs that have an extreme frequency e.g. for homogThreshold=10: < 10% or > 90%
#' Samples that didn't have vertical coverage at the SNV position (-1) are ignored in calc.
#'@param homogThreshold the proportion away from 1 that is still considered homogeneous enough. Range 0-1. Default = 0.1
#'@return the proportion of SNVs that are "nearly homogeneous" per sample
computeSnvFreqStatsThreshold <- function(snvFreqs, homogThreshold= 0.1) {
  homogThreshold = homogThreshold * 100
  highThresh <- max(100-homogThreshold,homogThreshold)
  lowThresh <- min(100-homogThreshold,homogThreshold)

  #Compute allele frequency distributional properties
  # for each sample, what is the proportion of SNVs that have an extreme frequency e.g. <10 or >90
  # Samples that didn't have vertical coverage at the SNV position (-1) are ignored in 10/90 calc.
  # if sample has only -1 values, it will have NaN as the freq
  freq_data_sample <- apply(snvFreqs,2,function(x) {
    x <- x[x!=-1] # makes an empty vector if only -1 values
    x <- x[!is.na(x)]
    sum(x < lowThresh | x > highThresh)/length(x)
  })
  return(freq_data_sample)
}

majorAllel <- function(d) {
  ma <- apply(d,1,function(x) {
    x <- x[x!=-1]
    x <- x[!is.na(x)]
    x[x <50] <- 0
    x[x >=50] <- 1
    return(median(x,na.rm = T))
  })
  return(ma)
}

