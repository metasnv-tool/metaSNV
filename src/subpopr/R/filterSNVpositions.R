
rmLowVarSNVs <- function(snvFreqs){
  #get the mean frequency of each SNV
  freq_means <- apply(snvFreqs,1,function(x) {
    # don't include -1 in calc of means
    a <- which(x==-1)
    if (length(a) > 0) {
      return(mean(x[-a]))
    } else {
      return(mean(x))
    }
  })


  # remove SNV positions with low variability among samples
  # i.e. those that have a mean frequency of >95% or <5%
  r <- c(which( freq_means > 95),which( freq_means < 5))
  if (length(r) > 0) {
    snvFreqs <- snvFreqs[-r,]
  }

  return(snvFreqs)
}

