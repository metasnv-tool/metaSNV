# remove rows and columns that are entirely NAs
rmNAfromDistMatrix <- function(distMa) {
  distMa <- distMa[rowSums(is.na(distMa))<nrow(distMa),colSums(is.na(distMa))<nrow(distMa)]
  # remove fewest rows and columns that will still minimise NAs
  while(any(is.na(distMa)) && nrow(distMa) > 0 ){
    sampleWithMaxNA <- names(sort(colSums(is.na(distMa)),decreasing = T)[1])
    distMa <- distMa[names(distMa) != sampleWithMaxNA, names(distMa) != sampleWithMaxNA]
  }
  return(distMa)
}