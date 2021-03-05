

# works on all results files at once
pyGetPlacingRelevantSubset <- function(outDir,metaSnvDir,scriptDir){
  path <- paste0(scriptDir, "/getGenotypingSNVSubset.py")
  if(!file.exists(path)){
    stop("Missing python script, expected: ",path)
  }
  system(paste("python",path,outDir,metaSnvDir,sep=" ")) # write results to outDir
}


#'@param posFile output from getPlacingRelevantSubset.py e.g. 537011_2.pos
#'@param minDepth (int) if vertical coverage at this position is less than 'x'minDepth' in a sample, then set the SNV frequency to -1 which will be NA later
pyConvertSNPtoAllelTable <- function(posFile, minDepth = 5, scriptDir){
  if(file.size(posFile) <= 1){ stop(paste("File is empty:",posFile))}
  path <- paste0(scriptDir, "/convertSNVtoAlleleFreq.py")
  if(!file.exists(path)){
    stop("Missing python script, expected: ",path)
  }
  system(paste("python",path,posFile,minDepth,sep=" "))
  #print(paste("python",path,posFile,minDepth,sep=" "))
}
