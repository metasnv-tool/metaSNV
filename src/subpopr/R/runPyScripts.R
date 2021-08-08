
# checks for python3 in the system path and returns the path to the executable
getPython3Path <- function(){
  py3Path <- Sys.which("python3")
  if(py3Path==""){ #not in path
    py3Path <- Sys.which("python")
  }

  if(py3Path==""){ #python and python3 not in path
    stop("Python ('python' or 'python3') is not in the system path. Aborting.")
  }

  #check that python is version 3
  pyV <- system2(py3Path,"-V",stdout = T,stderr = T)
  # should return e.g. "Python 2.7.16" or "Python 3.8.2"
  pyV <- sub(pattern = "Python ",replacement = "",x = pyV)

  if(pyV < 3){
    stop("Python is not version 3 or greater. Using ",py3Path," which is version ",pyV,". ",
         "Please update python or make python3 available on your system path. Aborting.")
  }
  return(py3Path)
}


# works on all results files at once
pyGetPlacingRelevantSubset <- function(outDir,metaSnvDir,scriptDir){
  path <- paste0(scriptDir, "/getGenotypingSNVSubset.py")
  if(!file.exists(path)){
    stop("Missing python script, expected: ",path)
  }
  py3Path <- getPython3Path()
  system(paste(py3Path,path,outDir,metaSnvDir,sep=" ")) # write results to outDir
}


#'@param posFile output from getPlacingRelevantSubset.py e.g. 537011_2.pos
#'@param minDepth (int) if vertical coverage at this position is less than 'x'minDepth' in a sample, then set the SNV frequency to -1 which will be NA later
pyConvertSNPtoAllelTable <- function(posFile, minDepth = 5, scriptDir){
  if(file.size(posFile) <= 1){ stop(paste("File is empty:",posFile))}
  path <- paste0(scriptDir, "/convertSNVtoAlleleFreq.py")
  if(!file.exists(path)){
    stop("Missing python script, expected: ",path)
  }
  py3Path <- getPython3Path()
  system(paste(py3Path,path,posFile,minDepth,sep=" "))
}
