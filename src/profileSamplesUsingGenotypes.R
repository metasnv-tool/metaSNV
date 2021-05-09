# Run as:
# Rscript profileSamplesUsingGenotypes.R

library("optparse")

option_list = list(
  make_option(c("-e", "--metaSnvSrcDir"), type="character", default=NULL,
              help="Path to directory containing metaSNV executables. Required.", metavar="character"),
  make_option(c("-m", "--metaSNVResultsDir"), type="character", default=NULL,
              help="Path to directory containing results from metaSNV.py run (e.g. contains folder called snpCaller). Required.", metavar="character"),
  make_option(c("-s", "--subpoprResultsDir"), type="character", default=NULL,
              help="Path to directory containing results from metaSNV_subpopr.py run. Required.", metavar="character"),
  make_option(c("-i", "--speciesID"), type="character", default=NULL,
              help="ID of species to use (e.g. 155864). Required.", metavar="character"),
  make_option(c("-o", "--outDir"), type="character", default="results",
              help="directory to create for results output [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

checkValue <- function(x,xStr){
  if(is.null(x)){
    print_help(opt_parser)
    stop("Value for ",xStr," must be supplied", call.=FALSE)
  }
}

checkValue(opt$metaSnvSrcDir,"metaSnvSrcDir")
checkValue(opt$metaSNVResultsDir,"metaSNVResultsDir")
checkValue(opt$subpoprResultsDir,"subpoprResultsDir")
checkValue(opt$speciesID,"speciesID")

checkDirExists <- function(x){
  if(!dir.exists(x)){
    stop("Required directory does not exist: ",x)
  }
}
checkDirExists(opt$metaSnvSrcDir)
checkDirExists(opt$metaSNVResultsDir)
checkDirExists(opt$subpoprResultsDir)

metaSnvSrcDir = opt$metaSnvSrcDir
metaSNVResultsDir=opt$metaSNVResultsDir
subpoprResultsDir=opt$subpoprResultsDir
speciesID=opt$speciesID
outDir=opt$outDir

 # metaSnvSrcDir="/g/scb2/bork/rossum/metaSNV2/metaSNV"
 # metaSNVResultsDir="/g/scb2/bork/rossum/metagenomes/mock/inSilico/subpopr/posCtrl/v2_ecoli/realv1_x20/metaSNV/outputs"
 # subpoprResultsDir="/g/scb2/bork/rossum/metagenomes/human/subspecGeoValidation/all_v4/subpopr/results_mimicCostea_408121/params.hr10.hs80.ps80.gs80/outputs/"
 # outDir="/g/scb2/bork/rossum/metagenomes/human/subspecGeoValidation/all_v4/subpopr/results_mimicCostea_408121/params.hr10.hs80.ps80.gs80/placing/ecoli/composedSamples/results1/"
 # speciesID="155864"


print("Getting required files...")

posPaths = list.files(path =subpoprResultsDir,
                      pattern = paste0(speciesID,"_.*hap_po.*\\..*"),
                      full.names = T)
hapMedPaths = list.files(path =subpoprResultsDir,
                      pattern = paste0(speciesID,"_.*hap.*median.*\\..*"),
                      full.names = T)
if(!all(file.exists(posPaths))){
  stop("Required files do not exist: ",posPaths)
}
if(!all(file.exists(hapMedPaths))){
  stop("Required files do not exist: ",hapMedPaths)
}

dir.create(path = outDir,recursive = T)
if(!dir.exists(outDir)){
  stop("Error creaing output dir: ",outDir)
}

resCopy <- sapply(c(posPaths,hapMedPaths),file.copy,to=outDir)

newPosPaths = list.files(path =outDir,
                      pattern = paste0(speciesID,"_.*hap_po.*\\..*"),
                      full.names = T)
if(!all(file.exists(newPosPaths))){
  stop("Error copying files? Required files do not exist: ",posPaths)
}

minDepth=1

print("Getting genotyping positions for all samples...")

source("/g/scb2/bork/rossum/metaSNV2/metaSNV/src/subpopr/R/runPyScripts.R")
#python $metaSnvSrcDir/src/subpopr/inst/getGenotypingSNVSubset.py  outDir $metaSNVResultsDir
pyGetPlacingRelevantSubset(outDir = outDir,
                           metaSnvDir = metaSNVResultsDir,
                           scriptDir = paste0(metaSnvSrcDir,"/src/subpopr/inst/"))

posPaths2 = list.files(path =outDir,
                      pattern = paste0(speciesID,"_[[:digit:]]+.pos$"),
                      full.names = T)

print("Compiling genotyping allele frequencies for all samples...")
tmp <- sapply(X = posPaths2,
       FUN = pyConvertSNPtoAllelTable,
       minDepth = minDepth,
       scriptDir = paste0(metaSnvSrcDir,"/src/subpopr/inst/"))

source(paste0(metaSnvSrcDir,"/src/subpopr/R/installOrLoadPackages.R"))
source(paste0(metaSnvSrcDir,"/src/subpopr/R/utils.R"))
source(paste0(metaSnvSrcDir,"/src/subpopr/R/writeSubpopsForAllSamples.R"))
source(paste0(metaSnvSrcDir,"/src/subpopr/R/profileSubpops.R"))

print("Compiling subspecies profiles for all samples...")
tmp <- useGenotypesToProfileSubpops(species=speciesID,
                             metaSNVdir=paste0(metaSNVResultsDir,"/"),
                             outDir=paste0(outDir,"/"))

print("Complete.")
