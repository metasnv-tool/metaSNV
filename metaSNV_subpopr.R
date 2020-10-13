#!/usr/bin/env Rscript

# qlogin -pe smp 12 -l h_vmem=10G # 10G memory per core is plenty
# REQUIRES PYTHON 3
# Rscript metaSNV_supopr.R -h

# clear the environment

rm(list=ls())

normalRun<-TRUE # use cmd line args
useExistingClustering<-TRUE
useExistingExtension<-TRUE
makeReports <- FALSE
makeGeneReports <- TRUE
doSpecies<-FALSE

ptm <- proc.time()
suppressPackageStartupMessages(library(futile.logger))
tmp <- flog.threshold(INFO) # assign to tmp to avoid NULL being returned and printed

if(normalRun){
  # Expectation is that this script will sit in the metaSNV directory,
  # which will include a directory ./src/subpopr
  
  # try to set the current working directory to the location of this file
  # works if this file is sourced() or has been called from cmd line (e.g. Rscript metaSNV_subpop.R [...])
  thisFile <- function() {
    cmdArgs <- commandArgs(trailingOnly = FALSE)
    needle <- "--file="
    match <- grep(needle, cmdArgs)
    if (length(match) > 0) {
      # Rscript
      return(normalizePath(sub(needle, "", cmdArgs[match])))
    } else {
      # 'source'd via R console
      return(normalizePath(sys.frames()[[1]]$ofile))
    }
  }
  
  scriptDir <- dirname(thisFile())
  
  
  # Parse params -------------------------------------------------------------
  
  suppressPackageStartupMessages(library(getopt))
  suppressPackageStartupMessages(library(optparse))
  
  option_list = list(
    make_option(c("-i", "--metaSnvResultsDir"), type="character",
                default=NULL,
                help="Path to directory that has the metaSNV results, used as input (required)",
                metavar="file path"),
    make_option(c("-o", "--outputDir"), type="character",
                default="results",
                help="Path to directory where subpopr results will be stored. \
                Default is \"./results/\"",
                metavar="file path"),
    make_option(c("-p", "--procs"), type="integer",
                default=1,
                help="Number of cores to use for parallel processing. \
                Default is 1.",
                metavar="integer"),
    make_option(c("-s", "--sampleSuffix"), type="character",
                default="",
                help="The constant suffix after the sample names in metaSNV's input bam files. \
                e.g. '.bam' or '.unique.sorted.bam'",
                metavar="string"),
    make_option(c("-a", "--speciesAbundance"), type="character",
                default="doNotRun",
                help="Path to file with species abundances (tsv, optional). \
                Rows are species, columns are samples. Column names must match file \
                names used as metaSNV input (bam files).",
                metavar="file path"),
    make_option(c("-m", "--isMotus"), type="logical",
                default=TRUE,
                help="Is the species abundance profile produced by mOTUs2? (TRUE or FALSE)",
                metavar="logical"),
    make_option(c("-g", "--geneAbundance"), type="character",
                default="doNotRun",
                help="Path to file with gene family abundances (tsv, optional). \
                Species abundances also required for gene correlation. Columns must be named. \
                First column must be named and contain gene family names. Subsequent column names must be sample IDs. \
                Columns must sum to 1.",
                metavar="file path"),
    make_option(c("-d", "--metadata"), type="character",
                default="doNotRun",
                help="Path to file with metadata csv for odds ratio \
                testing (optional)",
                metavar="file path"),
    make_option(c("-n", "--metadataSampleIDCol"), type="character",
                default="sampleID",
                help="Name of column with sample IDs in metadata csv for \
                odds ratio testing (optional)",
                metavar="character"),
    make_option(c("-r", "--createReports"), type="logical",
                default=TRUE,
                help="Whether or not to compile html summary reports (uses Rmarkdown) \
                (TRUE or FALSE). Default is TRUE.",
                metavar="logical"),
    make_option(c("-x", "--fixReadThreshold"), type="numeric",
                default=0.1,
                help="SNV locus filter: max proportion of reads with non-major allele \
                for locus to be considered to be used in defining clusters. (hr)\
                Default is 0.1 (i.e. 90% of reads have major allele)",
                metavar="numeric"),
    make_option(c("-y", "--fixSnvThreshold"), type="numeric",
                default=0.8,
                help="Sample filter: min proportion of SNVs where major \
                allele is sufficently abundant for sample to be used in \
                defining clusters. (hs) Default is 0.8 \
                (i.e. 80% of SNV loci have major allele with frequency > x )",
                metavar="numeric"),
    make_option(c("-z", "--genotypingThreshold"), type="numeric",
                default=0.8,
                help="Genotyping threshold: SNV allele must be more abundand \
                within the cluster by this many percentage points (as decimal <= 1). \
                (gs) Default 0.8",
                metavar="numeric"),
    make_option(c("-q", "--onlyDoSubspeciesDetection"), type="logical",
                default=FALSE,
                help="Whether to only do the first step of the pipeline \
                (just detect the presence/number of subspecies). Default is FALSE. \
                Only intended for troubleshooting.",
                metavar="logical")
  );
  
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);
  
  
}else{
  # TO RUN FROM WITHIN R WITHOUT OPTS ----------
  opt <- list()
  #scriptDir <- "/g/bork3/home/rossum/software/metaSNV2/metaSNV/"
  #setwd("/g/scb2/bork/rossum/metagenomes/human/subspecGeoValidation/all_v2/subpopr_v2")
  #opt$metadata <- "/g/scb2/bork/rossum/metagenomes/human/subspecGeoValidation/all_v2/metadata_allv2.csv"
  #opt$metaSnvResultsDir <- "/g/scb2/bork/rossum/metagenomes/human/subspecGeoValidation/all_v2/metaSNV/outputs_subspec/"
  #opt$speciesAbundance <- "doNotRun"
  #opt$geneAbundance <- "doNotRun"
  
  scriptDir <- "/Users/rossum/Dropbox/PostDocBork/subspecies/toolDevelopment/metaSNV/"
  workDir <- "/Volumes/KESU/scb2/metagenomes/human/subspecGeoValidation/all_v3/"
  setwd(paste0(workDir,"/subpopr"))
  opt$metadata <- paste0(workDir,"metadataForSubspeciesAnalysis.csv")
  opt$metaSnvResultsDir <- paste0(workDir,"/metaSNV/outputs/")
  opt$speciesAbundance <- paste0(workDir,"/motus20/motusForSelectedSamples.tsv") #"doNotRun"
  opt$geneAbundance <- paste0(workDir,"/geneContent/mapToPanGenomes/outputs/counts_unique_norm_tmp1k.tsv")
  opt$sampleSuffix <- ".subspec71.unique.sorted.bam"
  
  opt$procs <- 2
  opt$isMotus <- T
  opt$metadataSampleIDCol <- "sampleID"
  opt$outputDir <- "results_md_s49677-u"
  
  opt$fixReadThreshold <- 0.1
  opt$fixSnvThreshold <- 0.8
  opt$genotypingThreshold <- 0.9
  
  opt$onlyDoSubspeciesDetection<-FALSE
  
}

N.CORES <- opt$procs
SPECIES.ABUNDANCE.PROFILE<-opt$speciesAbundance
SPECIES.ABUND.PROFILE.IS.MOTUS<-opt$isMotus
KEGG.PATH <- opt$geneAbundance
METADATA.PATH <- opt$metadata
METADATA.COL.ID <- opt$metadataSampleIDCol

MAX.PROP.READS.NON.HOMOG <- opt$fixReadThreshold
MIN.PROP.SNV.HOMOG <- opt$fixSnvThreshold
SNV.SUBSPEC.UNIQ.CUTOFF <- opt$genotypingThreshold
CLUSTERING.PS.CUTOFF <- 0.8
DIST.METH.REPORTS <- "mann"

onlyDoSubspeciesDetection<-opt$onlyDoSubspeciesDetection

SAMPLE.ID.SUFFIX <- opt$sampleSuffix

toScreen <- TRUE # if TRUE, lots gets printed to screen, if FALSE, only goes to log file
printProgressBar <- TRUE

if(is.null(opt$metaSnvResultsDir)){
  print_help(opt_parser)
  stop("Path to metaSNV results must be supplied [-i]", call.=FALSE)
}
METASNV.DIR <- opt$metaSnvResultsDir
if (is.null(opt$outputDir)){
  print_help(opt_parser)
  stop("Path to output directory must be supplied [-o]", call.=FALSE)
}
OUT.DIR.BASE <- opt$outputDir

assert0to1<-function(x,nameOfParam){
  if(!is.numeric(x) | x < 0 | x > 1){
    stop("Param \"",nameOfParam ,"\"must be numeric and must be between 0 and 1")
  }
}

assert0to1(opt$fixReadThreshold,"fixReadThreshold")
assert0to1(opt$fixSnvThreshold,"fixSnvThreshold")
assert0to1(opt$genotypingThreshold,"genotypingThreshold")

#source(paste0(scriptDir,"/src/subpopr/inst/metaSNV_subpopr_SETTINGS.R"))
DIST.METH.REPORTS="mann" # what method to use for generating reports: either "mann" or "allele"
BAMS.TO.USE = NULL # default = NULL
ANALYSE.ALLELE.DISTANCES = F
USE.PACKAGE.PREDICTION.STRENGTH = FALSE # default = FALSE


SUBPOPR.DIR<-paste0(scriptDir,"/src/subpopr/")

SUBPOPR_RESULTS_DIR=paste0(OUT.DIR.BASE,"/params",
                           ".hr",MAX.PROP.READS.NON.HOMOG*100,
                           ".hs",MIN.PROP.SNV.HOMOG*100,
                           ".ps",CLUSTERING.PS.CUTOFF*100,
                           ".gs",SNV.SUBSPEC.UNIQ.CUTOFF*100,"/")
OUT.DIR=paste0(SUBPOPR_RESULTS_DIR,"/",basename(METASNV.DIR),"/")
dir.create(OUT.DIR, recursive = T, showWarnings = FALSE)

# Check input file existances ---------------------------------------

checkFile <- function(path, fileTypeName){
  if(!is.null(path) && !is.na(path) &&
     nchar(path) > 0 && path != "doNotRun"){
    if(!file.exists(path)){
      stop(fileTypeName, " file specified but does not exist: ",
           path)
    }
  }
}

checkFile(SPECIES.ABUNDANCE.PROFILE,"Species abundance")
checkFile(METADATA.PATH,"Metadata")
checkFile(KEGG.PATH, "Gene family abundance")
checkFile(METASNV.DIR, "MetaSNV output directory")

# Logging set up -----------------------------------------------------

logFile <- paste0(OUT.DIR,"/log.txt")
print(paste("Logging to:",logFile))

capture.output(print("Command was --------------------------------------------------"),
               file = logFile,append = FALSE)
capture.output(paste(commandArgs(trailingOnly = FALSE),collapse = " "),
               file = logFile,append = TRUE)

capture.output(print("Variable values --------------------------------------------------"),
               file = logFile,append = TRUE)
rm(option_list)
capture.output(ls.str(),file = logFile,
               append = TRUE) # print all variables (and values for strings)

capture.output(print("Run output --------------------------------------------------"),
               file = logFile,append = TRUE)
sink(file = logFile, append = TRUE,
     type = c("output", "message"), split = toScreen)




# Load library dependencies -------------------------------------------
print("Loading R libraries...")

#if(!is.null(LIB.DIR) && dir.exists(LIB.DIR)){
#  .libPaths(c(LIB.DIR))
#  print(paste0("Using R library directories:",paste(.libPaths(),collapse=" : ")))
#}

# REQUIRES CAIRO TO BE INSTALLED, EITHER THROUGH 'install.packages()' OR THROUGH 'conda install -c anaconda cairo'
# requires pandoc
suppressPackageStartupMessages(library(Cairo))
suppressPackageStartupMessages(library(fpc))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(kableExtra)) # to do: remove this from package
suppressPackageStartupMessages(library(rmarkdown)) # for report rendering

suppressPackageStartupMessages(library(coin)) # only used in phenotype assoc test part -- remove?
suppressPackageStartupMessages(library(questionr)) # only used in phenotype assoc test part -- remove?
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(batchtools))

#Error: pandoc version 1.12.3 or higher is required and was not found (see the help page ?rmarkdown::pandoc_available).
# throw and error if the required version of pandoc is not found
if(makeReports & !rmarkdown::pandoc_available(version = "1.12.3",error = F)){
  warning("Reports will not be generated because software 'pandoc' version 1.12.3 or higher",
          " is required and was not found. Recommended action: abort now, install/update pandoc",
          " https://pandoc.org/ ",
          " and restart script.")
  makeReports <- FALSE
}


# Load subpopr files --------------------------------------------------------------

print(paste0("Loading subpopr source files from ",SUBPOPR.DIR))

srcFiles <- list.files(path = paste0(SUBPOPR.DIR,"/R"),pattern = "*.R$",
                       full.names = TRUE,recursive = TRUE,ignore.case = TRUE)
if(length(srcFiles) < 1){
  stop("Source files for subpopr not found in:", SUBPOPR.DIR,
       " . Expected subdirectory called \"R\" with files matching: ",
       paste0(SUBPOPR.DIR,"/R/*.R"))
}

tmp <- sapply(srcFiles,source,echo = F,local = .GlobalEnv)
rmdDir <- paste0(SUBPOPR.DIR,"/inst/rmd/")
pyScriptDir <- paste0(SUBPOPR.DIR,"/inst/")
TAXA.NCBI.MOTU.MAP <- readRDS(paste0(SUBPOPR.DIR,"/data/TAXA.NCBI.MOTU.MAP.Rds"))

# Get species to analyse ----------------------------------------------------------------

print(paste0("Using files from ",METASNV.DIR,"/distances and ",METASNV.DIR,"/filtered/pop/"))

# get species IDs to analyse
specDist <- list.files(path=paste0(METASNV.DIR,"/distances"),pattern = '.*mann.dist$',full.names = T)
specSnpFreq <- list.files(path=paste0(METASNV.DIR,"/filtered/pop/"),pattern = '.*filtered.freq$',full.names = T)

specDist <- sub(basename(specDist) ,pattern = "\\..*$",replacement = "")
specSnpFreq <- sub(basename(specSnpFreq) ,pattern = "\\..*$",replacement = "")

species <- intersect(specDist,specSnpFreq)

if(!identical(specDist,specSnpFreq)){
  warning(paste("Species indicated in",
                paste0(METASNV.DIR,"/distances"),
                " (",length(specDist)," species) and",
                paste0(METASNV.DIR,"/filtered/pop/")," (",length(specSnpFreq),
                " species) are not identical. Using the intersection (",
                length(species),"species).Species not found in both locations: ",
                paste(setdiff(specDist,specSnpFreq),collapse = ",")))
  
}

if(length(specDist) == 0){
  stop(paste0("No appropriate files found in ", METASNV.DIR))
}

#species <- 657318 #537011 #"657317"


# Set up parallel processing ----------------------------------------

ncoresUsing <- min(N.CORES,length(species))

bpParam <- MulticoreParam(workers = min(N.CORES,length(species)),
                          jobname = "subpopr",
                          stop.on.error = FALSE,
                          threshold = "DEBUG",
                          log = TRUE,
                          progressbar = printProgressBar,
                          logdir = paste0(OUT.DIR,"/threadLogs"))
dir.create(paste0(OUT.DIR,"/threadLogs"), recursive = T, showWarnings = FALSE)

print(paste("Running subpopr on",length(species),"species using",ncoresUsing,"cores."))
print("Progress bar reflects the percentage of species analysed. Progression in time will not be linear.")

printBpError <- function(result){
  if(all(bpok(result))){
    return("") # blank prints "NULL"
  }else{
    print( paste("Error in ",length(which(!bpok(result)))," task(s)",
                 #paste0(which(!bpok(result)),collapse = ","),
                 ", see files with text \"Success: FALSE\" in ",paste0(OUT.DIR,"/threadLogs") ))
    #print(result[[which(!bpok(result))]])
    #above creates this: Error in result[[which(!bpok(result))]] : subscript out of bounds
  }
}

# Try to find subspecies #######################################################################
# (substructure/clustering) within species

runDefine <- function(spec){
  print("")
  flog.info("=== Assessing presence of subspecies in species %s ===", spec)
  #cat(dput(spec), file = paste0("logFile_", spec, ".txt"))
  defineSubpopulations(spec, metaSNVdir = METASNV.DIR, outDir = OUT.DIR,
                       maxPropReadsNonHomog = MAX.PROP.READS.NON.HOMOG,
                       minPropHomogSnvAllelesPerSample = MIN.PROP.SNV.HOMOG,
                       psCut = CLUSTERING.PS.CUTOFF,
                       uniqSubpopSnvFreqThreshold = SNV.SUBSPEC.UNIQ.CUTOFF,
                       bamFileNamesToUsePath = BAMS.TO.USE,
                       usePackagePredStrength = USE.PACKAGE.PREDICTION.STRENGTH)
}

if(!useExistingClustering){
  resultsPerSpecies <- BiocParallel::bptry(
    BiocParallel::bplapply(species, runDefine, BPPARAM = bpParam))
  printBpError(resultsPerSpecies)
  
  # this is error prone and output is not used anyway
  # resultsPerSpeciesFixed <- resultsPerSpecies
  # if(any(!bpok(resultsPerSpeciesFixed))){
  #   resultsPerSpeciesFixed[[which(!bpok(resultsPerSpeciesFixed))]] <- "Error"
  # }
  # resultsPerSpeciesDF <- cbind.data.frame(SpeciesID=species,
  #                                         ClusteringResult=unlist(resultsPerSpeciesFixed))
  # write.csv(x = resultsPerSpeciesDF,file = paste0(OUT.DIR,"/log_clusteringSummaryPerSpecies.csv"))
  saveRDS(resultsPerSpecies, file = paste0(OUT.DIR,"/log_clusteringSummaryPerSpecies.rds"))
}
# summarise the results from clustering
summariseClusteringResultsForAll(OUT.DIR,distMeth="mann")

allSubstruc <- list.files(path=OUT.DIR,
                          pattern = '_hap_out\\.txt$',full.names = T)
allSubstrucSpecies <- unique(sub(basename(allSubstruc) ,
                                 pattern = "_hap_out\\.txt$",replacement = ""))

print(paste0("Species with substructure: ",
             length(allSubstrucSpecies),"/",length(species)))

if(onlyDoSubspeciesDetection){
  combineAllSummaries(OUT.DIR)
  stop("Subpopr stopped due to 'onlyDoSubspeciesDetection' flag being true")
}

# Handle species with no subspecies #####################################################################
# for those species that did not cluster, generate a report so we can look into why

# get all species where no potential cluster medoids could be defined
noSubstruc1dir <- getClustMedoidDefnFailedDir(OUT.DIR)
noSubstruc1 <- list.files(path=noSubstruc1dir,
                          pattern = paste0(DIST.METH.REPORTS ,
                                           '_distMatrixUsedForClustMedoidDefns\\.txt$'),
                          full.names = T)
noSubstrucSpecies <- unique(sub(basename(noSubstruc1) ,
                                pattern = paste0("_",DIST.METH.REPORTS ,
                                                 "_distMatrixUsedForClustMedoidDefns\\.txt"),
                                replacement = ""))

if(makeReports){
  print("Compiling reports for species without clusters due to centroid failure")
  tmp <- BiocParallel::bptry(
    BiocParallel::bplapply(noSubstrucSpecies, BPPARAM = bpParam,
                           renderDetailedSpeciesReport,
                           metasnvOutDir = METASNV.DIR,
                           distMethod = DIST.METH.REPORTS ,
                           subpopOutDir = noSubstruc1dir,
                           bamSuffix = SAMPLE.ID.SUFFIX,
                           rmdDir = rmdDir ))
  printBpError(tmp)
}
# get all species where cluster medoids could be defined
# but clusters were not significant (PS values < threshold)
noSubstruc2dir <- getNoClusteringDir(OUT.DIR)
noSubstruc2 <- list.files(path=noSubstruc2dir,
                          pattern = paste0(DIST.METH.REPORTS ,'_distMatrixUsedForClustMedoidDefns\\.txt$'),
                          full.names = T)
noSubstrucSpecies <- unique(sub(basename(noSubstruc2) ,
                                pattern = paste0("_",DIST.METH.REPORTS ,"_distMatrixUsedForClustMedoidDefns\\.txt"),
                                replacement = ""))
if(makeReports){
  print("Compiling reports for species without clusters")
  tmp <- BiocParallel::bptry(
    BiocParallel::bplapply(noSubstrucSpecies,
                           renderDetailedSpeciesReport,
                           metasnvOutDir = METASNV.DIR,
                           distMethod = DIST.METH.REPORTS ,
                           subpopOutDir = noSubstruc2dir,
                           bamSuffix = SAMPLE.ID.SUFFIX,
                           rmdDir = rmdDir ))
  
  printBpError(tmp)
}
# Handle species with subspecies #######################################################################

# continue processing those species that could be used to define subspecies

# get all species with clustering/substructure

if(length(allSubstrucSpecies) == 0){
  stop(paste0("Substructure not detected in any species (",
              length(species)," tested). Aborting."))
}

# Genotype clusters #####################################################
# Try to determine genotypes for each cluster
# Use these genotypes to:
# 1) detect clusters in more samples
# 2) get abundances of these genotypes per sample (~subspecies abundace)

print("Genotyping clusters")
print("Identifying genotyping SNVs")

if(useExistingExtension){
  print("Using previously made .pos files and extension results")
}else{
  # creates *.pos files
  x <- tryCatch(expr =
                  pyGetPlacingRelevantSubset(outDir=OUT.DIR,
                                             metaSnvDir=METASNV.DIR,
                                             scriptDir = pyScriptDir),
                error = function(e){
                  print(paste("ERROR: ",e$message ))
                  print("Skipping subspecies genotyping.")
                }
  )
  
  
  # get all posFiles
  allPos <- list.files(path=OUT.DIR,pattern = '.*_.\\.pos$',full.names = T)
  doExtension<-TRUE
  if(length(allPos) == 0){
    warning("Genotyping failed. No *.pos files found. Not genotyping subspecies.")
    doExtension <- FALSE
  }
  
  if(doExtension){
    print("Compiling genotyping SNVs")
    #tmp <- foreach(pos=allPos) %dopar% pyConvertSNPtoAllelTable(posFile = pos)
    tmp <- BiocParallel::bptry(
      BiocParallel::bplapply(allPos, BPPARAM = bpParam,
                             pyConvertSNPtoAllelTable,
                             scriptDir = pyScriptDir))
    
    printBpError(tmp)
    
    print("Determining abundance of clusters using genotyping SNVs")
    #tmp <- foreach(spec=allSubstrucSpecies) %dopar% useGenotypesToProfileSubpops(spec, metaSNVdir=METASNV.DIR, outDir=OUT.DIR )
    tmp <- BiocParallel::bptry(
      BiocParallel::bplapply(allSubstrucSpecies, BPPARAM = bpParam,
                             useGenotypesToProfileSubpops,
                             metaSNVdir=METASNV.DIR,
                             outDir=OUT.DIR ))
    
    printBpError(tmp)
    
    summariseClusteringExtensionResultsForAll(resultsDir=OUT.DIR,distMeth="mann")
    
  }
}
# Compile detailed reports for species with subspecies/clusters --------------

runRend <- function(spec){
  print("")
  flog.info("Rendering report for species %s", spec)
  renderDetailedSpeciesReport(speciesID = spec,
                              subpopOutDir = OUT.DIR,
                              metasnvOutDir = METASNV.DIR,
                              distMethod = DIST.METH.REPORTS ,
                              bamSuffix = SAMPLE.ID.SUFFIX,
                              rmdDir = rmdDir)
}

if(makeReports){
  print("Compiling reports for species with clusters")
  tmp <- BiocParallel::bptry(
    BiocParallel::bplapply(allSubstrucSpecies,
                           BPPARAM = bpParam,
                           runRend))
  # if failed, try again...often it's just a timing conflict error from parallelising
  if(!all(bpok(tmp))){
    print(paste("Retrying compilation of",length(which(!bpok(tmp))),"failed reports"))
    tmp <- BiocParallel::bptry(
      BiocParallel::bplapply(X = allSubstrucSpecies,
                             BPREDO=tmp,
                             BPPARAM = bpParam,
                             runRend))
  }
  printBpError(tmp)
}

speciesToAssess <- list.files(path=OUT.DIR,pattern = '.*_extended_clustering_wFreq.tab$',full.names = F) %>%
  sub(pattern = "_extended_clustering_wFreq.tab",replacement = "")
if(length(speciesToAssess)>0){
  subpopFreqSumsStats <- assessSubpopCompleteness(speciesToAssess,subpoprOutDir = OUT.DIR)
  write.table(subpopFreqSumsStats,file=paste0(OUT.DIR,"/subpopFreqSumsStats.tsv"),sep = "\t",row.names = F,quote = F)
}

# Get subspecies abundances relative to whole community ---------------------------------------

if(doSpecies && !is.null(SPECIES.ABUNDANCE.PROFILE) &&
   SPECIES.ABUNDANCE.PROFILE != "doNotRun" &&
   file.exists(SPECIES.ABUNDANCE.PROFILE)){
  print("Calculating cluster abundances using species abundances...")
  
  tmp <- BiocParallel::bptry(
    BiocParallel::bplapply(allSubstrucSpecies, BPPARAM = bpParam,
                           useSpeciesAbundToCalcSubspeciesAbund,
                           speciesAbundanceProfileFilePath=SPECIES.ABUNDANCE.PROFILE,
                           outDir=OUT.DIR,
                           speciesProfileIsMotus = SPECIES.ABUND.PROFILE.IS.MOTUS))
  
  printBpError(tmp)
  abunds <- collectSubpopAbunds(OUT.DIR)
  if(is.null(abunds)){
    warning("Subspecies abundance calculations failed. ",
            "No expected results files exist (.*hap_coverage_extended_normed.tab). ",
            "See log files.")
  }
  
}else if(SPECIES.ABUNDANCE.PROFILE != "doNotRun"){
  print(paste0("Not running species abundance analysis.",
               " Required file not specified or does not exist: ",
               SPECIES.ABUNDANCE.PROFILE))
}



# Test metadata associations ##########
if(!is.null(METADATA.PATH) && file.exists(METADATA.PATH)){
  if(makeReports){
    
    doRendMd <- function(spec){
      flog.info("Rendering metadata association report for species %s", spec)
      renderTestPhenotypeAssocReport(speciesID = spec,
                                     subpopOutDir = OUT.DIR,
                                     categoryColumnNames = METADATA.COLS.TO.TEST, #"status",
                                     sampleIDColumnName = METADATA.COL.ID, #"ID",
                                     sampleExtension = SAMPLE.ID.SUFFIX, #".ULRepGenomesv11.unique.sorted.bam",
                                     metadataFile = METADATA.PATH,
                                     rmdDir = rmdDir)
    }
    print("Associating with metadata...")
    tmp <- BiocParallel::bptry(
      BiocParallel::bplapply(allSubstrucSpecies,
                             BPPARAM = bpParam,
                             doRendMd))
    # if failed, try again...often it's just a timing conflict error from parallelising
    if(!all(bpok(tmp))){
      print(paste("Retrying compilation of",length(which(!bpok(tmp))),"failed reports"))
      tmp <- BiocParallel::bptry(
        BiocParallel::bplapply(X = allSubstrucSpecies,
                               BPREDO=tmp,
                               BPPARAM = bpParam,
                               doRendMd))
    }
    printBpError(tmp)
  }
  summariseMetadataAssocResultsForAll(OUT.DIR)
}else if(METADATA.PATH != "doNotRun"){
  print(paste0("Not running phenotype/metadata association analysis.",
               " Required file not specified or does not exist: ",
               METADATA.PATH))
}

# Test for gene correlations ##########
if(!is.null(KEGG.PATH) && file.exists(KEGG.PATH) &&
   !is.null(SPECIES.ABUNDANCE.PROFILE) && #species abundances required for correlation with gene abundaces
   SPECIES.ABUNDANCE.PROFILE != "doNotRun" &&
   file.exists(SPECIES.ABUNDANCE.PROFILE)){
  print(paste("Testing for gene correlations for",length(allSubstrucSpecies),
              "species using",ncoresUsing,"cores"))
  
  geneFamilyType<-"Genes"
  # print("Correlating cluster and gene family abundances (Pearson & Spearman)...")
  # #pearson
  # tmp <- BiocParallel::bptry(
  #   BiocParallel::bplapply(allSubstrucSpecies, BPPARAM = bpParam,
  #                          correlateSubpopProfileWithGeneProfiles,
  #                          OUT.DIR,KEGG.PATH,
  #                          geneFamilyType=geneFamilyType))
  # 
  # # if failed, try again...often it's just a timing conflict error from parallelising
  # if(!all(bpok(tmp))){
  #   print(paste("Retrying",length(which(!bpok(tmp)))," failed computation of correlations..."))
  #   tmp <- BiocParallel::bptry(
  #     BiocParallel::bplapply(allSubstrucSpecies, BPPARAM = bpParam, #SerialParam() could redo with serial param if memory is issue
  #                            BPREDO=tmp,
  #                            correlateSubpopProfileWithGeneProfiles,
  #                            OUT.DIR,KEGG.PATH,
  #                            geneFamilyType=geneFamilyType))
  # }
  # printBpError(tmp)
  
  if(makeGeneReports){ #makeReports){
    
    print("Compiling gene content reports...")
    tmp <- BiocParallel::bptry(BiocParallel::bplapply(allSubstrucSpecies,
                                                      BPPARAM = bpParam, #SerialParam(), # for some reason, parallel fails here
                                                      renderGeneContentReport,
                                                      subpopOutDir = OUT.DIR,
                                                      geneFamilyAbundancesFile = KEGG.PATH,
                                                      bamSuffix= SAMPLE.ID.SUFFIX,
                                                      rmdDir = rmdDir))
    
    # if failed, try again...often it's just a timing conflict error from parallelising
    if(!all(bpok(tmp))){
      print(paste("Retrying compilation of",length(which(!bpok(tmp)))," failed gene content reports"))
      tmp <- BiocParallel::bptry(
        BiocParallel::bplapply(BPREDO=tmp,
                               allSubstrucSpecies,
                               BPPARAM = bpParam, #SerialParam(), # for some reason, parallel fails here
                               renderGeneContentReport,
                               subpopOutDir = OUT.DIR,
                               geneFamilyAbundancesFile = KEGG.PATH,
                               bamSuffix= SAMPLE.ID.SUFFIX,
                               rmdDir = rmdDir))
    }
    printBpError(tmp)
  }
  summariseGeneFamilyCorrelationResultsForAll(OUT.DIR,geneFamilyType)
}else if(KEGG.PATH != "doNotRun"){
  print(paste0("Not running gene content analysis.",
               " Required file not specified or does not exist: ",
               KEGG.PATH))
}


# Summarise results ##########
print("Summarising results...")
combineAllSummaries(OUT.DIR)
if(makeReports){
  renderResultsSummaryReport(OUT.DIR,rmdDir = rmdDir)
  print(paste0("Results summarised, see: ",OUT.DIR,"/resultsSummary.html"))
}else{
  print(paste0("Results summarised, see: ",OUT.DIR,"/summary_allResults.csv"))
}

print("Subpopr finished.")
print(proc.time() - ptm)

