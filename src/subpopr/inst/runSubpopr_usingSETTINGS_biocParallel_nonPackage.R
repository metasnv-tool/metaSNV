#!/usr/bin/env Rscript

# qlogin -pe smp 12 -l h_vmem=10G # 10G memory per core is plenty
# REQUIRES PYTHON (2 or 3)
# nohup Rscript runSubpopr.R

# mkdir subpoprLocalTest
# cp -r ~/Dropbox/PostDocBork/subspecies/toolDevelopment/subpopr/R ./subpoprLocalTest/
# cp -r ~/Dropbox/PostDocBork/subspecies/toolDevelopment/subpopr/inst ./subpoprLocalTest/
# cp -r ~/Dropbox/PostDocBork/subspecies/toolDevelopment/subpopr/data ./subpoprLocalTest/


# clear the environment
rm(list=ls())

# PARSE PARAMS -------------------------------------------------------------

library(getopt)
library(optparse)

option_list = list(
  #make_option(c("-s", "--settings"), type="character", default="SETTINGS.R",  # default=NULL,
  #            help="Settings file, default is %default in current directory", metavar="character")
  make_option(c("-s", "--settings"), type="character", default=NULL,
              help="Settings file", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$settings)){
  print_help(opt_parser)
  stop("Settings file must be supplied [-s]", call.=FALSE)
}

# LOAD SETTINGS ------------------------------------------------------------

# load settings for this run of subpopr
settingsFilePath<-opt$settings
if(!file.exists(settingsFilePath)){
  stop(paste("Error. Required settings file does not exist:",settingsFilePath))
}

source(settingsFilePath)

dir.create(OUT.DIR, recursive = TRUE, showWarnings = FALSE)
logFile <- paste0(OUT.DIR,"/log.txt")
sink(file = logFile, append = FALSE, type = c("output", "message"), split = TRUE)

makeReports <- TRUE

ls.str()  # print all variables (and values for strings)



# LOAD LIBRARY DEPENDENCIES -------------------------------------------

if(!is.null(LIB.DIR) && dir.exists(LIB.DIR)){
  .libPaths(c(LIB.DIR))
  print(paste0("Using R library directories:",paste(.libPaths(),collapse=" : ")))
}

#install.packages(c("fpc", "ape", "ggplot2", "gridExtra", "cluster","doParallel","readr","dplyr","tidyr","data.table"))
#install.packages('BiocManager')
#BiocManager::install(c("BiocParallel", "BatchJobs")) #batchtools ?

#install.packages("./doParallel_1.0.11.tgz", repos = NULL, type = .Platform$pkgType)
#install.packages("doParallel")

# subpopr dependnecies

# REQUIRES CAIRO TO BE INSTALLED, EITHER THROUGH 'install.packages()' OR THROUGH 'conda install -c anaconda cairo'
# requires pandoc
library(fpc)
library(ape)
library(ggplot2)
library(gridExtra)
library(cluster)
library(dplyr)
library(tidyr)
library(readr)
library(data.table)
library(kableExtra) # to do: remove this from package
library(rmarkdown) # for report rendering
#Error: pandoc version 1.12.3 or higher is required and was not found (see the help page ?rmarkdown::pandoc_available).
# throw and error if the required version of pandoc is not found
if(makeReports & !rmarkdown::pandoc_available(version = "1.12.3",error = F)){
  warning("Reports will not be generated because software 'pandoc' version 1.12.3 or higher",
          " is required and was not found. Recommended action: abort now, install/update pandoc",
          " https://pandoc.org/ ",
          " and restart script.")
  makeReports <- FALSE
}

library(ggrepel) # used?

library(doParallel) # aborts everything if error is encountered
library(BiocParallel) # doesn't die on error
library(batchtools)

library(coin) # only used in phenotype assoc test part -- remove?
library(questionr) # only used in phenotype assoc test part -- remove?

library(data.table) # only used in gene content corr calc -- remove?


# Load subpopr files --------------------------------------------------------------

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

ptm <- proc.time()

# Set up parallel processing ----------------------------------------

ncoresUsing <- min(N.CORES,length(species))
#cl <- makeCluster(min(N.CORES,length(species)), outfile=logFile)
#registerDoParallel(cl)  #alt: registerDoParallel(cores=2) crashes silently when trying to write png()
#ncoresUsing <- getDoParWorkers()

bpParam <- MulticoreParam(workers = min(N.CORES,length(species)),
                          stop.on.error = FALSE,
                          threshold = "DEBUG",
                          log = TRUE,
                          logdir = "logs")
dir.create("logs", recursive = T, showWarnings = FALSE)

print(paste("Running subpopr on",length(species),"species using",ncoresUsing,"cores"))

printBpError <- function(result){
  if(all(bpok(result))){
    return(NULL)
  }else{
    print("Error:")
    tail(attr(result[[which(!bpok(result))]], "traceback"))
  }
}

# Try to find subspecies #######################################################################
# (substructure/clustering) within species

#resultsPerSpecies <- foreach(spec=species, .combine='c') %dopar% {
#  print(spec)
#  cat(dput(spec), file = paste0("logFile_", spec, ".txt"))
#  defineSubpopulations(spec, metaSNVdir = METASNV.DIR, outDir = OUT.DIR,
#                                maxPropReadsNonHomog = MAX.PROP.READS.NON.HOMOG,
#                                minPropHomogSnvAllelesPerSample = MIN.PROP.SNV.HOMOG,
#                                psCut = CLUSTERING.PS.CUTOFF,
#                                uniqSubpopSnvFreqThreshold = SNV.SUBSPEC.UNIQ.CUTOFF,
#                                bamFileNamesToUsePath = BAMS.TO.USE)
#}
#resultsPerSpeciesDF <- cbind.data.frame(species, resultsPerSpecies)
#colnames(resultsPerSpeciesDF) <- c("SpeciesID","ClusteringResult")
#write.csv(x = resultsPerSpeciesDF,file = paste0(OUT.DIR,"/log_clusteringSummaryPerSpecies.csv"))
#quote = T, sep = ",",col.names = T,row.names = F)

runDefine <- function(spec){
  print(spec)
  #cat(dput(spec), file = paste0("logFile_", spec, ".txt"))
  defineSubpopulations(spec, metaSNVdir = METASNV.DIR, outDir = OUT.DIR,
                       maxPropReadsNonHomog = MAX.PROP.READS.NON.HOMOG,
                       minPropHomogSnvAllelesPerSample = MIN.PROP.SNV.HOMOG,
                       psCut = CLUSTERING.PS.CUTOFF,
                       uniqSubpopSnvFreqThreshold = SNV.SUBSPEC.UNIQ.CUTOFF,
                       bamFileNamesToUsePath = BAMS.TO.USE,
                       usePackagePredStrength = USE.PACKAGE.PREDICTION.STRENGTH)
}

resultsPerSpecies <- BiocParallel::bptry(BiocParallel::bplapply(species, runDefine, BPPARAM = bpParam))
resultsPerSpeciesDF <- cbind.data.frame(SpeciesID=species, ClusteringResult=unlist(resultsPerSpecies))
write.csv(x = resultsPerSpeciesDF,file = paste0(OUT.DIR,"/log_clusteringSummaryPerSpecies.csv"))

# summarise the results from clustering
summariseClusteringResultsForAll(OUT.DIR,distMeth="mann")

# Handle species with no subspecies #####################################################################
# for those species that did not cluster, generate a report so we can look into why

# get all species where no potential cluster medoids could be defined
noSubstruc1dir <- getClustMedoidDefnFailedDir(OUT.DIR)
noSubstruc1 <- list.files(path=noSubstruc1dir,pattern = paste0(DIST.METH.REPORTS ,'_distMatrixUsedForClustMedoidDefns\\.txt$'),full.names = T)
noSubstrucSpecies <- unique(sub(basename(noSubstruc1) ,pattern = paste0("_",DIST.METH.REPORTS ,"_distMatrixUsedForClustMedoidDefns\\.txt"),replacement = ""))

if(makeReports){
  # make report
  # tmp <- foreach(spec=noSubstrucSpecies) %dopar%
  #   renderDetailedSpeciesReport(speciesID = spec,metasnvOutDir = METASNV.DIR, distMethod = DIST.METH.REPORTS ,
  #                                        subpopOutDir = noSubstruc1dir,
  #                                        bamSuffix = SAMPLE.ID.SUFFIX,
  #                                       rmdDir = rmdDir)

  tmp <- BiocParallel::bptry(BiocParallel::bplapply(noSubstrucSpecies, BPPARAM = bpParam,
                                                    renderDetailedSpeciesReport,
                                                    metasnvOutDir = METASNV.DIR,
                                                    distMethod = DIST.METH.REPORTS ,
                                                    subpopOutDir = noSubstruc1dir,
                                                    bamSuffix = SAMPLE.ID.SUFFIX,
                                                    rmdDir = rmdDir ))
  printBpError(tmp)
}
# get all species where cluster medoids could be defined but clusters were not significant (PS values < threshold)
noSubstruc2dir <- getNoClusteringDir(OUT.DIR)
noSubstruc2 <- list.files(path=noSubstruc2dir,pattern = paste0(DIST.METH.REPORTS ,'_distMatrixUsedForClustMedoidDefns\\.txt$'),full.names = T)
noSubstrucSpecies <- unique(sub(basename(noSubstruc2) ,pattern = paste0("_",DIST.METH.REPORTS ,"_distMatrixUsedForClustMedoidDefns\\.txt"),replacement = ""))
if(makeReports){
  # tmp <- foreach(spec=noSubstrucSpecies) %dopar%
  #   renderDetailedSpeciesReport(speciesID = spec,metasnvOutDir = METASNV.DIR, distMethod = DIST.METH.REPORTS ,
  #                                        subpopOutDir = noSubstruc2dir,
  #                                        bamSuffix = SAMPLE.ID.SUFFIX,
  #                                       rmdDir = rmdDir)

  tmp <- BiocParallel::bptry(BiocParallel::bplapply(noSubstrucSpecies,
                                                    renderDetailedSpeciesReport,
                                                    metasnvOutDir = METASNV.DIR,
                                                    distMethod = DIST.METH.REPORTS ,
                                                    subpopOutDir = noSubstruc2dir,
                                                    bamSuffix = SAMPLE.ID.SUFFIX,
                                                    rmdDir = rmdDir ),
                             rmdDir = rmdDir)
  
  printBpError(tmp)
}
# Handle species with subspecies #######################################################################

# continue processing those species that could be used to define subspecies

# get all species with clustering/substructure
#allSubstruc <- list.files(path=OUT.DIR,pattern = '_[:digit:]+_hap_positions\\.tab$',full.names = T)
#allSubstrucSpecies <- unique(sub(basename(allSubstruc) ,pattern = "_[:digit:]+_hap_positions\\.tab$",replacement = ""))
allSubstruc <- list.files(path=OUT.DIR,pattern = '_hap_out\\.txt$',full.names = T)
allSubstrucSpecies <- unique(sub(basename(allSubstruc) ,pattern = "_hap_out\\.txt$",replacement = ""))


print(paste0("Species with substructure: ",length(allSubstrucSpecies),"/",length(species)))

if(length(allSubstrucSpecies) == 0){
  stopCluster(cl)
  stop(paste0("Substructure not detected in any species (",
              length(species)," tested). Aborting."))
}

# Genotype clusters #####################################################
# Try to determine genotypes for each cluster
# Use these genotypes to:
# 1) detect clusters in more samples
# 2) get abundances of these genotypes per sample (~subspecies abundace)

pyGetPlacingRelevantSubset(outDir=OUT.DIR, metaSnvDir=METASNV.DIR,scriptDir = pyScriptDir)

# get all posFiles
allPos <- list.files(path=OUT.DIR,pattern = '.*_.\\.pos$',full.names = T)

#tmp <- foreach(pos=allPos) %dopar% pyConvertSNPtoAllelTable(posFile = pos)
tmp <- BiocParallel::bptry(BiocParallel::bplapply(allPos, BPPARAM = bpParam,
                                                  pyConvertSNPtoAllelTable,
                                                  scriptDir = pyScriptDir))

printBpError(tmp)

#tmp <- foreach(spec=allSubstrucSpecies) %dopar% useGenotypesToProfileSubpops(spec, metaSNVdir=METASNV.DIR, outDir=OUT.DIR )
tmp <- BiocParallel::bptry(BiocParallel::bplapply(allSubstrucSpecies, BPPARAM = bpParam,
                                                  useGenotypesToProfileSubpops,
                                                  metaSNVdir=METASNV.DIR,
                                                  outDir=OUT.DIR ))

printBpError(tmp)

summariseClusteringExtensionResultsForAll(resultsDir=OUT.DIR,distMeth="mann")

if(makeReports){
  # Render reports for subspecies definitions ###########
  # tmp <- foreach(spec=allSubstrucSpecies) %dopar%
  #   renderDetailedSpeciesReport(speciesID = spec,
  #                                        subpopOutDir = OUT.DIR,
  #                                        metasnvOutDir = METASNV.DIR,
  #                                        distMethod = DIST.METH.REPORTS ,
  #                                        bamSuffix = SAMPLE.ID.SUFFIX,
  #                                       rmdDir = rmdDir)

  tmp <- BiocParallel::bptry(BiocParallel::bplapply(allSubstrucSpecies, BPPARAM = bpParam,
                                                    renderDetailedSpeciesReport,
                                                    subpopOutDir = OUT.DIR,
                                                    metasnvOutDir = METASNV.DIR,
                                                    distMethod = DIST.METH.REPORTS ,
                                                    bamSuffix = SAMPLE.ID.SUFFIX,
                                                    rmdDir = rmdDir))
  
  printBpError(tmp)
}
speciesToAssess <- list.files(path=OUT.DIR,pattern = '.*_extended_clustering_wFreq.tab$',full.names = F) %>%
  sub(pattern = "_extended_clustering_wFreq.tab",replacement = "")
if(length(speciesToAssess)>0){
  subpopFreqSumsStats <- assessSubpopCompleteness(speciesToAssess,subpoprOutDir = OUT.DIR)
  write.table(subpopFreqSumsStats,file=paste0(OUT.DIR,"/subpopFreqSumsStats.tsv"),sep = "\t",row.names = F,quote = F)
}

# Get subspecies abundances relative to whole community ---------------------------------------
if(file.exists(SPECIES.ABUNDANCE.PROFILE)){
  # get relative abundances of subspecies *across* species
  # (based on relative abundance within species and species abundances)
  # tmp <- foreach(spec=allSubstrucSpecies) %dopar%
  #   useSpeciesAbundToCalcSubspeciesAbund(spec,motuProfileFilePath=SPECIES.ABUNDANCE.PROFILE, outDir=OUT.DIR )
  tmp <- BiocParallel::bptry(BiocParallel::bplapply(allSubstrucSpecies, BPPARAM = bpParam,
                                                    useSpeciesAbundToCalcSubspeciesAbund,
                                                    speciesAbundanceProfileFilePath=SPECIES.ABUNDANCE.PROFILE,
                                                    outDir=OUT.DIR,
                                                    speciesProfileIsMotus = SPECIES.ABUND.PROFILE.IS.MOTUS))
  
  printBpError(tmp)
  abunds <- collectSubpopAbunds(OUT.DIR)

}else{
  warning(paste0("Not running species abundance analysis. Required file does not exist: ",SPECIES.ABUNDANCE.PROFILE))
}



# Test metadata associations ##########
if(file.exists(METADATA.PATH)){
  if(makeReports){
    # tmp <- foreach(spec=allSubstrucSpecies) %dopar%
    #   renderTestPhenotypeAssocReport(speciesID = spec,
    #                                           subpopOutDir = OUT.DIR,
    #                                           categoryColumnNames = METADATA.COLS.TO.TEST, #"status",
    #                                           sampleIDColumnName = METADATA.COL.ID, #"ID",
    #                                           sampleExtension = SAMPLE.ID.SUFFIX, #".ULRepGenomesv11.unique.sorted.bam",
    #                                           metadataFile = METADATA.PATH,
    #                                           rmdDir = rmdDir)


    tmp <- BiocParallel::bptry(BiocParallel::bplapply(allSubstrucSpecies, BPPARAM = bpParam,
                                                      renderTestPhenotypeAssocReport,
                                                      subpopOutDir = OUT.DIR,
                                                      categoryColumnNames = METADATA.COLS.TO.TEST, #"status",
                                                      sampleIDColumnName = METADATA.COL.ID, #"ID",
                                                      sampleExtension = SAMPLE.ID.SUFFIX, #".ULRepGenomesv11.unique.sorted.bam",
                                                      metadataFile = METADATA.PATH,
                                                      rmdDir = rmdDir))
    
    printBpError(tmp)
  }
  summariseMetadataAssocResultsForAll(OUT.DIR)
}else{
  warning(paste0("Not running phenotype/metadata association analysis. Required file does not exist: ",METADATA.PATH))
}

# Test for gene correlations ##########
if(file.exists(KEGG.PATH)){
  print(paste("Testing for gene correlations for",length(allSubstrucSpecies),"species using",getDoParWorkers(),"cores"))

  #tmp <- foreach(spec=allSubstrucSpecies) %dopar% correlateSubpopProfileWithGeneProfiles(spec,OUT.DIR,KEGG.PATH,geneFamilyType="Kegg", corrMethod="pearson")
  tmp <- BiocParallel::bptry(BiocParallel::bplapply(allSubstrucSpecies, BPPARAM = bpParam,
                                                    correlateSubpopProfileWithGeneProfiles,
                                                    OUT.DIR,KEGG.PATH,
                                                    geneFamilyType="Kegg",
                                                    corrMethod="pearson"))
  
  printBpError(tmp)
  #tmp <- foreach(spec=allSubstrucSpecies) %dopar% correlateSubpopProfileWithGeneProfiles(spec,OUT.DIR,KEGG.PATH,geneFamilyType="Kegg", corrMethod="spearman")
  tmp <- BiocParallel::bptry(BiocParallel::bplapply(allSubstrucSpecies, BPPARAM = bpParam,
                                                    correlateSubpopProfileWithGeneProfiles,
                                                    OUT.DIR,KEGG.PATH,
                                                    geneFamilyType="Kegg",
                                                    corrMethod="spearman"))
  
  printBpError(tmp)
  
  if(makeReports){
    # render report
    # tmp <- foreach(spec=allSubstrucSpecies) %dopar%
    #   renderGeneContentReport(speciesID = spec,
    #                                    subpopOutDir = OUT.DIR,
    #                                    geneFamilyAbundancesFile = KEGG.PATH,
    #                                   rmdDir = rmdDir)
    tmp <- BiocParallel::bptry(BiocParallel::bplapply(allSubstrucSpecies, 
                                                      BPPARAM = SerialParam(),#bpParam, # for some reason, parallel fails here
                                                      renderGeneContentReport,
                                                      subpopOutDir = OUT.DIR,
                                                      geneFamilyAbundancesFile = KEGG.PATH,
                                                      rmdDir = rmdDir))
    printBpError(tmp)
  }
  summariseGeneFamilyCorrelationResultsForAll(OUT.DIR)
}else{
  warning(paste0("Not running gene content analysis. Required file does not exist: ",KEGG.PATH))
}

combineAllSummaries(OUT.DIR)
if(makeReports){
  renderResultsSummaryReport(OUT.DIR,rmdDir = rmdDir)
}
#stopCluster(cl)



print(proc.time() - ptm)

