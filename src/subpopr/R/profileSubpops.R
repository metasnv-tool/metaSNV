
# TO DO
#' 1. where did All_SpecI_clusters come from?
#' 2. where did RefMGv9.taxid.annot.curated come from?


#' input files required:
#'   distanceMatrixFileMann = paste0("./data/",species,"_mann_distance_matrix.tab") #.filtered.mann.dist"
#'   distanceMatrixFileAllele = paste0("./data/",species,"_allele_distance_matrix.tab") #.filtered.allele.dist"
#'   freqCompFile = paste0("./data/",species,".snp.sample_filtered.freq") #.filtered.snvFreqs"
#'   [species].samples for filtering low coverage samples

#' for extension to all samples:
#' use paul_getPlacingRelevantSubset.py to create 537011_2.pos file from
# 1. [species]_[cluster]_hap_positions.tab (from writeGenotypeFreqs(...) )
# 2. SNPs_best_split_[X] from metaSNV output (snpCaller/called_SNPs.best_split_[X])
#
# then use paul_convertSNPtoAllelTable.py with 537011_2.pos to create 537011_2.pos.freq file


#testing
#defineSubpopulations(species = "357276",metaSNVdir = "/Users/rossum/Dropbox/PostDocBork/subspecies/toolDevelopment/subpopr/ES_PC/metaSNV/fr11_v1/minFilter/outputs_minFilter/",
#                      outDir = "/Users/rossum/Dropbox/PostDocBork/subspecies/toolDevelopment/subpopr/ES_PC/subpop/minFilter/fr11_v1/")

# defineSubpopulations(species = "657317",metaSNVdir = "/Users/rossum/Dropbox/PostDocBork/subspecies/toolDevelopment/testingSubpopr/test2_input/metaSNV/",
#                       outDir = "/Users/rossum/Dropbox/PostDocBork/subspecies/toolDevelopment/testingSubpopr/test2c_output/")
#
# defineSubpopulations(species = "101585",
# metaSNVdir = "~/Dropbox/PostDocBork/subspecies/toolDevelopment/testingSubpopr/inSilicoMock/desman/subpopr_results/metaSNV_results/outputs_minFilter/",
# outDir = "/Users/rossum/Dropbox/PostDocBork/subspecies/toolDevelopment/testingSubpopr/inSilicoMock/desman/subpopr_results/test")

# defineSubpopulations(species="1232452",
#                      metaSNVdir = "~/Dropbox/PostDocBork/subspecies/pancreaticCancer/2018-06_data/metaSNV/stool/outputs_default/",
#                      outDir = "~/Dropbox/PostDocBork/subspecies/pancreaticCancer/2018-06_data/subpopr/stoolDefaults/")

 # defineSubpopulations(species = "referenceGenome",
 #                       metaSNVdir = "~/Dropbox/PostDocBork/subspecies/toolDevelopment/testingSubpopr/inSilicoMock/mine/metaSNV/outputs_defaults/",
 #                       outDir = "/Users/rossum/Dropbox/PostDocBork/subspecies/toolDevelopment/testingSubpopr/inSilicoMock/mine/subpopr/resultsLocal/",
 #                      propHomogAllelesPerSampleCutoff = 0)


#' @return int the number of clusters or -1 if there is an insufficient number of samples in the metaSNV results used as input
#' @export
#' @param species int the speciesID
#' @param metaSNVdir location of the metaSNV output files
#' @param outDir location for the output files
#' @param randomSeed int used for cluster generation, default is random number
#' @param doFilterSamplesByAlleleDist when picking number of clusters and cluster medoids, should samples be
#' pre-filtered to only keep ones where most
#' (see minPropHomogSnvAllelesPerSample) SNV positions are dominated by 1 allele
#' @param maxPropReadsNonHomog within a sample, what proportion of reads must have the same SNV allele
#' for that position to considered as "homogeneous". Range 0-1. Default = 0.1.  E.g. 0.1 => for a SNV position to be considered
#' "nearly homogeneous" for further calculations, > 90% or < 10% of the reads in a sample must have the non-reference allele
#' @param minPropHomogSnvAllelesPerSample proportion of SNV positions that need to be dominated by one
#' allele for pre-filtering (see doFilterSamplesByAlleleDist)
#' @param psCut PS threshold value between 0-1 used to determined the number of clusters
#' @param uniqSubpopSnvFreqThreshold value between 0-1 that specifies how much more abundant a SNV must be
#' in samples representing a subspecies vs in other samples
#' @param bamFileNamesToUsePath path to a file with one bam file name per line. Only samples listed in this file will
#' be used for analysis. Note that SNV filtering is not re-done. To use all samples, leave as default value: NULL
defineSubpopulationsMannAndAllele <- function(species, metaSNVdir, outDir, randomSeed=NULL, minNumberOfSamplesToStart = 4,
                                 doFilterSamplesByAlleleDist = T,
                                 maxPropReadsNonHomog = 0.1,
                                 minPropHomogSnvAllelesPerSample = 0.8,
                                 psCut = 0.8,
                                 uniqSubpopSnvFreqThreshold=0.8,
                                 useAlleleDistances=T,
                                 bamFileNamesToUsePath=NULL){

  randomSeed <- if_else(is.null(randomSeed),false = randomSeed,true = sample(x = 4124:4613646,size = 1))
  # import data
  # old verison
  # distanceMatrixFileMann <- paste0(metaSNVdir,"/distances/",species,"_mann_distance_matrix.tab") #.filtered.mann.dist"
  # distanceMatrixFileAllele <- paste0(metaSNVdir,"/distances/",species,"_allele_distance_matrix.tab") #.filtered.allele.dist"
  # freqCompFile <- paste0(metaSNVdir,"/filtered/pop/",species,".snp.sample_filtered.freq") #.filtered.snvFreqs"


  # new metasnv_post version
  distanceMatrixFileMann <- paste0(metaSNVdir,"/distances/",species,".filtered.mann.dist")
  freqCompFile <- paste0(metaSNVdir,"/filtered/pop/",species,".filtered.freq") #.filtered.snvFreqs"

  if(!dir.exists(outDir)){
    dir.create(outDir)
  }

  # check for required files
  if(!dir.exists(metaSNVdir)){
    stop("Missing directory: ", metaSNVdir)
  }
  if( ! file.exists(distanceMatrixFileMann) ){
    stop("Missing file: ", distanceMatrixFileMann)
  }

  if( ! file.exists(freqCompFile) ){
    stop("Missing file: ", freqCompFile)
  }

  distMann <- read.table(distanceMatrixFileMann,header=T,row.names=1,check.names=F,strip.white = F,sep="\t")
  snvFreqs.filtered <- read.table(freqCompFile,header=T,row.names=1,check.names=F,strip.white = F,sep="\t")

  # sometimes there are blanks in the metaSNV distance output due to too many -1s
  # remove those here
  distMann <- rmNAfromDistMatrix(distMann)

  # verify that snvFreq sample names and distance matrix samples names match
  if(length(colnames(snvFreqs.filtered)) != length(colnames(distMann)) ||
     all(colnames(snvFreqs.filtered) != colnames(distMann))){
    samps <- intersect(colnames(distMann),colnames(snvFreqs.filtered))
    snvFreqs.filtered <- snvFreqs.filtered[samps]
    distMann <- distMann[samps,samps]
    warning(paste("Samples names do not match in the distance matrix and the filtered SNV file for species",species,".",
                  "Using the intersection, which is ",length(samps)," samples."))
    if(length(samps) < minNumberOfSamplesToStart){
      stop(paste0("Too few samples remain after selecting only those in the distance and SNP files. At least ",
                  minNumberOfSamplesToStart," are required for analysis. Aborting."))
    }
  }


  if(useAlleleDistances){
    distanceMatrixFileAllele <- paste0(metaSNVdir,"/distances/",species,".filtered.allele.dist")
    if( ! file.exists(distanceMatrixFileAllele) ){
      stop("Missing file: ", distanceMatrixFileAllele)
    }
    distAllele <- read.table(distanceMatrixFileAllele,header=T,row.names=1,check.names=F,strip.white = F,sep="\t")
    distAllele <- rmNAfromDistMatrix(distAllele)
    # verify that snvFreq sample names and distance matrix samples names match
    if(length(colnames(snvFreqs.filtered)) != length(colnames(distAllele)) ||
       all(colnames(snvFreqs.filtered) != colnames(distAllele))){
      samps <- intersect(colnames(distAllele),colnames(snvFreqs.filtered))
      snvFreqs.filtered <- snvFreqs.filtered[samps]
      distAllele <- distAllele[samps,samps]
      warning(paste("Samples names do not match in the distance matrix and the filtered SNV file for species",species,".",
                    "Using the intersection, which is ",length(samps)," samples."))
    }
  }

  if(!is.null(bamFileNamesToUsePath)){
    if( ! file.exists(bamFileNamesToUsePath) ){
      warning(paste0("Samples not subselected according to specified file as file does not exist: ", bamFileNamesToUsePath))
    }else{
      bamFileNamesToUse <- read.delim(file = bamFileNamesToUsePath,
                                      sep = "\n",header = F,as.is = T,
                                      col.names = c("bamNames"))[,"bamNames",T]
      if(length(bamFileNamesToUse) < minNumberOfSamplesToStart){
        stop(paste0("Insufficient samples would remain after selecting samples based on file :",bamFileNamesToUsePath," .",
                    "At least ", minNumberOfSamplesToStart, " samples are required. ",
                    "Only ",length(bamFileNamesToUse)," samples would be selected. Aborting."))
      }
      allSamps <- intersect(colnames(distMann),colnames(snvFreqs.filtered))
      sampsToKeep <- intersect(allSamps,bamFileNamesToUse)
      if(length(sampsToKeep) < minNumberOfSamplesToStart){
        stop(paste0("Insufficient samples remain after selecting samples based on file :",bamFileNamesToUsePath," .",
                    "Only ",length(sampsToKeep)," samples remain. Aborting. \n Maybe format is wrong? Example sample name: ",allSamps[1]))
      }
      warning(paste0("Samples subselected according to specified file: ", bamFileNamesToUsePath,". Samples remaining: ", length(sampsToKeep)))
      snvFreqs.filtered <- snvFreqs.filtered[sampsToKeep]
      distMann <- distMann[sampsToKeep,sampsToKeep]
    }
  }


  # need at least 4 samples for method to work (more required for it to be reasonable)
  if( ncol(snvFreqs.filtered) < minNumberOfSamplesToStart ){
    warning(paste0("Insufficient number of samples in metaSNV filtered SNV results for species: ",species, " (",ncol(distMann)," samples)"))
    return(paste0("Insufficient number of samples in metaSNV filtered SNV results (",ncol(distMann)," samples)"))
  }
  # need at least 4 samples for method to work (more required for it to be reasonable)
  if( is.null(distMann) || length(distMann) < 2 ||  ncol(distMann) < minNumberOfSamplesToStart ){
    warning(paste0("Insufficient number of samples in metaSNV dist matrix results for species: ",species, " (",ncol(distMann)," samples)"))
    return(paste0("Insufficient number of samples in metaSNV dist matrix results (",ncol(distMann)," samples)"))
  }

  # new version of metaSNV has SNV frequencies as [0,1]; old version (which code was initially written with) had them as [0,100]
  # if not -1, multiply by 100 (to do eventually is replace -1 with NA, but may have implications throughout code)
  if(max(snvFreqs.filtered,na.rm = T)>1){
    error("The max value in your SNV frequency file is > 1. Values are expected to be -1 or within [0 to 1]. Are you using an old version of metaSNV? ")
  }
  # change range to 0-100 to avoid changing code everywhere else
  snvFreqs.filtered[snvFreqs.filtered!=-1] <- snvFreqs.filtered[snvFreqs.filtered!=-1]*100

  snvFreqPlot(species,snvFreqs.filtered,
               outDir = getSnvFreqPlotDir(outDir),
               minPropHomogSnvAllelesPerSample,
              maxPropReadsNonHomog = maxPropReadsNonHomog)

  filePrefixMann=paste0(species,"_mann")
  filePrefixAllele=paste0(species,"_allele")

  # identify clusters from subset of data
  clustDf_m <- computeClusters(dist=distMann,
                               species = species,
                               doFilterSamplesByAlleleDist =  doFilterSamplesByAlleleDist,
                               minPropHomogSnvAllelesPerSample = minPropHomogSnvAllelesPerSample,
                               snvFreqs.filtered = snvFreqs.filtered,
                               filePrefix = filePrefixMann,
                               outDir = outDir,
                               #randomSeed = randomSeed,
                               maxPropReadsNonHomog = maxPropReadsNonHomog,
                               psCut = psCut)

  if(useAlleleDistances){
  clustDf_a <- computeClusters(dist=distAllele,
                               species = species,
                               doFilterSamplesByAlleleDist = doFilterSamplesByAlleleDist,
                               minPropHomogSnvAllelesPerSample = minPropHomogSnvAllelesPerSample,
                               snvFreqs.filtered = snvFreqs.filtered,
                               filePrefix = filePrefixAllele,
                               outDir = outDir,
                               #randomSeed = randomSeed,
                               maxPropReadsNonHomog = maxPropReadsNonHomog,
                               psCut = psCut)
  }
  if(is.null(clustDf_m) || is.character(clustDf_m) ){
    return(clustDf_m)
  }

  #if no mann sub populations, then we're done
  if(max(clustDf_m$clust) < 2){
    #print(paste("Species ",species," does not have any subpopulation structure with mann clust."))
    return("nClusters = 1")
  }

  # Compute variances and percentage explained by the clustering
  var_m <- variationExplainedByClusters(clustDf_m, snvFreqs.filtered, species, majorAllele = F,outDir)
  if(useAlleleDistances){
    var_a <- variationExplainedByClusters(clustDf_a, snvFreqs.filtered, species, majorAllele = T,outDir)
  }

  # compute subspecies' distinctive SNVs
  writeGenotypeFreqs(clustDf_m, snvFreqs.filtered, species, outDir,uniqSubpopSnvFreqThreshold)

  return(paste0("nClusters = ",max(clustDf_m$clust)))
}









#' @return int the number of clusters or -1 if there is an insufficient number of samples in the metaSNV results used as input
#' @export
#' @param species int the speciesID
#' @param metaSNVdir location of the metaSNV output files
#' @param outDir location for the output files
#' @param randomSeed int used for cluster generation, default is random number
#' @param minNumberOfSamplesToStart default is 100, anything less makes clustering unstable. Not advised to be below 100.
#' @param doFilterSamplesByAlleleDist when picking number of clusters and cluster medoids, should samples be
#' pre-filtered to only keep ones where most
#' (see minPropHomogSnvAllelesPerSample) SNV positions are dominated by 1 allele
#' @param maxPropReadsNonHomog within a sample, what proportion of reads must have the same SNV allele
#' for that position to considered as "homogeneous". Range 0-1. Default = 0.1.  E.g. 0.1 => for a SNV position to be considered
#' "nearly homogeneous" for further calculations, > 90% or < 10% of the reads in a sample must have the non-reference allele
#' @param minPropHomogSnvAllelesPerSample proportion of SNV positions that need to be dominated by one
#' allele for pre-filtering (see doFilterSamplesByAlleleDist)
#' @param psCut PS threshold value between 0-1 used to determined the number of clusters
#' @param uniqSubpopSnvFreqThreshold value between 0-1 that specifies how much more abundant a SNV must be
#' in samples representing a subspecies vs in other samples
#' @param bamFileNamesToUsePath path to a file with one bam file name per line. Only samples listed in this file will
#' be used for analysis. Note that SNV filtering is not re-done. To use all samples, leave as default value: NULL
#' @param usePackagePredStrength if TRUE uses fpc::prediction.strength to determine
#' best number of clusters. If FALSE uses local implementation used in Costea 2017. Local
#' implementation gives clusters of size = 1 a score of 0 (low),
#' the package implementation gives them a score of 1 (high).
defineSubpopulations <- function(species, distName = "mann",
                                 metaSNVdir, outDir,
                                 randomSeed=NULL,
                                 minNumberOfSamplesToStart = 100, # anything less makes clustering unstable
                                 doFilterSamplesByAlleleDist = T,
                                 maxPropReadsNonHomog = 0.1,
                                 minPropHomogSnvAllelesPerSample = 0.8,
                                 psCut = 0.8,
                                 uniqSubpopSnvFreqThreshold=0.8,
                                 bamFileNamesToUsePath = NULL,
                                 usePackagePredStrength = FALSE){

  #randomSeed <- ifelse(!is.null(randomSeed),yes = randomSeed,no = sample(x = 4124:4613646,size = 1))

  # new metasnv_post version
  distanceMatrixFile <- paste0(metaSNVdir,"/distances/",species,".filtered.",distName,".dist")
  freqCompFile <- paste0(metaSNVdir,"/filtered/pop/",species,".filtered.freq") #.filtered.snvFreqs"

  if(!dir.exists(outDir)){
    dir.create(outDir)
  }

  # check for required files
  if(!dir.exists(metaSNVdir)){
    stop("Missing directory: ", metaSNVdir)
  }
  if( ! file.exists(distanceMatrixFile) ){
    stop("Missing file: ", distanceMatrixFile)
  }

  if( ! file.exists(freqCompFile) ){
    stop("Missing file: ", freqCompFile)
  }

  distMa <- read.table(distanceMatrixFile,header=T,row.names=1,check.names=F,strip.white = F,sep="\t")
  snvFreqs.filtered <- read.table(freqCompFile,header=T,row.names=1,check.names=F,strip.white = F,sep="\t")

  # sometimes there are blanks in the metaSNV distance output due to too many -1s
  # remove those here
  distMa <- rmNAfromDistMatrix(distMa)

  # verify that snvFreq sample names and distance matrix samples names match
  if( length(colnames(snvFreqs.filtered) ) != length(colnames(distMa)) |
      any(colnames(snvFreqs.filtered) != colnames(distMa)) ){
    samps <- intersect(colnames(distMa),colnames(snvFreqs.filtered))
    snvFreqs.filtered <- snvFreqs.filtered[samps]
    distMa <- distMa[samps,samps]
    warning(paste("Samples names do not match in the distance matrix and the filtered SNV file for species",species,".",
                  "Using the intersection, which is ",length(samps)," samples."))
    if(length(samps) < minNumberOfSamplesToStart){
      stop(paste0("Too few samples remain after selecting only those in the distance and SNP files. At least ",
                  minNumberOfSamplesToStart," are required for analysis. Aborting."))
    }
  }

  if(!is.null(bamFileNamesToUsePath)){
    if( ! file.exists(bamFileNamesToUsePath) ){
      warning(paste0("Samples not subselected according to specified file as file does not exist: ", bamFileNamesToUsePath))
    }else{
      bamFileNamesToUse <- read.delim(file = bamFileNamesToUsePath,
                                      sep = "\n",header = F,as.is = T,
                                      col.names = c("bamNames"))[,"bamNames",T]
      if(length(bamFileNamesToUse) < minNumberOfSamplesToStart){
        stop(paste0("Insufficient samples would remain after selecting samples based on file :",bamFileNamesToUsePath," .",
                    "At least ", minNumberOfSamplesToStart, " samples are required. ",
                    "Only ",length(bamFileNamesToUse)," samples would be selected. Aborting."))
      }
      allSamps <- intersect(colnames(distMa),colnames(snvFreqs.filtered))
      sampsToKeep <- intersect(allSamps,bamFileNamesToUse)
      if(length(sampsToKeep) < minNumberOfSamplesToStart){
        stop(paste0("Insufficient samples remain after selecting samples based on file :",bamFileNamesToUsePath," .",
                    "Only ",length(sampsToKeep)," samples remain. Aborting. \n Maybe format is wrong? Example sample name: ",allSamps[1]))
      }
      warning(paste0("Samples subselected according to specified file: ", bamFileNamesToUsePath,". Number of samples remaining: ", length(sampsToKeep)))
      snvFreqs.filtered <- snvFreqs.filtered[,sampsToKeep]
      distMa <- distMa[sampsToKeep,sampsToKeep]
    }
  }


  # need at least 4 samples for method to work (more required for it to be reasonable)
  if( ncol(snvFreqs.filtered) < minNumberOfSamplesToStart ){
    warning(paste0("Insufficient number of samples in metaSNV filtered SNV results for species: ",species, " (",ncol(distMa)," samples)"))
    return(paste0("Insufficient number of samples in metaSNV filtered SNV results (",ncol(distMa)," samples)"))
  }
  # need at least 4 samples for method to work (more required for it to be reasonable)
  if( is.null(distMa) || length(distMa) < 2 ||  ncol(distMa) < minNumberOfSamplesToStart ){
    warning(paste0("Insufficient number of samples in metaSNV dist matrix results for species: ",species, " (",ncol(distMa)," samples)"))
    return(paste0("Insufficient number of samples in metaSNV dist matrix results (",ncol(distMa)," samples)"))
  }

  # new version of metaSNV has SNV frequencies as [0,1]; old version (which code was initially written with) had them as [0,100]
  # if not -1, multiply by 100 (to do eventually is replace -1 with NA, but may have implications throughout code)
  if(max(snvFreqs.filtered,na.rm = T)>1){
    error("The max value in your SNV frequency file is > 1. Values are expected to be -1 or within [0 to 1]. Are you using an old version of metaSNV? ")
  }
  # change range to 0-100 to avoid changing code everywhere else
  snvFreqs.filtered[snvFreqs.filtered!=-1] <- snvFreqs.filtered[snvFreqs.filtered!=-1]*100

  snvFreqPlot(species,snvFreqs.filtered,
              outDir = getSnvFreqPlotDir(outDir),
              minPropHomogSnvAllelesPerSample,
              maxPropReadsNonHomog = maxPropReadsNonHomog)

  filePrefix=paste0(species,"_",distName)

  # identify clusters from subset of data
  clustDf <- computeClusters(dist=distMa,
                               species = species,
                               doFilterSamplesByAlleleDist =  doFilterSamplesByAlleleDist,
                               minPropHomogSnvAllelesPerSample = minPropHomogSnvAllelesPerSample,
                               snvFreqs.filtered = snvFreqs.filtered,
                               filePrefix = filePrefix,
                               outDir = outDir,
                               randomSeed = randomSeed,
                               maxPropReadsNonHomog = maxPropReadsNonHomog,
                               psCut = psCut,
                               usePackagePredStrength = usePackagePredStrength)

  if(is.null(clustDf) || is.character(clustDf) ){
    return(clustDf)
  }

  #if no sub populations, then we're done
  if(length(unique(clustDf$clust)) <= 1){
    return("nClusters = 1")
  }

  # Compute variances and percentage explained by the clustering
  varExp <- variationExplainedByClusters(clustDf, snvFreqs.filtered, species, filePrefix, majorAllele = F,outDir)

  # compute subspecies' distinctive SNVs
  writeGenotypeFreqs(clustDf, snvFreqs.filtered, species, outDir,uniqSubpopSnvFreqThreshold)

  return(paste0("nClusters = ",length(unique(clustDf$clust)) ))
}




# useGenotypesToProfileSubpops(species = "357276",
#                               motuProfileFilePath ="/Users/rossum/Dropbox/PostDocBork/subspecies/pancreaticCancer/motus2/all_samples.motusv2.relabund.renamedSamples.tsv",
#                               outDir = "/Users/rossum/Dropbox/PostDocBork/subspecies/toolDevelopment/subpopr/ES_PC/metaSNV/fr11_v1/minFilter/",
#                               metaSNVdir = "/Users/rossum/Dropbox/PostDocBork/subspecies/toolDevelopment/subpopr/ES_PC/metaSNV/fr11_v1/minFilter/outputs_minFilter")

#' @export
useGenotypesToProfileSubpops <- function(species, metaSNVdir, outDir){

  # use the genotype profiles to classify all samples
  # reads in output from writeGenotypeFreqs

  # use paul_getPlacingRelevantSubset.py to create 537011_2.pos file from
  # 1. [species]_[cluster]_hap_positions.tab (from writeGenotypeFreqs(...) )
  # 2. SNPs_best_split_[X] from metaSNV output (snpCaller/called_SNPs.best_split_[X])
  # then use paul_convertSNPtoAllelTable.py with 537011_2.pos to create 537011_2.pos.freq file
  # see http://blog.quantitations.com/tutorial/2012/11/17/including-python-code-in-an-r-package/
  # TODO add sample names to python script outputs

  # no distinctive genotyping positions for this species, perhaps the cutoff was bad
  if(!file.exists(paste0(outDir,species,'_hap_freq_median.tab'))){
    warning(paste0("No distinctive genotyping positions for species ",species,", perhaps the cutoff was bad. ",
                   "See ",paste0(outDir,species,'_hap_out.txt')," for details. ",
                  "(Required file does not exist: ", paste0(outDir,species,'_hap_freq_median.tab)')))
    return(NULL)
  }

  if(!file.exists(paste0(metaSNVdir,'/all_samples'))){
    stop(paste0("File with paths to samples does not exist. This should have been created by metaSNV.",
                   "Expected path: ",paste0(metaSNVdir,'/all_samples')))
    return(NULL)
  }


  # metaSNV makes a file called all_samples, this is required because the columns in the SNPcall output files
  # don't have headers and the columns are in the order indicated in this file
  ss <- read.table(paste0(metaSNVdir,'/all_samples'),header=F,check.names=F)
  sampleNamesInMetaSNVorder <- as.character(ss$V1)
  sampleNamesInMetaSNVorder <- basename(sampleNamesInMetaSNVorder)

  extendedClusteringFreqs <- writeSubpopsForAllSamples(species, sampleNamesInMetaSNVorder, outDir)

  if(is.null(extendedClusteringFreqs) || nrow(extendedClusteringFreqs) == 0){
    warning(paste0("Species does not have cluster frequencies. Species: ",species))
  }

  return(NULL)
}

#'get subspecies abundance using species abundance and relative abundance of subspecies
useSpeciesAbundToCalcSubspeciesAbund <- function(species, speciesAbundanceProfileFilePath, outDir,speciesProfileIsMotus){

  extendedClusteringFreqs <- read.table(paste(outDir,species,'_extended_clustering_wFreq.tab',sep=''),sep='\t')

  if(is.null(extendedClusteringFreqs) || nrow(extendedClusteringFreqs) == 0){
    warning(paste0("Species does not have cluster frequencies. Species: ",species))
    return(NULL)
  }

  if(!file.exists(speciesAbundanceProfileFilePath)){
    stop(paste0("Cannot compute overall subspecies abundances. ",
                "Required species abundance file not found: ",
                speciesAbundanceProfileFilePathß))
    return(NULL)
  }
  # we have the relative abundance of the species' subpops in each cluster
  # now we scale the abundance of these subpops by the abundance of the species

  if(speciesProfileIsMotus){
    writeSubpopAbundMotusProfile(species, clusterFreqs = extendedClusteringFreqs,
                                 outDir, motuProfileFilePath = speciesAbundanceProfileFilePath)

  }else{
    writeSubpopAbundSpeciesAbund(species, clusterFreqs = extendedClusteringFreqs,
                                 outDir=outDir,speciesProfileFilePath= speciesAbundanceProfileFilePath)
  }
  return(NULL)
}

