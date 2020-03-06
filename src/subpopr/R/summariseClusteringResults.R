summariseClusteringResults <- function(rdsPath,resultsDir){
  res <- readRDS(rdsPath)
  speciesID <- res$species
  speciesName <- getSpeciesName(speciesID)
  # get number of clusters
  nClusters <- nrow(res$clusinfo)
  # get PS value
  psVal <- res$psVals[nClusters]
  # get stability of number of clusters
  nClusScore <- res$stabilityAssessment$numClusStabScore
  # get stability of each cluster
  clusScores <- res$stabilityAssessment$clusMembStabScores
  # get number of samples in each cluster
  clusSizes <- res$clusinfo[,"size"]
  names(clusSizes) <- names(clusScores)
  # this approach needed so relative path respects subfolders (like "noClustering/")
  reportAbsPath <- paste0(dirname(rdsPath),"/",speciesID,"_detailedSpeciesReport.html")
  reportRelPath <- sub(pattern = resultsDir,replacement = "./",x = reportAbsPath)

  return(list(speciesID = speciesID,
              speciesName = speciesName,
              #numberOfSamples = res$numberOfSamplesTotal,
              numberOfSamplesUsedForClusterDetection = nrow(as.matrix(res$dist)),
              numberOfClusters = nClusters,
              predictionStrengthValue = round(psVal,digits = 4),
              confidenceInNumberOfClusters = nClusScore,
              confidencePerCluster = paste(clusScores,collapse = "-"),
              clusterSizes = paste(clusSizes,collapse = "-"),
              detailedClusteringResultsFile = reportRelPath))
}

getAllClusteringRdsPaths <- function(resultsDir,distMeth="mann"){
  # for each clustering RDS in main result dir and 2 others
  suffix <- paste0('_',distMeth,'_clusteringResult.rds')
  allRdsFiles <- list.files(recursive = T,
                            path = paste0(resultsDir,"/"),
                            pattern = paste0('.*',suffix),
                            full.names = T)
  #allRdsFiles <- as.character(allRdsFiles)
  speciesNames <- as.character(sapply(allRdsFiles, FUN =
                                        function(x){sub(pattern = suffix,
                                                        replacement = "",
                                                        x = basename(x))} ))
  names(allRdsFiles) <- speciesNames
  return(allRdsFiles)
}

summariseClusteringResultsForAll <- function(resultsDir,distMeth="mann"){
  allRdsFiles <- getAllClusteringRdsPaths(resultsDir,distMeth="mann")
  clustSum <- lapply(allRdsFiles, FUN = summariseClusteringResults, resultsDir=resultsDir)
  df <- do.call(rbind.data.frame,clustSum)
  write.csv(x = df,file = paste0(resultsDir,"/summary_clustering.csv"),quote = T)
  #df$speciesID <- row.names(df)
  saveRDS(df, paste0(resultsDir,"/summary_clustering.rds"))

  #df$species <- row.names(df)

  # allRdsFiles <- getAllClusteringRdsPaths(resultsDir,distMeth="mann")
  # library(purrr)
  # library(tidyr)
  # df <- purrr::map_dfr(.x = allRdsFiles,.f = summariseClusteringResults,.id = "species")
  # df <- tidyr::separate(data = df, col = confidencePerCluster,into = paste0("clusterConfidence_",1:3),fill="right")
  # df <- tidyr::separate(data = df, col = clusterSizes,into = paste0("clusterSize_",1:3),fill="right")
  # readr::write_tsv(x = df,path = paste0(resultsDir,"/summary_clustering.tsv"))
}

# return a list of strings with: (1) species ID, (2) whether the extension worked, (3) new cluster sizes
summariseClusteringExtensionResults <- function(resultsDir,speciesID){

  # handle case where there was only 1 cluster
  initClus <- paste0(getNoClusteringDir(resultsDir),speciesID,"_mann_clustering.tab")
  if(file.exists(initClus)){
    subSpecFreqDf <- read.table(initClus,row.names = 1,as.is = T)
    clusNames <- unique(subSpecFreqDf$clust)
    if(length(clusNames) == 1){
      extWorked <- "No clusters"
      sizesStr <- "NA"
      return(list(speciesID=speciesID,
                  ClusterGenotyping=extWorked,
                  GenotypedClusterSizes=sizesStr))
    }
  }

  # handles cases where there was >1 cluster and extension may have failed
  extClus <- paste0(resultsDir,speciesID,"_extended_clustering.tab")
  if(file.exists(extClus)){
    extWorked <- "Succeeded"
    subSpecFreqDf <- read.table(extClus,row.names = 1,as.is = T)
    sizes <- table(subSpecFreqDf$clust)
    sizesStr <- paste(sizes,collapse = "-")
  }else{
    extWorked <- "Failed"
    sizesStr <- "NA"
  }
  return(list(speciesID=speciesID,
              ClusterGenotyping=extWorked,
              GenotypedClusterSizes=sizesStr))
}

summariseClusteringExtensionResultsForAll <- function(resultsDir,distMeth="mann"){
  allSpeciesRDSs <- getAllClusteringRdsPaths(resultsDir,distMeth="mann")
  allSpeciesNames <- unique(sort(names(allSpeciesRDSs)))
  clustExtSum <- lapply(allSpeciesNames, FUN = summariseClusteringExtensionResults,resultsDir=resultsDir)
  df <- do.call(rbind.data.frame,clustExtSum)
  saveRDS(df, paste0(resultsDir,"/summary_clusteringExtension.rds"))
  write.csv(x = df,file = paste0(resultsDir,"/summary_clusteringExtension.csv"),quote = T,row.names = F)
}



# return a list of strings with: (1) species ID, (2) whether the metadata was tested (3) if any sig assocs were found (4) the html report file
summariseMetadataAssocResults <- function(resultsDir,speciesID,sigCutOff=0.05){

  suffixes <- c("_metadataMultivarOddsRatio_extended.csv","_metadataMultivarOddsRatio_notExtended.csv")
  oddsResFiles <- paste0(resultsDir,speciesID,suffixes)
  anySig <- F
  if(!any(sapply(oddsResFiles,file.exists))){
    anySig <- NA
    extWorked <- "no test results"
    reportRelPath <- NA
  }else{
    for(oddsResFile in oddsResFiles){
      if(file.exists(oddsResFile)){
        extWorked <- "tests performed"
        df <- read.csv(oddsResFile)
        # remove row where predictor column == "(Intercept)"
        df <- df[df$predictor != "(Intercept)",]
        anySig <- anySig | any(df$p < sigCutOff)
      }
    }
    reportRelPath <- paste0("./",speciesID,"_testSubspecMultiPhenoAssoc.html")
    #reportAbsPath <- paste0(resultsDir,"/",reportRelPath)
    #reportRelPath <- ifelse(file.exists(reportAbsPath),reportRelPath,paste("not compiled:",reportAbsPath))
  }
  return(list(speciesID=speciesID,
              assocWithMetadataTested=extWorked,
              anySignifAssocWithMetadata=anySig,
              detailedMetadataAssocResultsFile=reportRelPath))
}

summariseMetadataAssocResultsForAll <- function(resultsDir,distMeth="mann"){
  allSpeciesRDSs <- getAllClusteringRdsPaths(resultsDir,distMeth="mann")
  allSpeciesNames <- unique(sort(names(allSpeciesRDSs)))
  mdAssoc <- lapply(allSpeciesNames,
                        FUN = summariseMetadataAssocResults,
                        resultsDir=resultsDir)
  df <- do.call(rbind.data.frame,mdAssoc)
  saveRDS(df, paste0(resultsDir,"/summary_metadataAssoc.rds"))
  write.csv(x = df,file = paste0(resultsDir,"/summary_metadataAssoc.csv"),quote = T,row.names = F)
}


# return a list of strings with: (1) species ID, (2) whether the metadata was tested (3) if any sig assocs were found (4) the html report file
summariseGeneFamilyCorrelationResults <- function(resultsDir,speciesID,sigCutOff=0.05){
  #resultsDir <- "/Volumes/KESU/scb2/bork/rossum/subspecies/testingSubpopr/inSilicoMock/mine/smallerTestSet/smallGenomes/04_subpopr_nonPack/params.hr5.hs80.ps80/defaults/"
  #speciesID <- "refGenome2clus"

  suffixes <- c("_corrKegg-spearman.tsv","_corrKegg-pearson.tsv")
  resFiles <- paste0(resultsDir,speciesID,suffixes)
  anySig <- F
  if(!any(sapply(resFiles,file.exists))){
    anySig <- NA
    corrWorked <- "no test results"
    reportRelPath <- NA

  }else{
    for(resFile in resFiles){
      if(file.exists(resFile)){
        corrWorked <- "tests performed"
        df <- read.table(resFile,header = T)
        anySig <- anySig | any(df$q.valueBH < sigCutOff)
        rm(df)
      }
    }
    reportRelPath <- paste0("./",speciesID,"_geneContentReport.html")
    #reportAbsPath <- paste0(resultsDir,"/",reportRelPath)
    #reportRelPath <- ifelse(file.exists(reportAbsPath),reportRelPath,paste("not compiled:",reportAbsPath))
  }

  return(list(speciesID=speciesID,
              geneFamCorrTested=corrWorked,
              anySignifGeneFamCorrs=anySig,
              detailedGeneFamCorrResultsFile=reportRelPath))
}

summariseGeneFamilyCorrelationResultsForAll <- function(resultsDir,distMeth="mann"){
  allSpeciesRDSs <- getAllClusteringRdsPaths(resultsDir,distMeth="mann")
  allSpeciesNames <- unique(sort(names(allSpeciesRDSs)))
  mdAssoc <- lapply(allSpeciesNames,
                    FUN = summariseGeneFamilyCorrelationResults,
                    resultsDir=resultsDir)
  df <- do.call(rbind.data.frame,mdAssoc)
  saveRDS(df, paste0(resultsDir,"/summary_geneFamilyCorrAssoc.rds"))
  write.csv(x = df,file = paste0(resultsDir,"/summary_geneFamilyCorrAssoc.csv"),quote = T,row.names = F)
}


combineAllSummaries <- function(resultsDir,distMeth="mann"){

  clustSum <- readRDS(paste0(resultsDir,"/summary_clustering.rds"))
  fullSummary <- clustSum

  if(file.exists(paste0(resultsDir,"/summary_clusteringExtension.rds"))){
    clustExtSum <- readRDS(paste0(resultsDir,"/summary_clusteringExtension.rds"))
    fullSummary <- merge.data.frame(fullSummary,clustExtSum,by = "speciesID",sort = T)
  }

  if(file.exists(paste0(resultsDir,"/summary_metadataAssoc.rds"))){
    mdAssocSum <- readRDS(paste0(resultsDir,"/summary_metadataAssoc.rds"))
    fullSummary <- merge.data.frame(fullSummary,mdAssocSum,by = "speciesID",sort = T)
  }

  if(file.exists(paste0(resultsDir,"/summary_geneFamilyCorrAssoc.rds"))){
    geneFamCorrSum <- readRDS(paste0(resultsDir,"/summary_geneFamilyCorrAssoc.rds"))
    fullSummary <- merge.data.frame(fullSummary,geneFamCorrSum,by = "speciesID",sort = T)
  }

  saveRDS(fullSummary, paste0(resultsDir,"/summary_allResults.rds"))
  write.csv(x = fullSummary,file = paste0(resultsDir,"/summary_allResults.csv"),quote = T,row.names = F)

}

createAllSummaries<-function(resultsDir){
  summariseClusteringResultsForAll(resultsDir)
  summariseClusteringExtensionResultsForAll(resultsDir)
  summariseMetadataAssocResultsForAll(resultsDir)
  summariseGeneFamilyCorrelationResultsForAll(resultsDir)
  combineAllSummaries(resultsDir)
}