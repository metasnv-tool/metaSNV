
#'@return dataframe with clustering info, or string if clustering failed with description of cause
#'@param dist matrix of distances
#'@param maxPropReadsNonHomog within a sample, what proportion of reads must have the same SNV allele
#' for that position to considered as "homogeneous". Range 0-1. Default = 0.1.  E.g. 0.1 => for a SNV position to be considered
#' in further calculations, > 90% or < 10% of the reads in a sample must have the non-reference allele
#'@param minPropHomogSnvAllelesPerSample within a sample, what proportion of SNV positions need to be "nearly homogeneous" for
#'this sample to be used in cluster medoid determination. Range 0-1. Default = 0.8. E.g. 0.8 => only samples with 80% of their SNVs
#'having nearly the same allele across all samples will be used for medioid cluster determination.
#'@param psCut PS threshold value between 0-1 used to determined the number of clusters
#'@param skipClustAssign if TRUE (default) don't add 'impure' samples to clusters
#'@param removeOutlierSamples if TRUE (default) samples are removed if their mean dissimilarity
#'is more than 5 times the standard deviation more or less than the mean of all samples' mean dissimilarities
#'@param maxNoutliers only used if removeOutlierSamples = TRUE; if there are more than this number of
#'"outliers" then don't remove any. Use 'Inf' to always remove outliers.
#'@param usePackagePredStrength if TRUE uses fpc::prediction.strength to determine
#' best number of clusters. If FALSE uses local implementation used in Costea 2017. Local
#' implementation gives clusters of size = 1 a score of 0 (low),
#' the package implementation gives them a score of 1 (high).
computeClusters <- function(species, dist, doFilterSamplesByAlleleDist,
                            snvFreqs.filtered, filePrefix,
                            outDir, randomSeed,
                            minPropHomogSnvAllelesPerSample = 0.8,
                            maxPropReadsNonHomog=0.1,
                            psCut=0.8,skipClustAssign = TRUE,
                            removeOutlierSamples=TRUE, maxNoutliers = 5,
                            usePackagePredStrength = FALSE){

  numberOfSamplesTotal <- nrow(as.matrix(dist))
  distPostBasicQC <- dist
  # get samples for determining number of clusters and their medoids

  # if there are a small number of samples that are very unusual, remove them
  # here because they can really mess up the clustering metrics & assignments
  x <- removeOutliersFromDistMatrixMinDissim(distMatrix = as.matrix(dist),
                                       maxNoutliers = maxNoutliers,species = filePrefix)
  dist <- x$distMatrix
  outliersRemoved <- x$outliersRemoved

  # remove samples where majority of SNVs are not either very high or very low abundance
  snvFreqPropManyThresh <- computeSnvFreqStats(snvFreqs = snvFreqs.filtered)
  # write this to file after clustering, unless there's an problem (see `if` below)
  homogThreshold <- if_else(maxPropReadsNonHomog > 0.5,
                            1 - maxPropReadsNonHomog, maxPropReadsNonHomog) # e.g. if it's 0.90 -> 0.10
  snvFreqPropThresh <- computeSnvFreqStatsThreshold(snvFreqs = snvFreqs.filtered, homogThreshold = homogThreshold)

  if(doFilterSamplesByAlleleDist & minPropHomogSnvAllelesPerSample > 0){

    distForMedoids <- onlyKeepSamplesWithDistinctFreqs(snvFreqPropByCut = snvFreqPropThresh,
                                                       distMatrix = as.matrix(dist),
                                                       minPropHomogSnvAllelesPerSample = minPropHomogSnvAllelesPerSample)

    if(is.null(nrow(as.matrix(distForMedoids))) || nrow(as.matrix(distForMedoids)) < 6 ){

      warning(paste0(filePrefix,": After removing samples that do not have extreme SNV frequencies, ",
                  "insufficient samples (< 6) remain to pick the number of clusters and cluster medoids. ",
                  "Cutoff for SNV allele frequency: ", maxPropReadsNonHomog," . ",
                  "Cutoff for proportion of SNV positions passing allele frequency cutoff per sample: ",
                  minPropHomogSnvAllelesPerSample,". ",
                  "For values see ",outDir,"/",species,"_freq_composition.tab",
                  ". Aborting for ",filePrefix,"."))

      write.table(snvFreqPropManyThresh,paste(getClustMedoidDefnFailedDir(outDir),species,'_freq_composition.tab',sep=''),sep='\t',quote=F)
      return(paste0("After removing samples that do not have extreme SNV frequencies, ",
             "insufficient samples (< 6) remain to pick the number of clusters and cluster medoids. (n samples = ",
             nrow(as.matrix(distForMedoids)),")" ))
    }
    if( length(unique(unlist(distForMedoids))) <= 1 ){
      msg <- paste0(filePrefix,": After removing samples that do not have extreme SNV frequencies, ",
                     "all values in the distance matrix are equivalent ( all = ",unique(unlist(distForMedoids)),")",
                     ". Aborting for ",filePrefix,".")
      write.table(snvFreqPropManyThresh,paste(getClustMedoidDefnFailedDir(outDir),species,'_freq_composition.tab',sep=''),sep='\t',quote=F)
      warning(msg)
      return(msg)
    }

  }else{
    distForMedoids <- dist
  }

  if(!is.null(nrow(as.matrix(distForMedoids))) && nrow(as.matrix(distForMedoids)) < 100){
    warning(paste0(filePrefix,": fewer than 100 samples being used in determining number of clusters. Results may not be robust. ",
                   "Number of samples used: ",nrow(as.matrix(distForMedoids))))
  }

  # identify number of clusters and their medoids and assign samples to clusters
  clustering <- getClusteringResult(distForMedoids, filePrefix, outDir, randomSeed,
                                    psCut=psCut,
                                    usePackagePredStrength = usePackagePredStrength)
  clustering$numberOfSamplesTotal <- numberOfSamplesTotal
  clustering$distMatrixPostBasicQC <- as.matrix(distPostBasicQC)
  clustering$outliersRemoved <- outliersRemoved
  clustering$species <- species
  if(is.null(clustering)){
    # this often happends when there are a small number of samples in the dist matrix (e.g. <= 6)
    outDir <- getClustMedoidDefnFailedDir(outDir)
  }

  # write results to "noSubstructure" directory
  if( max(clustering$clustering) == 1){
    outDir <- getNoClusteringDir(outDir)
  }

  saveRDS(object = clustering,file = paste(outDir,filePrefix,'_clusteringResult.rds',sep=''))

  # write this now that we know which outDir it should go to
  write.table(snvFreqPropManyThresh,paste(outDir,species,'_freq_composition.tab',sep=''),sep='\t',quote=F)

  # assign any remaining samples to clusters identified from sample subset
  # I'm skipping this (skipClustAssign = TRUE) because I don't want to dilute
  # the clusters before determining distinctive SNVs
  if(skipClustAssign){
    clustDf <- getClustDf(clustering = clustering, dist = dist) #distForMedoids
  }else{
    clustDf <- assignSamplesToClusters(dist = dist, clustering = clustering,
                                       outDir = outDir)
  }
  write.table(clustDf,paste(outDir,filePrefix,'_clustering.tab',sep=''),sep='\t',quote=F)

  # pcoa
  pcoa <- computePCoA(dist,clustDf = clustDf,filePrefix = filePrefix,snvFreqPropThresh, outDir)
  if(!is.null(pcoa)){
    plotPCoAs(pcoa,clustDf,filePrefix,snvFreqPropThresh, outDir, minPropHomogSnvAllelesPerSample)
  }else{
    return(paste0("Error during pcoa - check distance matrix (after optional filtering): ",species))
  }

  if(!is.na(clustering$failureReason)){
    return(clustering$failureReason)
  }

  return(clustDf)
}


#' @return dist distance object with only samples that have distinct frequencies of SNVs
#' @param distMatrix distance matrix can that be coerced to dist object
#' @param minPropHomogSnvAllelesPerSample the proportion of SNVs within a sample that must have a mostly homogeneous allele
onlyKeepSamplesWithDistinctFreqs <- function(snvFreqPropByCut, distMatrix, minPropHomogSnvAllelesPerSample = 0.8){
  #Remove samples that are combinations of multiple alleles per snv position
  # Samples that didn't have vertical coverage at the SNV position (-1) are ignored in 10/90 calc.
  samplesWithEnoughHomogSNVs <- names(snvFreqPropByCut[snvFreqPropByCut >= minPropHomogSnvAllelesPerSample])
  samplesWithEnoughHomogSNVs <- samplesWithEnoughHomogSNVs[!is.na(samplesWithEnoughHomogSNVs)] # can get NAs if there's a NaN or NA in snvFreqs
  samplesWithEnoughHomogSNVs <- samplesWithEnoughHomogSNVs[samplesWithEnoughHomogSNVs %in% rownames(distMatrix)]
  mann_dist <- as.dist(distMatrix[samplesWithEnoughHomogSNVs,samplesWithEnoughHomogSNVs])
  return(mann_dist)
}

#Get prediction strength. This is modified from fpc::prediction.strength
# to work directly on the distance matrix
#'@return NULL if too few samples to cluster from
predStrengthCustom <- function (distance, Gmin = 2, Gmax = 10, M = 50,
                                classification = "centroid", cutoff = 0.8, nnk = 1,
                                clustermethod = cluster::pam, ...)
{
  dist <- as.matrix(distance)
  n <- nrow(dist)
  nf <- c(floor(n*0.5), n - floor(n*0.5))
  indvec <- clcenters <- clusterings <- jclusterings <- classifications <- list()
  prederr <- list()

  for (k in Gmin:Gmax) {
    prederr[[k]] <- numeric(0)
    for (l in 1:M) {
      nperm <- sample(n, n)
      indvec[[l]] <- list()
      indvec[[l]][[1]] <- nperm[1:nf[1]]
      indvec[[l]][[2]] <- nperm[(nf[1] + 1):n]
      for (i in 1:2) {
        if(identical(clustermethod, cluster::pam) ){
          clusterings[[i]] <- as.vector(cluster::pam(as.dist(dist[indvec[[l]][[i]],indvec[[l]][[i]]]), k, diss=TRUE))
        }else{
          clusterings[[i]] <- as.vector(clustermethod(as.dist(dist[indvec[[l]][[i]],indvec[[l]][[i]]]), k, diss=TRUE, ...))$result
        }

        jclusterings[[i]] <- rep(-1, n)
        jclusterings[[i]][indvec[[l]][[i]]] <- clusterings[[i]]$clustering
        centroids <- clusterings[[i]]$medoids
        j <- 3 - i
        classifications[[j]] <- fpc::classifdist(distance, jclusterings[[i]],
                                                 method = classification, centroids = centroids,
                                                 nnk = nnk)[indvec[[l]][[j]]]
      }

      ps_f <- matrix(0, nrow = 2, ncol = k)

      for (i in 1:2) {
        for (kk in 1:k) {
          nik <- sum(clusterings[[i]]$clustering == kk)
          if (nik > 1) {
            a <- which(clusterings[[i]]$clustering[1:(nf[i] - 1)] == kk)
            ps_f[i,kk] <- sum(outer(classifications[[i]][a],classifications[[i]][a],'=='))-length(a)
            ps_f[i,kk] <- ps_f[i, kk]/(nik * (nik - 1))
          }
        }
      }

      ps <- ps_f

      prederr[[k]][l] <- mean(c(min(ps[1, ]), min(ps[2,
                                                     ])))
    }
  }
  mean.pred <- numeric(0)
  if (Gmin > 1)
    mean.pred <- c(1)
  if (Gmin > 2)
    mean.pred <- c(mean.pred, rep(NA, Gmin - 2))
  for (k in Gmin:Gmax) mean.pred <- c(mean.pred, mean(prederr[[k]]))
  optimalk <- max(which(mean.pred > cutoff))
  out <- list(predcorr = prederr, mean.pred = mean.pred, optimalk = optimalk,
              cutoff = cutoff, method = clusterings[[1]]$clustermethod,
              Gmax = Gmax, M = M)
  class(out) <- "predstr"
  out
}

getMaxNumClustersToTry <- function(distObj,warn=T,filePrefix = "",
                                   defaultMaxNumClusters=10,
                                   minClusterSize = 5) {
  Gmax <- defaultMaxNumClusters
  n <- nrow(as.matrix(distObj))
  nf <- c(floor(n*0.5), n - floor(n*0.5))
  maxNclus <- floor(n/minClusterSize)

  if(min(min(nf)-1,maxNclus) < Gmax){
    Gmax <- min( min(nf)-1 , maxNclus)
    if(warn){
      warning(paste0("After filtering, species ",filePrefix," only has ",n,
                     " samples to identify clusters from.",
                     " Only testing prediction strength up to ",Gmax,
                     " clusters instead of up to ",defaultMaxNumClusters))
    }
  }
  return(Gmax)
}

#' @param usePackagePredStrength if TRUE uses fpc::prediction.strength to determine
#' best number of clusters. If FALSE uses local implementation used in Costea 2017. Local
#' implementation gives clusters of size = 1 a score of 0 (low),
#' the package implementation gives them a score of 1 (high).
getClusPredStrengthResult <- function(distObj,warn=T,defaultMaxNumClusters=10,
                                      psCut=0.8,minClusterSize = 5,filePrefix = "",
                                      usePackagePredStrength = FALSE){
  Gmax <- getMaxNumClustersToTry(defaultMaxNumClusters = defaultMaxNumClusters,
                                 distObj = distObj,filePrefix = filePrefix,
                                 minClusterSize = minClusterSize, warn = warn)
  if(Gmax <= 1){
    res <- NULL
  }else{
    #resC <- predStrengthCustom(distDistinct,species=filePrefix,cutoff = psCut)
    # custom version ^ was used originally because the method didn't support
    # a dist matrix as input, now it does. The new method gives different results, especially
    # when the dataset is small. This isn't due to the difference in clustering method implementation
    # (cluster::pam -> fpc::claraCBI).
    # the results look very wrong from the package version when there is a small number of clusters
    # e.g. 10 clusters for 22 samples
    # this is because a cluster of size 1 has a good score in the new implementation and a bad score in the local
    # implementation
    # the implementation in the original paper is better represented by the new package version
    # but the package implementation is also unstable, even when using 100% of the input
    # but perhaps this reflects real uncertainty

    if(usePackagePredStrength){
      res <- fpc::prediction.strength(xdata = distObj,distances = T,cutoff = psCut,
                                      Gmax = Gmax,clustermethod = claraCBI, #pamLike=T, error
                                      Gmin = 2, M = 50,classification = "centroid")  # defaults
    }else{
      res <- predStrengthCustom(distObj,Gmax = Gmax,Gmin = 2,cutoff = psCut,
                                clustermethod = cluster::pam)
    }

  }
  return(res)
}

#'@param dist distance matrix of all samples
#'@param snvFreqPropThresh
#'@param filePrefix string to use as begining of output files e.g. 12345_mann
#'@param randomSeed int seed for determination of number of clusters
#'@param psCut numeric value between 0-1 that is used as cut off for PS
#'value to determine the number of clusters
#'@param usePackagePredStrength if TRUE uses fpc::prediction.strength to determine
#' best number of clusters. If FALSE uses local implementation used in Costea 2017. Local
#' implementation gives clusters of size = 1 a score of 0 (low),
#' the package implementation gives them a score of 1 (high).
#'@return an object of class "pam" representing the clustering. Created by cluster::pam().
#'or 1 if no clusters. null if there's an error
getClusteringResult <- function(distDistinct, filePrefix, outDir, randomSeed,
                                psCut = 0.8, minClusterSize = 3,
                                usePackagePredStrength = FALSE){

  #Get predictions strength for different number of clusters
  #set.seed(randomSeed)

  res <- getClusPredStrengthResult(distDistinct, warn = TRUE, psCut = psCut,
                                   minClusterSize = minClusterSize,
                                   filePrefix = filePrefix,
                                   usePackagePredStrength = usePackagePredStrength,
                                   defaultMaxNumClusters = 15)

  # we used a reduced data set to get the optimal number of clusters
  # next we use that info to actually cluster the data

  if(is.null(res)){
    # this often happends when there are a small number of samples in the dist matrix (e.g. <= 6)
    warning(paste0("Cluster medoid definition failed for: ",filePrefix))

    outDirFailed <- getClustMedoidDefnFailedDir(outDir)
    write.table(as.matrix(distDistinct),
                paste(outDirFailed,filePrefix,'_distMatrixUsedForClustMedoidDefns.txt',sep=''),sep='\t',quote=F)
    #return(NULL)
    numClusters <- 1
    failureReason <- "Cluster medoid definition failed"
  }else{
    numClusters <- res[["optimalk"]]
    failureReason <- NA
  }

  # used to be: clustering <- cluster::pam(distDistinct, numClusters, diss=TRUE)
  # using claraCBI to match what's used in the fpc::prediction.strength
  if(usePackagePredStrength){
    clustering <- fpc::claraCBI(data = distDistinct,k = numClusters,
                                diss = TRUE,usepam = TRUE)$result
  }else{
    clustering <- cluster::pam(distDistinct, numClusters, diss=TRUE)
  }


  clustering$numberOfSamplesUsedForClusterDetection <- nrow(as.matrix(distDistinct))
  clustering$dist <- distDistinct
  clustering$psVals <- res[["mean.pred"]]
  clustering$psValsAll <- res[["predcorr"]]
  clustering$notes <- c()
  clustering$numClusters <- numClusters
  # # doesn't work nicely -- often get weird results -- maybe it would work better if it was directly applied to data, not to dist matrix
  # if(numClusters > 1){
  #   # might want to use this later
  #   clustering$fuzzyClusterResult <- cluster::fanny(x = distDistinct,k = numClusters,diss = T)
  #   # second closest number of clusters
  #   # nextPS <- sort(x,decreasing = T)[2]
  #   # if(nextPS > psCut){
  #   #   nextBestNclust <- which(x == nextPS)
  #   #   clustering$nextBestPs <- list()
  #   #   clustering$nextBestPs$psVal <- nextPS
  #   #   clustering$nextBestPs$fuzzyClusterResult <- cluster::fanny(x = distDistinct,k = nextBestNclust,diss = T)
  #   # }
  # }
  # plot(fuzzyClusterResult,which=1,main = "Fuzzy clustering result")
  # toPlot <- fuzzyClusterResult$membership[names(sort(fuzzyClusterResult$clustering)),]
  # heatmap(toPlot,
  #         scale = "none",Rowv = NA,
  #         #Rowv = order(fuzzyClusterResult$clustering),
  #         RowSideColors = c("red","blue","green")[fuzzyClusterResult$clustering[rownames(toPlot)]])

  # Assess stability of clustering result
  nSamples <- length(labels(distDistinct))

  if(nSamples >= 10){ # need at least 10 samples to do this
    minSamplesToUse <- 10
    lowProp <- max(0.3, ceiling(10/nSamples*10)/10)
    # assess clustering stability
    subsampleProportions<-as.list(seq(from=lowProp,to=1,by=0.1))
    clusNumStabilityIter <- 10 # increases time a lot so 10 is enough
    # assess clustering stability for the number of clusters
    clusNumStability <- getClusNumStability(subsampleProportions = subsampleProportions,
                                            nIterClusStability = clusNumStabilityIter,
                                          distObj=distDistinct, psCut = psCut)
    clusNumStabilityPlots <- getClusNumStabilityPlots(clusNumStability)
    # and assess clustering stability for the cluster membership
    clusMembStability <- getClusMembStability(subsampleProportions = subsampleProportions,
                                              numClusters = numClusters,
                                              distObj=distDistinct)
    clusMembStabilityPlots <- getClusMembStabPlots(clusMembStability)
    clusteringStabilityAssessment <- summariseClusteringStability(nClusStability = clusNumStability,
                                                                  clusMembStability = clusMembStability,
                                                                  numClusters = numClusters)
    clusteringStabilityAssessment$species <- filePrefix

    clustering$stabilityAssessment <- clusteringStabilityAssessment
  }

  # don't count clusters that are too small, they're not useful for us
  clusterSizes <- table(clustering$clustering)
  if(min(table(clustering$clustering)) < minClusterSize){
    # remove clusters smaller than the minimum cluster size, we don't want to use them
    # going forward
    bigEnough <- names(clusterSizes[clusterSizes >= minClusterSize])
    tooSmall <- names(clusterSizes)[!names(clusterSizes) %in% bigEnough]
    clustering$notes <- c(clustering$notes,
                          paste("ids of clusters that are too small:",paste(tooSmall,collapse = ",")))
    clustering$clustering <- clustering$clustering[!clustering$clustering%in%tooSmall]
    clustering$medoids <- clustering$medoids[as.numeric(bigEnough)]
    clustering$id.med <- clustering$id.med[as.numeric(bigEnough)]
    clustering$objective <- NULL
    clustering$isolation <- NULL
    warning(paste0("Species: ",filePrefix,": removed ",length(tooSmall),
                   " clusters because they were too small (",
                  "<",minClusterSize," samples). Num clusters remaining: ",length(bigEnough)))
  }

  numClusters <- length(unique((clustering$clustering)))
  # if 1 cluster, write to a different directory and set value in clustering object
  if(numClusters <= 1){
    #warning(paste0("No substructure for: ",filePrefix," dist. Optimal number of clusters is 1. See files in: ",outDirNoSubStructure))
    clustering$clustering <- rep(x = 1, times = length(clustering$clustering))
    outDir <- getNoClusteringDir(outDir)
  }
  clustering$failureReason <- failureReason


  # for record keeping, in case filtering was done
  write.table(as.matrix(distDistinct),paste(outDir,"/",filePrefix,'_distMatrixUsedForClustMedoidDefns.txt',sep=''),sep='\t',quote=F)

  png(filename = paste(outDir,"/",filePrefix,'_distMatrixUsedForClustMedoidDefns_heatmap.png',sep=''),
      res = 200,width = 40,height = 40,units = "cm")
    heatmap(as.matrix(distDistinct),distfun = as.dist,scale = "none",margins = c(20,20))
  dev.off()

  write.table(res[["mean.pred"]],paste(outDir,filePrefix,'_PS_values.tab',sep=''),sep='\t',quote=F)

  if(exists("clusNumStabilityPlots") & exists("clusMembStabilityPlots")){
    saveClustStabilityPlots(outDir, filePrefix, clusNumStabilityPlots, clusMembStabilityPlots)
  }

  return(clustering)
}


#'@param clustering cluster object produced by getClusteringResult()
#'@return data.frame with samples as row names and
#'1 column called 'clust' which is the clusterID that the sample belongs to
getClustDf <- function(clustering,dist){
  # need to do it this way or row names won't be present
  if (length(unique(clustering$clustering)) == 1) {#There's no clustering here
    df <- data.frame(row.names=colnames(dist),clust=rep(1,ncol(dist)))
  }else{
    df <- data.frame(clust=clustering$clustering)
  }
  return(df)
}

#'@param dist distance matrix of all samples
#'@param clustering cluster object produced by getClusteringResult()
#'@param filePrefix string to use as begining of output files e.g. 12345_mann
#'@param outdir
#'@return data.frame with samples as row names and
#'1 column called 'clust' which is the clusterID that the sample belongs to
assignSamplesToClusters <- function(dist, clustering, outDir){

  if (max(clustering$clustering) == 1) {#There's no clustering here
    df <- data.frame(row.names=colnames(dist),clust=rep(1,ncol(dist)))
  }else {
    #get cluster assignments for all the samples that were left out from cluster definition step
    assign <- rep(-1,ncol(dist))
    names(assign) <- colnames(dist)
    assign[names(clustering$clustering)] <- clustering$clustering
    centroids <- clustering$medoids
    # if any samples were not included in the medoid clustering id step then
    # assign them to a cluster now
    if(-1 %in% assign){
      clust <- fpc::classifdist(dist, assign, method = "centroid", centroids = centroids)
    }else{
      clust <- assign
    }

    df <- data.frame(clust)
  }
  return(df)
}

getEig <- function(pca){
  eig <- pca[["values"]][["Eigenvalues"]]
  eig[eig<0] <- 0
  eig <- eig/sum(eig)*100
  return(eig)
}

getPCoADF <- function(pca,clustDf,freq_data_sample){
  pcoa_df <- data.frame(pca[["vectors"]][,1:2])
  pcoa_df$propFreqHomog <- freq_data_sample[rownames(pcoa_df)]
  pcoa_df$clust <- clustDf[rownames(pcoa_df),]
  return(pcoa_df)
}

computePCoA <- function(dist,clustDf,filePrefix,freq_data_sample, outDir){

  #Get PCOA projection

  # try catch HERE
  pca <- tryCatch({
    return(ape::pcoa(dist))
  }, error = function(e) {
    print(paste0("Error during pcoa for ",filePrefix," - check distance matrix (after optional filtering):", e))
    return(NULL)
  })

  if(is.null(pca) || ncol(pca[["vectors"]]) < 2){
    return(NULL)
  }
  pcoa_df <- getPCoADF(pca,clustDf,freq_data_sample)
  write.table(pcoa_df,paste(outDir,filePrefix,'_pcoa_proj.tab',sep=''),sep='\t',quote=F)

  return(pca)
}

plotPCoAsBoth <- function(pca_m,clustDf_m,pca_a,clustDf_a,freq_data_sample,species, outDir, minPropHomogSnvAllelesPerSample = 0.8){

  eig_m <- getEig(pca_m)
  pcoa_df_mann <- getPCoADF(pca_m,clustDf_m,freq_data_sample)

  eig_a <- getEig(pca_a)
  pcoa_df_allele <- getPCoADF(pca_a,clustDf_a,freq_data_sample)

  png(paste(outDir,species,'_clustering.png',sep=''),res = 200,width = 12,height = 18,units = "cm")
  gridExtra::grid.arrange(
               ggplot2::ggplot(pcoa_df_mann,ggplot2::aes(x=Axis.1,y=Axis.2)) + ggplot2::geom_point(ggplot2::aes(color=propFreqHomog)) +
                 ggplot2::ggtitle(species) + ggplot2::scale_color_gradient2(midpoint=minPropHomogSnvAllelesPerSample) +
                 ggplot2::xlab(sprintf("PC1: %3.2f%%",eig_m[1])) + ggplot2::ylab(sprintf("PC2: %3.2f%%",eig_m[2])),
               ggplot2::ggplot(pcoa_df_mann,ggplot2::aes(x=Axis.1,y=Axis.2)) + ggplot2::geom_point(ggplot2::aes(color=as.factor(clust))) +
                 ggplot2::ggtitle(species),
               ggplot2::ggplot(pcoa_df_allele,ggplot2::aes(x=Axis.1,y=Axis.2)) + ggplot2::geom_point(ggplot2::aes(color=propFreqHomog)) +
                 ggplot2::ggtitle(species) + ggplot2::scale_color_gradient2(midpoint=minPropHomogSnvAllelesPerSample) +
                 ggplot2::xlab(sprintf("allele PC1: %3.2f%%",eig_a[1])) + ggplot2::ylab(sprintf("allele PC2: %3.2f%%",eig_a[2])),
               ggplot2::ggplot(pcoa_df_allele,ggplot2::aes(x=Axis.1,y=Axis.2)) + ggplot2::geom_point(ggplot2::aes(color=as.factor(clust))) +
                 ggplot2::ggtitle(species))
  dev.off()
}

plotPCoAs <- function(pcoa,clustDf,filePrefix,freq_data_sample,outDir, minPropHomogSnvAllelesPerSample = 0.8){

  eig <- getEig(pcoa) # only for axis % var explained
  pcoa_df <- getPCoADF(pcoa,clustDf,freq_data_sample)
  maxFrq <- round(max(pcoa_df$propFreqHomog,na.rm = T),2)
  minFrq <- round(min(pcoa_df$propFreqHomog,na.rm = T),2)
  nClusters <- length(levels(as.factor(pcoa_df$clust)))
  png(paste0(outDir,filePrefix,'_clustering.png'),res = 300,width = 12,height = 18,units = "cm")
  speciesName <- sub(filePrefix,pattern = "_mann",replacement = "")
  gridExtra::grid.arrange(

    ggplot2::ggplot(pcoa_df,ggplot2::aes(x=Axis.1,y=Axis.2)) +
      ggplot2::geom_point(ggplot2::aes(fill=as.factor(clust)),alpha=0.5,
                          colour="black", # for point outlines
                          shape = 21) + # shapes >= 21 have black outlines
      ggplot2::ggtitle(speciesName)+
      scale_fill_brewer(name="Subspecies",palette = "Dark2")+#palette = 3,type = "qual")+
      #scale_fill_viridis_d("Subspecies dominant",alpha=0.7)+
      ggplot2::xlab(sprintf("PC1: %3.1f%%",eig[1])) +
      ggplot2::ylab(sprintf("PC2: %3.1f%%",eig[2]))+
      theme_bw()+theme(legend.position = "bottom",legend.box = "vertical"),

    ggplot(pcoa_df,aes(x=Axis.1,y=Axis.2)) +
      geom_point(aes(fill=round(propFreqHomog,2),
                     shape = as.factor(clust)),size=2,alpha=0.7, colour="black") +
      #ggplot2::ggtitle(paste0(speciesName, " (range:",minFrq,"-",maxFrq ,")")) +
      ggplot2::ggtitle(speciesName)+
      scale_shape_manual(name="Subspecies dominant",
                         values = rep(c(22,23,24,25),times=ceiling(nClusters/4)),
                         na.translate=T,na.value=21) +
      # ggplot2::scale_fill_gradient2(midpoint=minPropHomogSnvAllelesPerSample,
      #                                limits=c(0,1),
      #                                low="firebrick4",mid="lightgoldenrod1",high = "deepskyblue4",
      #                                guide = guide_colourbar(barwidth = 10),
      #                                name="Proportion of SNVs \nwhere allele frequency \nis <10% or >90%") +
      #ggplot2::scale_fill_viridis_c(direction = -1,option = "plasma",alpha=0.7,
      #                                 "Proportion of SNVs \nwhere allele frequency \nis <10% or >90%") +
      #                              guide = guide_colourbar(barwidth = 10,
      ##                                                      #label.theme = element_text(size=)
      ##                                                      ))+
      ggplot2::scale_fill_gradientn(values = c(0,minPropHomogSnvAllelesPerSample,1),
                                    limits=c(0,1),
                                    colors=c("firebrick4","lightgoldenrod1", "darkblue"),
                                    guide = guide_colourbar(barwidth = 10),
                                    name="Proportion of SNV loci \nwhere allele is approx. \n\"homogeneous\" (fixed) ") +
      ggplot2::xlab(sprintf("PC1: %3.1f%%",eig[1])) +
      ggplot2::ylab(sprintf("PC2: %3.1f%%",eig[2]))+
      theme_bw()+theme(legend.position = "bottom",legend.box = "vertical")

  )
  dev.off()
}
