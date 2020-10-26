
#'
#'@param sampleNames character vector of all samples names used in metaSNV - order is important
#'@return data frame of clusters and abundance/frequency within samples or NULL if no genotyping positions were identified
#'@param minGenotypeAbundance minimum mean abundance of genotyping SNVs for sample to be assigned to cluster
#'e.g. if sample X has -1 for 30% of the genotyping SNVs, then it is not assigned to a cluster
#'and not included in results
writeSubpopsForAllSamples <- function(species,sampleNames, outDir,
                                      maxPropUncalledSNV = 0.2,
                                      minGenotypeAbundance = 80){

  # use paul_getPlacingRelevantSubset.py to create 537011_2.pos file from
  # 1. [species]_[cluster]_hap_positions.tab (from writeGenotypeFreqs(...) )
  # 2. SNPs_best_split_[X] from metaSNV output (snpCaller/called_SNPs.best_split_[X])
  # then use paul_convertSNPtoAllelTable.py to creatse 537011_2.pos.freq file
  all_hap <- list.files(path=outDir,paste(species,'_.\\.pos\\.freq$',sep=''),full.names = T)
  if(length(all_hap)==0){
    warning(paste0("Can't find ",species,".*_pos.freq files. ",
                   "Did you run pyGetPlacingRelevantSubset(...) &",
                   " pyConvertSNPtoAllelTable(...) & ",
                   "useGenotypesToProfileSubpops(spec, ",
                   "metaSNVdir=METASNV.DIR, outDir=OUT.DIR )"))
  }

  all_freq <- NULL

  # for each cluster .pos file for this species
  for (d in all_hap) {
    fullData <- read.table(d,header=F,row.names=1,sep='\t')
    colnames(fullData) <- sampleNames # TODO fix this to be less fragile
    fullData <- makeNA(fullData) # change -1 from metaSNV to NA

    # get e.g. 537011_2
    spec_hap <- strsplit(basename(d) ,'\\.')[[1]][1]
    # get the cluster number e.g. ...._2
    uscoreSplit <- strsplit(spec_hap,'_')[[1]]
    cluster <- as.numeric(uscoreSplit[length(uscoreSplit)])
    if(!is.numeric(cluster)){
      stop(paste("Species-cluster name did not conform to expectations. ",
                 "Expected something like: 537011_1 or ref_mOTU_v2_2227_1 where final number after underscore is the cluster.",
                 "Instead got: ", spec_hap))
    }

    hap_info <- read.table(paste(outDir,spec_hap,'_hap_positions.tab',sep=''),
                           header=T,row.names=1,sep='\t',as.is = T)

    # subset full file of SNPs to just get those SNPs specific to this cluster (as opposed to all clusters for this specI)
    # metaSNV ids also change if a annotation file was used or not
    # e.g.
    # 515620.PRJNA29073.CP001104:-:1219674:A  == {species}:{annotation}:{position}:{allele}
    # 515620.PRJNA29073.CP001104:-:1219674:C
    ids_snpCaller_Processed <- sapply(X = strsplit(x = rownames(fullData), split = ":"),
                                      function(id){paste(id[c(1,3,4)],collapse=":")})
    rownames(fullData) <- ids_snpCaller_Processed
    # also change "A>T" --> "T"
    # e.g.
    # 515620.PRJNA29073.CP001104:-:1219674:T>A:.
    # 515620.PRJNA29073.CP001104:-:1219674:T>C:.
    ids_filtered_Processed <- sapply(X=strsplit(x = as.character(hap_info$posId), split = ":"),
                                     function(id){paste(c(id[c(1,3)],
                                                          sub(pattern = ".>",replacement = "",x = id[4])),
                                                        collapse=":")})
    hap_info$posId <- ids_filtered_Processed

    data <- fullData[hap_info$posId,]

    #Flip relevant positions
    if (sum(hap_info[["flip"]]) > 0) {
      flipSNVs <- hap_info[hap_info$flip==TRUE,"posId"]
      data[flipSNVs,] <- 100-data[flipSNVs,]
    }


    # get the columns (samples) where less than 20% of values are NA
    data <- data[,which(apply(data,2,function(x) {sum(!is.na(x))}) >= maxPropUncalledSNV * nrow(data))]

    if( is.null(dim(data)) ){
      #warning(paste("Insufficient data for",species, "cluster",cluster))
      next()
    } else if(is.null(data) || ncol(data) == 0 || nrow(data) == 0){
      #warning(paste("No data for",species, "cluster",cluster))
      next()
    }

    #Compute the frequency for this genotype
    # old way was just to take the median of the SNV frequencies
    #hap_freq <- data.frame(apply(data,2,median,na.rm=T))
    hap_freq <- data.frame(apply(data,2,mean,na.rm=T))

    # new way -- take the percent of genotyping SNVs that are present
    #minAlleleAbundance <- 80 # this percentage of reads must have the allele for it to be considered present
    # hap_freq <- data.frame(apply(data,2,function(x){
    #   genotypePresent <- x > minAlleleAbundance
    #   genotypePresent <- genotypePresent[!is.na(genotypePresent)] # remove NAs
    #   # these are the SNV positions that were not covered with sufficient depth in the sample
    #   prevalenceOfGenotype <- sum(genotypePresent,na.rm = T)/length(genotypePresent)
    #   return(prevalenceOfGenotype*100)
    # }))
    colnames(hap_freq) <- 'freq'
    hap_freq$Cluster <- cluster
    hap_freq$Sample <- rownames(hap_freq)

    all_freq <- rbind(all_freq,hap_freq)
  }

  if(is.null(all_freq)){
    write(file=paste(outDir,species,'_extended_clustering_stat.txt',sep=''),
          x = paste0("Species ",species, ": 0/",length(all_hap)," clusters had usable placing data."),append = T)
    return(NULL)
  }

  #Get samples that could be quantified in all genotypes
  all_freq <- subset(all_freq,!is.na(all_freq$freq))
  t <- table(all_freq$Sample)
  n <- names(which(t == max(t)))
  all_freq <- subset(all_freq, Sample %in% n)

  #This is a stupid way of doing this transformation but i cannot be bothered to do it better
  full <- data.frame(row.names=rownames(subset(all_freq,Cluster==1)))
  for (c in unique(all_freq$Cluster)) {
    s <- subset(all_freq,c==Cluster)
    s <- s[,1,drop=F]
    colnames(s) <- c
    full <- cbind(full,s)
  }

  #not used
  #dd <- read.table(paste(outDir,species,'_hap_freq_median.tab',sep=''),sep='\t',row.names=NULL,header=T)
  #dd <- subset(dd,i==1)

  #Remove samples that don't seem to have a coherent assignment. Not much we can do about these.
  # 120 would mean that SNVs from multiple genotypes were observed at high abundance
  nSampNotClearClus <- sum(rowSums(full)<80 | rowSums(full)>120)
  full <- full[which(rowSums(full)>=80 & rowSums(full)<=120),]

  nSampLowPresence <- sum(rowSums(full)<80)
  nSampMultiPresence <- sum(rowSums(full)<=120)
  full <- full[which(rowSums(full)>=80 & rowSums(full)<=120),]


  if(nSampNotClearClus > 0){
    write(file=paste(outDir,species,'_extended_clustering_stat.txt',sep=''),
          x = paste0("Species ",species, ": ",
                     nSampNotClearClus,' out of ', nrow(full)+nSampNotClearClus ,
                     ' samples rejected due to incoherent subpecies assignment. ',
                     "Number of samples where summed abundance of clusters was < 80%: ",
                     nSampLowPresence,
                     ". Number of samples where summed abundance of clusters was > 120%:",
                     nSampMultiPresence),append = T)
  }



  write.table(full,paste(outDir,species,'_extended_clustering_wFreq.tab',sep=''),sep='\t',quote=F)

  cluster <- data.frame(apply(full,1,function(x){
    assignment <- which(x>minGenotypeAbundance)
    if(length(assignment)==1){return(assignment)}
    return(NA_integer_)
  }))
  colnames(cluster) = 'clust'

  write.table(cluster,paste(outDir,species,'_extended_clustering.tab',sep=''),sep='\t',quote=F)

  return(full)

}

