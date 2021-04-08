
#'
#'@param sampleNames character vector of all samples names used in metaSNV - order is important
#'@return data frame of clusters and abundance/frequency within samples or NULL if no genotyping positions were identified
#'@param minGenotypeAbundance minimum mean abundance of genotyping SNVs for sample to be assigned to cluster
#'e.g. if sample X has -1 for 30% of the genotyping SNVs, then it is not assigned to a cluster
#'and not included in results
writeSubpopsForAllSamples <- function(species,sampleNames, outDir,
                                      maxPropUncalledSNV = 0.2,
                                      minGenotypeAbundance = 80){

  # use getGenotypingSNVSubset.py to create 537011_2.pos file from
  # 1. [species]_[cluster]_hap_positions.tab (from writeGenotypeFreqs(...) )
  # 2. SNPs_best_split_[X] from metaSNV output (snpCaller/called_SNPs.best_split_[X])
  # then use convertSNVtoAlleleFreq.py to creatse 537011_2.pos.freq file
  all_hap <- list.files(path=outDir,paste(species,'_.*\\.pos\\.freq$',sep=''),full.names = T)
  if(length(all_hap)==0){
    warning(paste0("Can't find ",species,".*_pos.freq files. ",
                   "Did you run pyGetPlacingRelevantSubset(...) &",
                   " pyConvertSNPtoAllelTable(...) & ",
                   "useGenotypesToProfileSubpops(spec, ",
                   "metaSNVdir=METASNV.DIR, outDir=OUT.DIR )"))
  }

  all_freq <- NULL # median
  all_freq_summary <- NULL

  # for each cluster .pos file for this species
  for (d in all_hap) {
    fullData <- read.table(d,header=F,row.names=1,sep='\t')
    if(length(sampleNames) != ncol(fullData) ){
      stop(paste0("File: ",d," does not have expected number of columns (each column is a sample).",
                  " According to all_samples, it should have ",length(sampleNames)," columns."))
    }
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
    if (sum(hap_info[["flip"]]) > 0) { #if any flips
      flipSNVs <- hap_info[hap_info$flip==TRUE,"posId"]
      data[flipSNVs,] <- 100-data[flipSNVs,]
    }

    # get the columns (samples) where less than 20% of values are NA
    data <- data[,which(apply(data,2,function(x) {sum(!is.na(x))}) >= maxPropUncalledSNV * nrow(data))]

    if( is.null(dim(data)) ){
      warning(paste("Insufficient data for",species, "cluster",cluster))
      next()
    } else if(is.null(data) || ncol(data) == 0 || nrow(data) == 0){
      warning(paste("No data for",species, "cluster",cluster))
      next()
    }

    #Compute the frequency for this genotype
    # just use the median of the SNV frequencies (though I print more info here)
    hap_freq_median <- data.frame(apply(data,2,median,na.rm=T))
    hap_freq_summary <- data.frame(Sample = colnames(data),
                                   Cluster = cluster,
                                   mean=apply(data,2,mean,na.rm=T),
                                   median=apply(data,2,median,na.rm=T),
                                   standardDeviation=apply(data,2,sd,na.rm=T),
                                   prevalence=apply(data,2,function(x){x<-x[!is.na(x)]; sum(x>0)/length(x)}),
                                   prevalenceGte5=apply(data,2,function(x){x<-x[!is.na(x)]; sum(x>=5)/length(x)}),
                                   n0 = apply(data,2,function(x){x<-x[!is.na(x)]; sum(x==0)}),
                                   n100 = apply(data,2,function(x){x<-x[!is.na(x)]; sum(x==100)}),
                                   nNoCoverage=apply(data,2,function(x){sum(is.na(x))}))

    colnames(hap_freq_median) <- 'freq' # dataframe with 1 column and samples as row names
    hap_freq_median$Cluster <- cluster
    hap_freq_median$Sample <- rownames(hap_freq_median)
    all_freq <- rbind(all_freq,hap_freq_median)


    all_freq_summary <- rbind(all_freq_summary,hap_freq_summary)

  }

  if(is.null(all_freq)){
    write(file=paste(outDir,species,'_extended_clustering_stat.txt',sep=''),
          x = paste0("Species ",species, ": 0/",length(all_hap)," clusters had usable placing data."),append = T)
    return(NULL)
  }

  write.table(all_freq_summary,row.names = F,
              paste(outDir,species,'_extended_clustering_abundanceSummaryStats.tsv',sep=''),
              sep='\t',quote=F)


  rmNAandSpread <- function(df){
    #Get samples that could be quantified in all genotypes
    df <- subset(df,!is.na(df$freq))
    t <- table(df$Sample)
    n <- names(which(t == max(t)))
    df <- subset(df, Sample %in% n)

    #This is a stupid way of doing this transformation but i cannot be bothered to do it better
    df_wide <- data.frame(row.names=rownames(subset(df,Cluster==1)))
    for (c in unique(df$Cluster)) {
      s <- subset(df,c==Cluster)
      s <- s[,1,drop=F]
      colnames(s) <- c
      df_wide <- cbind(df_wide,s)
    }
    return(df_wide)
  }

  full <- rmNAandSpread(all_freq)

  #not used
  #dd <- read.table(paste(outDir,species,'_hap_freq_median.tab',sep=''),sep='\t',row.names=NULL,header=T)
  #dd <- subset(dd,i==1)

  #Remove samples that don't seem to have a coherent assignment. Not much we can do about these.
  # 120 would mean that SNVs from multiple genotypes were observed at high abundance
  nSampNotClearClus <- sum(rowSums(full)<80 | rowSums(full)>120)
  nSampLowPresence <- sum(rowSums(full)<80)
  nSampMultiPresence <- sum(rowSums(full)>120)

  filtered <- full[which(rowSums(full)>=80 & rowSums(full)<=120),]

  if(nSampNotClearClus > 0){
    write(file=paste(outDir,species,'_extended_clustering_stat.txt',sep=''),
          x = paste0("Species ",species, ": ",
                     nSampNotClearClus,' out of ', nrow(full) ,
                     ' samples rejected due to incoherent subpecies assignment. ',
                     "Number of samples where summed abundance of clusters was < 80%: ",
                     nSampLowPresence,
                     ". Number of samples where summed abundance of clusters was > 120%:",
                     nSampMultiPresence),append = T)
  }

  # remove samples that have high abundance but low prevalence of genotyping SNVs
  # these are likely not well characterised by the original dataset
  # if the subspecies has over 30% abundance, then 75% of its genotyping SNV positions should be covered
  badSamples <- unique(all_freq_summary[all_freq_summary$median>30 & all_freq_summary$prevalence<0.75,"Sample"])
  if( length(badSamples) > 0){
    write(file=paste(outDir,species,'_extended_clustering_stat.txt',sep=''),
          x = paste0("Species ",species, ": ",
                     length(badSamples),' out of ', nrow(filtered) ,
                     ' samples rejected due to extreme mismatch between ',
                     'median abundance of genotyping SNVs (>30%) and prevalence of genotyping SNVs (<75%).')
          ,append = T)
    filtered <- filtered[!row.names(filtered) %in% badSamples,]
  }

  write.table(full,paste(outDir,species,'_extended_clustering_wFreq_unfiltered.tab',sep=''),sep='\t',quote=F)
  write.table(filtered,paste(outDir,species,'_extended_clustering_wFreq.tab',sep=''),sep='\t',quote=F)

  cluster <- data.frame(apply(filtered,1,function(x){
    assignment <- which(x>minGenotypeAbundance)
    if(length(assignment)==1){return(assignment)}
    return(NA_integer_)
  }))
  colnames(cluster) = 'clust'

  write.table(cluster,paste(outDir,species,'_extended_clustering.tab',sep=''),sep='\t',quote=F)

  return(filtered)

}

