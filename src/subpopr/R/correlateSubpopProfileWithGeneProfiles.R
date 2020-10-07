# geneAbundancePath <- "/g/scb2/bork/rossum/metagenomes/human/subspecGeoValidation/all_v3/geneContent/mapToPanGenomes/outputs/counts_unique_norm_sumByNog.tsv"
# geneFamilyType <- "TMP1"
# species <- "742765"
# SAMPLE.ID.SUFFIX<-".subspec71.unique.sorted.bam"
# outDir = "/g/scb2/bork/rossum/metagenomes/human/subspecGeoValidation/all_v3/subpopr/results_md_s257745/params.hr10.hs80.ps80.gs80/outputs/"
#' @exportMethod
#' @return number of specifically correlated genes
correlateSubpopProfileWithGeneProfiles <- function(species,outDir,geneAbundancePath,
                                                   geneFamilyType){
  # library(data.table) # for faster correlations
  allClustAbund <- read.table(paste(outDir,species,'_allClust_relativeAbund.tab',sep=''),sep='\t')
  speciesAbund <- rowSums(allClustAbund)
  
  allClustAbund <- mutate(allClustAbund, sampleName = rownames(allClustAbund))

  # remove clusters only see in less than 3 samples
  allClustAbundInMin3Samples <- allClustAbund %>%
    gather(key = "cluster",value = "abundance",-sampleName)%>%
    filter(abundance > 0) %>%
    group_by(cluster)  %>% filter( n()>= 3) %>% ungroup()

  if(nrow(allClustAbundInMin3Samples) == 0  ){
    print(paste0("No subspecies seen in >= 3 samples for species:",species))
    return(NULL)
  }

  allClustAbund <- allClustAbundInMin3Samples
  
  # figure out which samples we need and only load data for those samples
  subspeciesSamples <- unique(allClustAbund$sampleName)
  
  # get sample names from the gene family profiles
  geneFamilyColumnNames <- read_tsv(file = geneAbundancePath,n_max=1,
                     comment = "#",trim_ws = T,col_names = T) %>% colnames()
  geneAbundSamples <- geneFamilyColumnNames[-1]
  
  # only use samples that have gene abundance and subspecies frequency data
  samplesToUse <- intersect(unique(allClustAbund$sampleName), # subspeciesSamples
                            unique(geneAbundSamples)) # gene abund samples
  geneAbundColNamesToUse <- samplesToUse
  # first try to fix with most common problem -- presence of extra stuff after file name e.g. [sampleID].subspec71.unique.sorted.bam
  if(length(samplesToUse) == 0 & exists("SAMPLE.ID.SUFFIX")){
    namesWithoutSuffix <- geneAbundSamples
    namesWithSuffix <- paste0(namesWithoutSuffix,SAMPLE.ID.SUFFIX)
    samplesToUse <- intersect(namesWithSuffix, subspeciesSamples)
    if(length(samplesToUse)>0){
      warning("For species '",species,"': No overlapping sample IDs between clustering and gene abundance profiles. ",
              "Fixed by adding sample suffix '",SAMPLE.ID.SUFFIX,"' to IDs in gene abundance file. ",
              "Now, out of ",
              length(subspeciesSamples)," samples with SNV data, ",
              length(samplesToUse)," also have gene abundance data.")
      names(geneAbundSamples) <- namesWithSuffix
      geneAbundColNamesToUse <- geneAbundSamples[samplesToUse]
      newColumnNames <- namesWithSuffix
      names(newColumnNames) <- namesWithoutSuffix
      fixColNames <- T
      #geneFamilyProfiles$sampleName <- namesWithSuffix
    }
  }
  
  if(length(samplesToUse) == 0){
    stop(paste("No overlapping sample IDs between clustering and gene family abundance profiles.",
               " Are IDs formatted differently? Example clustering sample IDs: ",
               paste(head(n=3,unique(allClustAbund$sampleName)),collapse=", "),
               ". Example gene family sample IDs:",
               paste(head(n=3,unique(geneAbundSamples)),collapse=", "), sep=""))
  }
  
                     
  f <- function(x,pos){ 
    # only keep the samples (columns) that have this species
    x <- x[,colnames(x) %in% c(colnames(x)[1],geneAbundColNamesToUse) ] 
    # only keep genes that are present
    filter(x,rowSums(x[,-1])>0)
  }
  geneFamilyProfiles <- read_tsv_chunked(file = geneAbundancePath,
                                                 comment = "#",trim_ws = T,col_names = T,
                                                 chunk_size = 10000,
                                                 callback = DataFrameCallback$new(f)
  ) %>% rename(geneFamily = 1) # rename the first column  
  
  
  # geneFamilyProfiles <- read_tsv(geneAbundancePath,comment = "#",col_names = T) %>%
  #   rename(geneFamily = 1) %>%
  #   gather(key = "sampleName",value = "abundance",-geneFamily)

  if(fixColNames){
    colnames(geneFamilyProfiles)[-1] <- newColumnNames[colnames(geneFamilyProfiles)[-1]]  
  }
  
  geneFamilyProfiles <- geneFamilyProfiles %>%
    gather(key = sampleName, value = abundance, -geneFamily) %>% 
    group_by(sampleName) %>%
    mutate(geneAbundNorm = abundance/sum(abundance)) %>%
    select(-abundance) %>%
    ungroup()
  
  allClustAbund <- filter(allClustAbund, sampleName %in% samplesToUse) %>%
    rename(clustAbund = abundance) %>%
    mutate(cluster = sub(x = cluster, pattern = "^X",replacement = "")) %>% # remove X prefix that was introduced when cluster was column name
    complete(cluster, nesting(sampleName), fill = list(clustAbund = 0)) #add 0 abundances for clusters not seen in samples

  nClus <- allClustAbund$cluster %>% unique() %>% length()
  print(paste0("Species ",species,": Testing ",nClus," clusters for correlation with ",
               geneFamilyProfiles %>% filter(geneAbundNorm>0 & geneFamily !="-1") %>%
                 pull(geneFamily) %>% unique() %>% length()," geneFamily groups in ",
               length(unique(allClustAbund$sampleName))," samples"))

  # use data table for faster corr computation
  # this is fast enough but could re-write to avoid this join and use less memory
  fullData <- inner_join(geneFamilyProfiles, allClustAbund,by="sampleName") %>% as.data.table()
  
  pseudocount <- min(geneFamilyProfiles$geneAbundNorm[geneFamilyProfiles$geneAbundNorm>0],na.rm = T)/1000

  cor.test.mod <- function(x,y,method,exact,geneFamily,cluster){
    # if no variance in either value, don't compute correlation
    if(diff(range(x)) == 0){return(NULL)}
    if(diff(range(y)) == 0){return(NULL)}

    # log 10 transform for pearson to reduce weight of rare, large values
    if(method == "pearson"){
      x <- log10(x+pseudocount)
      y <- log10(y+pseudocount)
    }
    
    res <- cor.test(x,y,method=method,exact=exact)
    if(!is.numeric(res$statistic)) {return(NULL)}
    res$parameter <- F # NULL causes error for data.table
    res$data.name <- F # save space
    res$method <- method
    if(method == "pearson"){
      if(is.null(res$conf.int)){
        # if conf int is not defined for one of the groups (e.g. if n == 3)
        # then data.table will fail b/c it will expect numeric values for conf.int.low and high
        res$conf.int <- c(as.numeric(NA),as.numeric(NA))
      }
      res$conf.int.low <- res$conf.int[1]
      res$conf.int.high <- res$conf.int[2]
      res$conf.int <- F # otherwise returns 2 rows - one each for high and low conf value
    }
    res$nObs <- length(x)
    return(res)
  }
  
  doAndSaveCorr <- function(corrMethod){
      corr <- fullData[, cor.test.mod(geneAbundNorm,clustAbund,
                                  method=corrMethod,
                                  exact = !(corrMethod=="spearman"),
                                  geneFamily,cluster ),
                    by = c("geneFamily","cluster")]
      
      corr <- select(corr, -data.name, -parameter)
      corr$q.valueBH <- p.adjust(corr$p.value,method = "BH")
      write_delim(x = corr,path = paste0(outDir,"/",species,"_corr",geneFamilyType,"-",corrMethod,".tsv"),delim = "\t")
      write_delim(x = corr[q.valueBH<0.05],path = paste0(outDir,"/",species,"_corr",geneFamilyType,"-",corrMethod,"-qSig.tsv"),delim = "\t")
      return(corr)
  }
  #geneFamilyType <- "TMP"
  corrsList <- map(c("spearman"="spearman","pearson"="pearson"),doAndSaveCorr)
  
  subspeciesSpecificGenesDf <- selectSubspeciesSpecificGenes(corrP = corrsList$pearson %>% as_tibble(), 
                                                                  corrS = corrsList$spearman %>% as_tibble())
  
  write_delim(x = subspeciesSpecificGenesDf,path = paste0(outDir,"/",species,"_corr",geneFamilyType,"-clusterSpecificGenes.tsv"),delim = "\t")
  
  if(nrow(subspeciesSpecificGenesDf) == 0){
    flog.info(paste("No cluster-specific genes found for species",speciesID))
    return(0)
  }
  specificGenes <- unique(subspeciesSpecificGenesDf$geneFamily)
  fullData %>% 
    filter(geneFamily %in% specificGenes) %>% 
    write_delim(path = paste0(outDir,"/",species,"_corr",geneFamilyType,"-clusterSpecificGeneAbundances.tsv"),delim = "\t")
  
  return(length(specificGenes))
}


#' Filter the correlation data to detect genes that are specific to (one or more) subspecies
#'@param corrP the pearson correlation results from correlateSubpopProfileWithGeneProfiles()
#'@param corrS the spearman correlation results from correlateSubpopProfileWithGeneProfiles()
#'@param minObs the minimal number of data points to consider a correlation ok
#'@param statCutoff the cut off value for q value
#'@param maxBadCorrR the highest correlation R value a non-match subspecies can have
#'@param minCorrPearson the minimal Pearson correlation R value a matched subspecies can have
#'@param minCorrSpearman the minimal Spearman correlation R value a matched subspecies can have
selectSubspeciesSpecificGenes <- function(corrP,corrS, minObs=10, statCutoff=0.05,
                                                 maxBadCorrR=0.2, 
                                                 minCorrPearson = 0.8, minCorrSpearman = 0.6){
  
   #corrP = read_tsv("/g/scb2/bork/rossum/metagenomes/human/subspecGeoValidation/all_v3/subpopr/results_md_s257745/XX_params.hr10.hs80.ps80.gs80_FULL-RESULTS/outputs/457424_corrGenes-pearson.tsv")
   #corrS = read_tsv("/g/scb2/bork/rossum/metagenomes/human/subspecGeoValidation/all_v3/subpopr/results_md_s257745/XX_params.hr10.hs80.ps80.gs80_FULL-RESULTS/outputs/457424_corrGenes-spearman.tsv") 
  # maxBadCorrR = 0.2; minObs=10; statCutoff=0.05;maxBadCorrR=0.2;minCorrPearson = 0.8; minCorrSpearman = 0.6
  
  #merge spearman and pearson results
  corrS$conf.int <- NA
  corrS$conf.int.low <- NA
  corrS$conf.int.high <- NA
  corr <- rbind(corrP,corrS) %>% 
    mutate(statSigLgl = q.valueBH < statCutoff,
           statisticallySignificant = if_else(q.valueBH < statCutoff,true = paste0("q < ",statCutoff),false = paste0("q >= ",statCutoff)),
           #cluster = paste0("Cluster ",cluster),
           cluster=factor(cluster)) %>% 
    rename(correlationR = estimate)

  # - q value < 0.05 for both Spearman and Pearson
  # - minimum of 10 observations (data points)
  # - correlation R statistic > 0.8 for Pearson
  # - correlation R statistic > 0.6 for Spearman
  # - correlated with >= 1 subspecies
  # - correlation R statistic for at least one other subspecies < 0.2 for both Spearman and Pearson
  # - for every subspecies, the gene is either strongly correlated or strongly not correlated 
  
  # min Pearson correlation 0.8 - computed on log10‐transformed relative abundances
  # min Spearman correlation 0.6
  
  subspeciesSpecificGenesDf <- corr %>% 
    mutate(corrPass = case_when( method == "pearson" & correlationR >= minCorrPearson & statSigLgl & nObs >= minObs ~ TRUE,
                                    method == "spearman" & correlationR >= minCorrSpearman & statSigLgl & nObs >= minObs ~ TRUE,
                                    TRUE ~ FALSE), 
          corrVeryBad = correlationR < maxBadCorrR ) %>% 
    group_by(geneFamily, cluster) %>% 
    summarise(geneIsCorrelated = all(corrPass),
              geneIsNotCorrelated = all(corrVeryBad)) %>% 
    ungroup() %>% 
    group_by(geneFamily) %>% 
    # correlation is specific to one or more subspecies (but not to all)
    filter( xor(geneIsCorrelated,geneIsNotCorrelated) & # gene should be strongly correlated or definitely not correlated with all subspecies
            sum(geneIsCorrelated) >= 1 & # at least 1 subspecies should pass
            sum(geneIsNotCorrelated) >= 1 , # at least 1 subspecies should fail
           .preserve = T) %>% 
    ungroup()
    
  # return the number of genes
  return(subspeciesSpecificGenesDf)
}




