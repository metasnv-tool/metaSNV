

# keggProfilesPath <- "/Users/rossum/Dropbox/PostDocBork/subspecies/toolDevelopment/subpopr/ES_PC/mapped/igc/incomplete/kegg/keggIncompleteRenamed.txt"
# geneAbundancePath <- keggProfilesPath
# species <- "357276"
#outDir = "/Users/rossum/Dropbox/PostDocBork/subspecies/toolDevelopment/subpopr/ES_PC/metaSNV/fr11_v1/minFilter/"
#' @exportMethod
#' @return dataframe with correlation scores
correlateSubpopProfileWithGeneProfiles <- function(species,outDir,geneAbundancePath,geneFamilyType,corrMethod="pearson"){
  # library(data.table) # for faster correlations

  allClustAbund <- read.table(paste(outDir,species,'_allClust_relativeAbund.tab',sep=''),sep='\t')

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

  # geneFamilyProfiles <- read_tsv(geneAbundancePath,comment = "#",col_names = T,
  #                                cols(.default = col_number(),X1 = col_character())
  #                                ) %>%
  #   rename(geneFamily = X1) %>%
  #   gather(key = "sampleName",value = "abundance",-geneFamily)
  geneFamilyProfiles <- read_tsv(geneAbundancePath,comment = "#",col_names = T)
  colnames(geneFamilyProfiles)[1] <- "geneFamily"
  geneFamilyProfiles <- geneFamilyProfiles %>%
    gather(key = "sampleName",value = "abundance",-geneFamily)


  # only use samples that have gene abundance and subspecies frequency data
  subspeciesSamples <- unique(allClustAbund$sampleName)
  geneAbundSamples <- unique(geneFamilyProfiles$sampleName)
  samplesToUse <- intersect(unique(allClustAbund$sampleName), # subspeciesSamples
                            unique(geneFamilyProfiles$sampleName)) # gene abund samples

  # first try to fix with most common problem -- presence of extra stuff after file name e.g. [sampleID].subspec71.unique.sorted.bam
  if(length(samplesToUse) == 0 & exists("SAMPLE.ID.SUFFIX")){
    namesWithoutSuffix <- geneFamilyProfiles$sampleName
    namesWithSuffix <- paste0(namesWithoutSuffix,SAMPLE.ID.SUFFIX)
    samplesToUse <- intersect(namesWithSuffix, subspeciesSamples)
    if(length(samplesToUse)>0){
      warning("For species '",species,"': No overlapping sample IDs between clustering and gene abundance profiles. ",
              "Fixed by adding sample suffix '",SAMPLE.ID.SUFFIX,"' to IDs in gene abundance file. ",
              "Now, out of ",
              length(subspeciesSamples)," samples with SNV data, ",
              length(samplesToUse)," also have gene abundance data.")
      geneFamilyProfiles$sampleName <- namesWithSuffix
    }
  }

  samplesToUse <- intersect(unique(allClustAbund$sampleName), # subspeciesSamples
                            unique(geneFamilyProfiles$sampleName)) # gene abund samples

  if(length(samplesToUse) == 0){
    stop(paste("No overlapping sample IDs between clustering and gene family abundance profiles.",
               " Are IDs formatted differently? Example clustering sample IDs: ",
               paste(head(n=3,unique(allClustAbund$sampleName)),collapse=", "),
               ". Example gene family sample IDs:",
               paste(head(n=3,unique(geneFamilyProfiles$sampleName)),collapse=", "), sep=""))
  }

  geneFamilyProfiles <- geneFamilyProfiles %>%
    filter(sampleName %in% samplesToUse) %>%
    group_by(sampleName) %>%
    mutate(koAbundNorm = abundance/sum(abundance)) %>% # should we do this?
    select(-abundance) %>%
    ungroup()

  allClustAbund <- filter(allClustAbund, sampleName %in% samplesToUse) %>%
    rename(clustAbund = abundance) %>%
    mutate(cluster = sub(x = cluster, pattern = "^X",replacement = "")) %>% # remove X prefix that was introduced when cluster was column name
    complete(cluster, nesting(sampleName), fill = list(clustAbund = 0)) #add 0 abundances for clusters not seen in samples

  nClus <- allClustAbund$cluster %>% unique() %>% length()
  print(paste0("Species ",species,": Testing ",nClus," clusters for correlation with ",
               geneFamilyProfiles %>% filter(koAbundNorm>0 & geneFamily !="-1") %>%
                 pull(geneFamily) %>% unique() %>% length()," geneFamily groups in ",
               length(unique(allClustAbund$sampleName))," samples"))

  fullData <- full_join(geneFamilyProfiles, allClustAbund,by="sampleName") %>% as.data.table()
  # for speed optimisation later
  # print(system.time(c <- left_join(geneFamilyProfiles, allClustAbund,by="sampleName")))
  # print(system.time(a <- as.data.table(geneFamilyProfiles)[as.data.table(allClustAbund),on="sampleName",allow.cartesian=TRUE])) # this is the fastest, but not by much
  # print(system.time(b <- as.data.table(allClustAbund)[as.data.table(geneFamilyProfiles),on="sampleName",allow.cartesian=TRUE]))

  cor.test.mod <- function(x,y,method,exact,geneFamily,cluster){
    # if no variance in either value, don't compute correlation
    if(diff(range(x)) == 0){return(NULL)}
    if(diff(range(y)) == 0){return(NULL)}

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

  corr <- fullData[, cor.test.mod(koAbundNorm,clustAbund,method=corrMethod,exact = !(corrMethod=="spearman"),geneFamily,cluster ),
                    by = c("geneFamily","cluster")]

  corr <- select(corr, -data.name, -parameter)
  corr$q.valueBH <- p.adjust(corr$p.value,method = "BH")
  write_delim(x = corr,path = paste0(outDir,"/",species,"_corr",geneFamilyType,"-",corrMethod,".tsv"),delim = "\t")
  write_delim(x = corr[q.valueBH<0.1],path = paste0(outDir,"/",species,"_corr",geneFamilyType,"-",corrMethod,"-qSig1.tsv"),delim = "\t")

  return(corr)
}
