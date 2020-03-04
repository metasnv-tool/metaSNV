

#'if covered the whole subpop space, sum of subpop abundances as profiled using genotyping positions
#'should be near 1
#'@import dplyr
#'@param test3_output/537011_extended_clustering_wFreq.tab
#'@export
assessSubpopCompletenessPerSpecies <- function(species,subpoprOutDir){

  subpopFreqFile <- paste0(subpoprOutDir,species,"_extended_clustering_wFreq.tab")
  subpopFreq <- read.table(file = subpopFreqFile, header = T) %>% mutate(sampleName = row.names(.))

  subpopFreq %>%
    gather(key = "cluster", value="freq", -sampleName) %>%
    arrange(sampleName) %>%
    group_by(sampleName) %>%
    summarise(freqSum = sum(freq),
              nClus = n()) %>%
    #filter(freqSum != 100 ) %>%
    arrange( desc(abs(freqSum-100)) ) %>%
    mutate(species = species)

}

#' returns stats on the summed frequency of all clusters/subpops per sample
#' was e.g. lt100 = the proportion of samples where the sum of the cluster frequencies was less than 100
#'@import dplyr
#'@param subpoprOutDir directory where subpopr output is saved
#'@param species array vector of species ids
assessSubpopCompleteness <- function(speciesToAssess,subpoprOutDir){

  freqSums <- lapply(speciesToAssess,assessSubpopCompletenessPerSpecies, subpoprOutDir = subpoprOutDir) %>% do.call(rbind,.)
  freqSumsStats <- freqSums %>%
    group_by(species) %>%
    summarise( nClus = max(nClus),
               nSamples = n(),
               eq100 = length(freqSum[freqSum == 100])/nSamples,
               gt100 = length(freqSum[freqSum > 100])/nSamples,
               gt110 = length(freqSum[freqSum > 110])/nSamples,
               gt120 = length(freqSum[freqSum > 120])/nSamples,
               lt100 = length(freqSum[freqSum < 100])/nSamples,
               lt90 = length(freqSum[freqSum < 90])/nSamples,
               lt80 = length(freqSum[freqSum < 80])/nSamples,
               lt50 = length(freqSum[freqSum < 50])/nSamples
               ) %>%
#    mutate_if(is.numeric,funs(prop = ./nSamples)) %>%
    mutate(warningFlag = (eq100 < 0.8 | gt100 > 0.05 | lt90 > 0.05 | gt120 != 0 | lt50 != 0 ) )  %>%
    merge(x=.,y=TAXA.NCBI.MOTU.MAP,by.x="species",by.y="ncbiTaxID",all.x=T) %>%
    arrange(desc(warningFlag), eq100)


  return(freqSumsStats)
}


assessSubpopExclusivityPerSpecies <- function(species,subpoprOutDir){

  subpopFreqFile <- paste0(subpoprOutDir,species,"_extended_clustering_wFreq.tab")
  subpopFreq <- read.table(file = subpopFreqFile, header = T) %>% mutate(sampleName = row.names(.))

  subpopFreq %>%
    gather(key = "cluster", value="freq", -sampleName) %>%
    arrange(sampleName) %>%
    group_by(sampleName) %>%
    summarise(isExclusive = (100 %in% freq) && sum(freq) == 100 ,
              nClus = n()) %>%
    summarise(nSamples = n(),
              nClus = as.integer(max(nClus)),
              propExclusive = sum(isExclusive)/nSamples
              ) %>%
    mutate(species = species)

}


assessSubpopExclusivity <- function(species,subpoprOutDir,specIDir){

  exclu <- lapply(species,assessSubpopExclusivityPerSpecies, subpoprOutDir = subpoprOutDir) %>% do.call(rbind,.)
  excluStats <- exclu %>%
    merge(x=.,y=TAXA.NCBI.MOTU.MAP,by.x="species",by.y="ncbiTaxID",all.x=T) %>%
    arrange(propExclusive)

  return(excluStats)
}



