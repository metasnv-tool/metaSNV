
#' scale the abundance of subpops by the abundance of the species
#' @param species the numeric species ID
#' @param clusterFreqs frequencies of the clusters/subpops,
#' with each row a sample and each column a cluster/subpop. Output from [writeSubpopsForAllSamples()] (species,"_extended_clustering_wFreq.tab")
#' @param motuProfileFilePath path to motu profiler relative abundance output e.g. all_samples.motusv2.relabund.tsv
writeSubpopAbundMotusProfile <- function(species, clusterFreqs, outDir,motuProfileFilePath){

  #Read the mOTU abundance table
  #columns are samples, rows are specI clusters, values are abundance
  motuProfile <- parseMotu2Profile(motuProfileFilePath)

  #accomodate package and non-package
  taxaMap <- ifelse("subpopr" %in% tolower(.packages()),
                    subpopr::TAXA.NCBI.MOTU.MAP, # from /data folder
                    TAXA.NCBI.MOTU.MAP)

  # get the specI cluster ID / mOTU IDs that corresponds to the species being analysed
  ncbiToMotus <- taxaMap %>%
    filter(ncbiTaxID == species) %>%
    select(ncbiTaxID, ref_mOTU_cluster) %>%
    distinct()

  if(nrow(ncbiToMotus) > 1){
    stop(paste0("More than one mOTU match for species ", species,". mOTUs: ", paste(ncbiToMotus$ref_mOTU_cluster,collapse = ",")))
  }
  # get the corresponding mOTU id
  motuID <- ncbiToMotus$ref_mOTU_cluster[1]

  if( ! motuID %in% motuProfile$mOTU){
    stop(paste0("Species not found in mOTU profile. Species: ",species,"; mOTU id: ",motuID))
  }
  # get the species' mOTU profile
  motuProfileThisSpecies <- motuProfile[motuProfile$mOTU == motuID,,drop=T]

  if(nrow(clusterFreqs) == 0){
    stop(paste0("Species does not have cluster frequencies. Species: ",species))
  }

  # only use samples that have species abundance and subspecies frequency data
  samplesToUse <- intersect(names(motuProfileThisSpecies), rownames(clusterFreqs))
  if(length(samplesToUse) == 0){
    stop(paste("No overlapping sample IDs between clustering and mOTU abundance profiles. Are IDs formatted differently? Example clustering IDs: ",
               paste(tail(n=3,rownames(clusterFreqs)),collapse=", "),
               ". Example mOTU profile IDs:",
               paste(tail(n=3,names(motuProfileThisSpecies)),collapse=", "), sep=""))
  }

  #Subset subpop/cluster data to samples that also have taxonomic profile
  clusterFreqs <- clusterFreqs[samplesToUse,]
  motuProfileThisSpecies <- motuProfileThisSpecies[samplesToUse]
  motuProfileThisSpecies <- as.numeric(unlist(motuProfileThisSpecies))

  #Get the coverages and scale with subpop frequency
  allClustAbund <- clusterFreqs / 100 * motuProfileThisSpecies
  write.table(allClustAbund,paste(outDir,species,'_allClust_relativeAbund.tab',sep=''),sep='\t',quote=F)

  # for each cluster/subpop (for reverse compatibility)
  for (x in 1:ncol(clusterFreqs)) {
    clusterNFreq <- clusterFreqs[,x,drop=F]

    #Get the coverages and scale with subpop frequency
    clustAbund <- clusterNFreq[1] / 100 * motuProfileThisSpecies

    write.table(clustAbund,paste(outDir,species,'_clust_',x,'_hap_coverage_extended_normed.tab',sep=''),sep='\t',quote=F)
  }

}

#' scale the abundance of subpops by the abundance of the species
#' @param species the numeric species ID
#' @param clusterFreqs frequencies of the clusters/subpops,
#' with each row a sample and each column a cluster/subpop. Output from [writeSubpopsForAllSamples()] (species,"_extended_clustering_wFreq.tab")
#' @param speciesProfileFilePath path to profile of relative abundances of spcies,
#' 1st columns must have names that match the names of the species used in the
#' reference database the samples were mapped against for SNV calling
writeSubpopAbundSpeciesAbund <- function(species, clusterFreqs, outDir,speciesProfileFilePath){

  #profileFilePath <- "/Volumes/KESU/scb2/bork/rossum/subspecies/testingSubpopr/inSilicoMock/mine/smallerTestSet/smallGenomes/01_createData/outputs/speciesAbundances.tsv"
  #clusterFreqs <- read.table("/Volumes/KESU/scb2/bork/rossum/subspecies/testingSubpopr/inSilicoMock/mine/smallerTestSet/smallGenomes/04_subpopr_nonPack/params.hr5.hs80.ps80/defaults/refGenome2clus_extended_clustering_wFreq.tab")
  #Read the  abundance table
  #columns are samples, rows are species, values are abundance
  rawProfile <- read.delim(comment.char = "#",
                           file = speciesProfileFilePath,header = F,sep="\t",as.is = T) # header = F because "-" in sample names is translated to "."
  colnames(rawProfile) <- rawProfile[1,]
  colnames(rawProfile)[1] <- "species"
  rawProfile <- rawProfile[-1,]

  if( ! species %in% rawProfile$species){
    stop(paste0("Species not found in the species abundance profile. Species: ",species))
  }

  speciesProfile <- rawProfile[rawProfile$species == species, ,drop=F] # drop=T]

  if(nrow(speciesProfile) > 1){
    stop("More than one species in the abundance profile matched for species named: ", species)
  }

    if(nrow(clusterFreqs) == 0){
    stop(paste0("Species does not have cluster frequencies. Species: ",species))
  }

  # only use samples that have species abundance and subspecies frequency data
  samplesToUse <- intersect(names(speciesProfile), rownames(clusterFreqs))
  if(length(samplesToUse) == 0){
    stop(paste("No overlapping sample IDs between clustering and apecies abundance profiles. Are IDs formatted differently? Example clustering IDs: ",
               paste(tail(n=3,rownames(clusterFreqs)),collapse=", "),
               ". Example species profile IDs:",
               paste(tail(n=3,names(speciesProfile)),collapse=", "), sep=""))
  }

  #Subset subpop/cluster data to samples that also have taxonomic profile
  clusterFreqs <- clusterFreqs[samplesToUse,]
  speciesProfile <- rawProfile[rawProfile$species == species, ,drop=T]
  speciesProfile <- speciesProfile[samplesToUse]
  speciesProfile <- as.numeric(unlist(speciesProfile))

  #Get the coverages and scale with subpop frequency
  allClustAbund <- clusterFreqs / 100 * speciesProfile
  write.table(allClustAbund,paste(outDir,species,'_allClust_relativeAbund.tab',sep=''),sep='\t',quote=F)

  # for each cluster/subpop (for reverse compatibility)
  for (x in 1:ncol(clusterFreqs)) {
    clusterNFreq <- clusterFreqs[,x,drop=F]

    #Get the coverages and scale with subpop frequency
    clustAbund <- clusterNFreq[1] / 100 * speciesProfile

    write.table(clustAbund,paste(outDir,species,'_clust_',x,'_hap_coverage_extended_normed.tab',sep=''),sep='\t',quote=F)
  }

}