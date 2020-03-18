
getTaxaMap <- function(){
  # accomodate package and non-package
  # ifelse() does not work here because returns same shape as test
  if("subpopr" %in% tolower(.packages())){
    taxaMap <- subpopr::TAXA.NCBI.MOTU.MAP
  }else{
    taxaMap <- TAXA.NCBI.MOTU.MAP
  }
}


#'
#'@return data frame with 2 columns: key & value
#'
getSpeciesTaxonomy <- function(speciesID){

  taxaMap <- getTaxaMap()

  #get taxa info
  if(speciesID %in% taxaMap$ncbiTaxID){
    taxonomyDf <- taxaMap %>% filter(ncbiTaxID == speciesID) %>%
      select(ncbiTaxID,everything(),-isRep) %>%
      gather() %>%
      distinct()
    return(taxonomyDf)
  }else if(speciesID %in% taxaMap$ref_mOTU_cluster){
    taxonomyDf <- taxaMap %>% filter(ref_mOTU_cluster == speciesID) %>%
      select(ref_mOTU_cluster,everything(),-isRep) %>%
      gather() %>%
      distinct()
    return(taxonomyDf)
  }else{
    return(NULL)
  }
}

getSpeciesName <- function(speciesID){
  taxonomyDf <- getSpeciesTaxonomy(speciesID)
  if(is.null(taxonomyDf)){
    return(paste("Taxonomy unknown for ",speciesID))
  }else{
    return(taxonomyDf[taxonomyDf$key=="refOTUs","value",T])
  }
}
