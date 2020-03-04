
#' @exportMethod
collectSubpopAbunds <- function(resultsDir, writeDir=NULL){
  if(is.null(writeDir)){
    writeDir <- resultsDir
  }
  allClustAbunds <- list.files(paste0(resultsDir,"/"),paste('.*hap_coverage_extended_normed.tab',sep=''),full.names = T)

  if(length(allClustAbunds)==0){
    warning("No cluster abundances to summarise.")
    return(NULL)
  }
  subpopAbunds <- data.frame()


  # for each species and each cluster .pos file
  for (d in allClustAbunds) {
    # get e.g. 537011
    species <- strsplit(basename(d) ,'_')[[1]][1]
    # get the cluster number
    cluster <- strsplit(basename(d) ,'_')[[1]][3]

    abunds <- read.table(d,header=T,row.names=1,sep='\t',col.names = c("X1"))
    subpopAbunds <- rbind.data.frame(subpopAbunds,
                                     data.frame("sampleName"=row.names(abunds),
                                                "species"=rep(species,nrow(abunds)),
                                                "subpop"=rep(cluster,nrow(abunds)),
                                                "abundance"=abunds$X1))

  }
  subpopAbunds <- subpopAbunds[order(subpopAbunds$sampleName),]
  write.table(subpopAbunds,file=paste0(writeDir,"/subpopAbunds.tsv"),sep = "\t",row.names = F,quote = F)

  return(subpopAbunds)
}
