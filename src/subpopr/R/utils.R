
rowMedians = function(d,na.rm = T) {
  apply(d,1,function(x) median(x[x!=-1],na.rm = na.rm))
}

rowMeans = function(d,na.rm = T) {
  apply(d,1,function(x) mean(x[x!=-1],na.rm = na.rm))
}

rowMax = function(d,na.rm = T) {
  apply(d,1,function(x) max(x,na.rm = na.rm))
}

makeNA <- function(DT) {
  DT[DT==-1] <- NA
#  for (j in seq_len(ncol(DT)))
#    set(DT,which(DT[[j]]==-1),j,NA)
  return(DT)
}

# https://github.com/tidyverse/dplyr/issues/361
sample_n_groups = function(tbl, size, replace = FALSE, weight = NULL) {
  # regroup when done
  grps = tbl %>% groups %>% lapply(as.character) %>% unlist
  # check length of groups non-zero
  keep = tbl %>% summarise() %>% ungroup() %>% sample_n(size, replace, weight)
  # keep only selected groups, regroup because joins change count.
  # regrouping may be unnecessary but joins do something funky to grouping variable
  tbl %>% right_join(keep, by=grps) %>% group_by_(.dots = grps)
}

getClustMedoidDefnFailedDir <- function(outDir){
  outDirFailed <- paste0(outDir,"/clustMedoidDefnFailed/")
  if(!dir.exists(outDirFailed)){
    dir.create(outDirFailed)
  }
  return(outDirFailed)
}


getNoClusteringDir <- function(outDir){
  outDirNoSubStructure <- paste0(outDir,"/noClustering/")
  if(!dir.exists(outDirNoSubStructure)){
    dir.create(outDirNoSubStructure)
  }
  return(outDirNoSubStructure)
}

getSnvFreqPlotDir <- function(outDir){
  outDirSnvFreqPlot <- paste0(outDir,"/snvFreqPlots/")
  if(!dir.exists(outDirSnvFreqPlot)){
    dir.create(outDirSnvFreqPlot)
  }
  return(outDirSnvFreqPlot)
}
