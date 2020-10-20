
# dplyr version
#profileFilePath <- "/Users/rossum/Dropbox/PostDocBork/subspecies/pancreaticCancer/motus2/all_samples.motusv2.relabund.tsv"
# parseMotu2Profile <- function(profileFilePath){
#   rawProfile <- read.table(file = profileFilePath, col_names = T)
#
#   parsedProfile <- rawProfile %>%
#     # handle mOTUs that were not taxonomically classified
#     mutate(X1 = if_else(condition = X1 == -1 ,true = "Unclassified [Unclassified]", false = X1 )) %>%
#     tidyr::separate(col = X1,into=c("speciesName","mOTU"),sep=" \\[") %>%
#     mutate(mOTU = sub(mOTU,pattern = "]",replacement = ""))
#
#   return(parsedProfile)
# }
#

#profileFilePath <- "/Users/rossum/Dropbox/PostDocBork/subspecies/pancreaticCancer/motus2/all_samples.motusv2.relabund.tsv"
parseMotu2Profile <- function(profileFilePath){
  rawProfile <- read.delim(comment.char = "#",
                           file = profileFilePath,header = F,sep="\t",as.is = T) # header = F because "-" in sample names is translated to "."

  # get sample names
  colnames(rawProfile) <- rawProfile[1,]
  colnames(rawProfile)[1] <- "X"
  rawProfile <- rawProfile[-1,]

  rawProfile[rawProfile$X==-1,"X"] <- "Unclassified [Unclassified]"

  # format, eg: "unknown Clostridiales [meta_mOTU_v2_7800]"
  tax <- do.call(rbind,strsplit(x = rawProfile$X,split = "[ ]* \\["))
  tax <- data.frame(tax)
  if(ncol(tax)==1){ #in case the formatting was off
    tax$var2 <- tax[,1]
  }
  colnames(tax) <- c("speciesName","mOTU")
  tax$mOTU <- sub(x = tax$mOTU, pattern = "]",replacement = "")

  parsedProfile <- cbind.data.frame(tax,rawProfile[,colnames(rawProfile) != "X"])
  return(parsedProfile)
}



