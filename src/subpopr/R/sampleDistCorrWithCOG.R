# writeSampleDistMatrices <- function(snvFreq,species){
#
#   mann_dist <- mann.dist(snvFreqs)
#   write.table(as.matrix(mann_dist),paste('output/',species,'_mann_distance_matrix.tab',sep=''),sep='\t',quote=F)
#
#   allele_dist <- allel.dist(snvFreqs)
#   write.table(as.matrix(allele_dist),paste('output/',species,'_allele_distance_matrix.tab',sep=''),sep='\t',quote=F)
#
# }


writeSampleDistCorrWithCOGgeneSNVs <- function(snvFreqs,mann_dist,allele_dist,cogs,species, inDir){

  cogs <- read.table(paste0(inDir,'/COG_genes.tab'),header=F,sep='\t',check.names=F,row.names=1)

  cogs$genome <- sapply(strsplit(rownames(cogs),'\\.'),'[[',1)
  cogs$length <- cogs[["V4"]]-cogs[["V3"]]

  genes_of_snps <- sapply(strsplit(rownames(snvFreqs),':'),'[[',2)
  COG_snps <- snvFreqs[which (genes_of_snps %in% rownames(cogs)),]

  cog_mann_dist <- mann.dist(COG_snps)

  m <- mantel.randtest(cog_mann_dist,mann_dist)
  mantel <- c(mantel,m[["obs"]])

  write.table(mantel,paste('extra/',species,'_mann_mantel.cor',sep=''),sep='\t',quote=F,row.names=F,col.names=F)


  cog_allele_dist <- allel.dist(COG_snps)

  m <- mantel.randtest(cog_allele_dist,allele_dist)
  mantel <- c(mantel,m[["obs"]])

  write.table(mantel,paste('extra/',species,'_allele_mantel.cor',sep=''),sep='\t',quote=F,row.names=F,col.names=F)


}
