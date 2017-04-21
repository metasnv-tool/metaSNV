###############################################
#  metaSNP : Step III   `Post-processing`     #
###############################################

# This code is part of the metagenomic SNP calling pipeline (metaSNP)
# Copyright (c) 2016 Robin Munch
# Licenced under the GNU General Public License (see LICENSE) 
#
# This script computes SNP profile distances between sample pairs.
# - distances are computed over the SNP space.
# - default distance measure: manhattan distance.
# - PCoA visualization to gain a quick overview of the SNP space.

## INVOCATION (manually):
#
# Example: 
#   $ R --no-save --file=/path/2/src/computeDistances.R --args  infile.freq outfolder/
#
# Loop (all): 
#   for i in $(ls project_dir/filtered/pop/*.freq); \
#       do echo R --no-save --file=../src/subspeciesClustering.R \
#       --args $i project_dir/distances/  \
#   ;done
#############


## Functions

majorAllel = function(d) {
  apply(d,1,function(x) {
    x <- x[x!=-1]
    x[x <50] <- 0
    x[x >=50] <- 1
    return(median(x))
  })
}

rowMedians = function(d) {
  apply(d,1,function(x) median(x[x!=-1]))
}

rMeans = function(d) {
  apply(d,1,function(x) mean(x[x!=-1]))
}

getPnPs = function(s,d) {
  s <- as.character(s$V2)
  s <- substr(s,1,1)
  m <- apply(d,2,function(x) {
    a <- s[which(x > 50)]
    return(sum(a=='N')/sum(a=='S'))
  })
  
  return(m)
}


#############
## libraries
library(ape) #pcoa
library(ggplot2)
library(gridExtra)
#############

#############
## Arguments 
#argv <- c("/g/bork3/home/costea/R-3.0.2/bin/R", "--no-save", "--file=src/computeDistances.R", "--args", "tutorial/filtered/pop/420247.filtered.freq", "tutorial/distances/")
argv <- commandArgs(trailingOnly = F) # All arguments including --file
scriptPath <- dirname(sub("--file=","",argv[grep("--file",argv)])) ##sub(pattern, replacement, x, ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)
dist_location <- paste(scriptPath,'/distances/',sep='')

####
# Source required packages
source(paste(dist_location,'Strain_functions.R',sep=''))

args <- commandArgs(trailingOnly = TRUE) #--args
print(args)

file <- args[1]
outf <- args[2]

print(file)
print(outf)

data <- read.table(file,header=T,row.names=1,check.names=F)

species <- strsplit(sapply(strsplit(file,"\\/"),tail,1),'\\.')[[1]][1]

mantel <- NULL

print(dim(data))
colnames(data)


#### MetaSNP First (data)
#Remove mostly constant positions
freq_means <- apply(data,1,function(x) {
  a <- which(x==-1)
  if (length(a) > 0) {
    return(mean(x[-a]))
  } else {
    return(mean(x))
  }
})
#Remove mostly constant positions
r <- c(which( freq_means > 0.95),which( freq_means < 0.05))
if (length(r) > 0) {
  data <- data[-r,]
}

## Remove low resolution Samples:
lowRes <- which(colSums(data==-1)>(nrow(data)*0.5))
if (length(lowRes) > 0) {
  #Remove samples with more than half -1's
  print('Removing samples because they do not have enough resolution')
  print(length(lowRes))
  #data <- data[,-lowRes]
}

# Compute allele frequency distributional properties 
freq_data_sample_20_80 <- apply(data,2,function(x) {
  x <- x[x!=-1]
  sum(x < 0.20 | x > 0.80)/length(x)
})
freq_data_sample_10_90 <- apply(data,2,function(x) {
  x <- x[x!=-1]
  sum(x < 0.10 | x > 0.90)/length(x)
})
freq_data_sample_5_95 <- apply(data,2,function(x) {
  x <- x[x!=-1]
  sum(x < 0.5 | x > 0.95)/length(x)
})

freq_data_sample <- data.frame(freq_data_sample_20_80,freq_data_sample_10_90,freq_data_sample_5_95)
write.table(freq_data_sample,paste(outf,species,'_freq_composition.tab',sep=''),sep='\t',quote=F, col.names = NA)

# Requirement:
# !! <freq_euc.so> !! # 
#------COMPILATION--------#
# The compiling is trivial:
# R CMD COMPILE freq_euc.c
# R CMD SHLIB freq_euc.o
#------COMPILATION--------#
mann_dist <- mann.dist(data*100)
write.table(as.matrix(mann_dist),paste(outf,species,'_mann_distance_matrix.tab',sep=''),sep='\t',quote=F, col.names = NA)

allele_dist <- allel.dist(data*100)
write.table(as.matrix(allele_dist),paste(outf,species,'_allele_distance_matrix.tab',sep=''),sep='\t',quote=F, col.names = NA)


############################################################################


#### Vizualize the SNP space #####

##  PCoA Plots for each Species

#Get PCOA projection
pca <- pcoa(mann_dist)
eig <- pca[["values"]][["Eigenvalues"]]
eig[eig<0] <- 0
eig_m <- eig/sum(eig)*100

pcoa_df <- data.frame(pca[["vectors"]][,1:3])
pcoa_df$f <- freq_data_sample[rownames(pcoa_df),]$freq_data_sample_10_90
#pcoa_df$samples <- df[rownames(pcoa_df),]

write.table(pcoa_df,paste(outf,species,'_mann_pcoa_proj.tab',sep=''),sep='\t',quote=F, col.names=NA)

## Vizualize data as 2D PCoA

pdf(paste(outf,species,'_mann_PCoA.pdf',sep=''))

grid.arrange(ggplot(pcoa_df,aes(x=Axis.1,y=Axis.2)) + geom_point() + ggtitle(paste("Species:",species))
             + xlab(sprintf("PCoA1: %3.2f%%",eig_m[1])) + ylab(sprintf("PCoA2: %3.2f%%",eig_m[2]))  )

dev.off()