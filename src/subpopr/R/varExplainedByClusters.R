
variationExplainedByClusters <- function(clustDf,snvFreqs, species, filePrefix, majorAllele, outDir){

  c <- clustDf$clust
  if (max(c) == 1) {#This doesn't cluster, so don't bother.
    var <- 0
    filePrefix <- paste0(filePrefix,'_noClust')
  } else {

    if(majorAllele){
      snvFreqs[snvFreqs < 50 & snvFreqs!=-1] <- 0
      snvFreqs[snvFreqs >=50] <- 1
    }

    #Total sum of squares
    sst <- apply(snvFreqs,1,function(x) {x <- x[x!=-1 & !is.na(x)]; sum((x-mean(x))^2)})
    #Replace missing with mean
    ddd1 <- t(apply(snvFreqs,1,function(x) {x1 <- x[x!=-1 & !is.na(x)]; m <- mean(x1); x[x==-1 | !is.na(x)] <- m; return(x)}))

    #Sum of squares between
    ssb <- rep(0,length(sst))
    classes = unique(c)
    for (cl in classes) {
      d1 <- ddd1[,which(c==cl),drop=FALSE]
      ssb <- ssb + sum(c==cl)*((apply(d1,1,mean)-apply(ddd1,1,mean))^2)
    }
    var <- sum(ssb)/sum(sst)*100
  }
  write.table(var,paste(outDir,filePrefix,'_variation.val',sep=''),sep='\t',quote=F,row.names=F,col.names=F)
  return(var)
}


