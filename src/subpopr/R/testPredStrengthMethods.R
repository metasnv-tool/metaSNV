library(ggplot2)
library(dplyr)
library(tidyr)


getNclus <- function(i, params, drawHeatmap=FALSE,nOutliersToAdd=NULL){
  nSamp <- params[i,1]
  diffBetweenWithin <-params[i,2]
  sizeClust1 <- ceiling(nSamp*0.6)
  sizeClust2 <- nSamp - sizeClust1

  rangeDelim <- (1-diffBetweenWithin)/2
  dissimBetweenClusters <- runif(n = 1000,min = 1-rangeDelim,max = 1) # high values
  dissimWithinClusters <- runif(n = 1000,min = 0.001,max = rangeDelim) # low values

  topLeft <- matrix(sample(x = dissimWithinClusters,size = (sizeClust1*sizeClust1),replace = T),ncol = sizeClust1)
  bottomRight <- matrix(sample(x = dissimWithinClusters,size = (sizeClust2*sizeClust2),replace = T),ncol = sizeClust2)

  topRight <- matrix(sample(x = dissimBetweenClusters,size = (sizeClust1*sizeClust2),replace = T),nrow = sizeClust1)
  bottomLeft <- matrix(sample(x = dissimBetweenClusters,size = (sizeClust2*sizeClust1),replace = T),nrow = sizeClust2)

  # this doesn't give a symmatric matrix, but I'm just taking one triangle
  # when I make it a dist obj so it's ok
  x <- rbind(cbind(topLeft,topRight),
             cbind(bottomLeft,bottomRight))

  if(!is.null(nOutliersToAdd) & nOutliersToAdd>0){
    for(i in 1:nOutliersToAdd){
      x<-rbind(x,sample(x = dissimBetweenClusters,size = ncol(x),replace = T))
      x<-cbind(x,sample(x = dissimBetweenClusters,size = nrow(x),replace = T))
    }
  }

  diag(x)<-0 # doesn't change anything
  distObj <- as.dist(x)

  res<-list()
  resCustom <- predStrengthCustom(distObj,Gmin = 2,Gmax = 10,cutoff = 0.8,
                     clustermethod = cluster::pam)

  resPackage <- fpc::prediction.strength(xdata = distObj,distances = T,cutoff = 0.8,Gmax = 10,
                          clustermethod = fpc::claraCBI, #pamLike=T, error
                          Gmin = 2, M = 50,classification = "centroid")

  if(drawHeatmap){
    pheatmap::pheatmap(as.matrix(distObj),border_color = NA,
                       cluster_rows = FALSE,cluster_cols = FALSE)
  }

  #fclust::Fclust(X = as.matrix(distObj),distance=T,k=2)
  #fclust::Fclust(X = as.matrix(distObj),distance=T,k=3)

  res[["custom"]] <- resCustom$predcorr
  res[["package"]] <- resPackage$predcorr

  nClus <- list(nclus_custom=resCustom$optimalk,nclus_package=resPackage$optimalk)

  return(nClus)
}

runTests<-function(){

  nSampleTrys <- seq(20,150,by=10) # number of samples
  diffBetweenWithinTrys <- seq(-0.4,0.8,by=0.1) #Size of gap in dissimilarities between clusters

#res1 <- sapply(169:180,getNclus)

params <- expand.grid(nSamples=nSampleTrys,
                      distinctiveness=diffBetweenWithinTrys)

res1_withOutliers <- sapply(1:nrow(params),getNclus,params = params, nOutliersToAdd = 3)
res1_noOutliers <- sapply(1:nrow(params),params = params, getNclus)
res1<-res1_withOutliers

res <- data.frame(t(res1),stringsAsFactors = F)
res <- cbind(params,res)

res <- tidyr::gather(res,key = "method",value = "nClus",-nSamples,-distinctiveness)
res$nSamples <- as.numeric(res$nSamples)
res$nClus <- as.numeric(res$nClus)

res %>% mutate(NumClusters=factor(nClus,ordered = T)) %>%
  ggplot(aes(x=nSamples,y=method,
               fill=NumClusters))+
  geom_tile(height=0.5,color="grey50")+
  #geom_text(aes(label=nClus))+
  facet_grid(distinctiveness~.)+
  theme(legend.position = "top")


resThresh <- res %>%
  group_by(method,distinctiveness) %>%
  filter(nClus == 2) %>% slice(1)

resThresh %>% ggplot(aes(x=distinctiveness,y=nSamples,color=method))+
  geom_point()+geom_line()+ylim(c(0,NA))+
  xlim(c(min(diffBetweenWithinTrys),max(diffBetweenWithinTrys)))+
  xlab("Size of gap in dissimilarities between clusters")+
  ylab("Minimum number of samples \nto get correct number of clusters")

# outliers really throw both methods off, unless there are lots of samples
params[172,]
getNclus(i = 172,drawHeatmap = F,nOutliersToAdd = 0)
getNclus(i = 172,drawHeatmap = F,nOutliersToAdd = 3)


# From the plots above, looks like we need ~100 samples

# Which method is better at dealing with outliers?
params[177,]
paramsBigSamples <- params[params$nSamples>100,]
getNclus(i = 179,params = paramsBigSamples, drawHeatmap = T,nOutliersToAdd = 3)
getNclus(i = 179,params = paramsBigSamples, drawHeatmap = T,nOutliersToAdd = 10)

}
