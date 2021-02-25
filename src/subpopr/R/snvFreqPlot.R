
snvFreqPlot<-function(species,snvFreqs,outDir,
                      minPropHomogSnvAllelesPerSample=0.8,
                      maxPropReadsNonHomog = 0.1){

  maxPropReadsNonHomog = maxPropReadsNonHomog * 100

  cutoffPropSNV <- minPropHomogSnvAllelesPerSample
  cutoffHomog <- maxPropReadsNonHomog

  cutoffsLow <- seq(from = 0, to = 49, by = 1)
  cutoffsHigh <- seq(from = 100, to = 50, by = -1)
  cutoffHomogIndex <- cutoffHomog + 1

  totalSNVs <- nrow(snvFreqs)

  getPropPassingCutoff <- function(i, direction){
    cutoffIndex <- combos[i,1]
    sampleID <- combos[i,2]
    freqs <- snvFreqs[,sampleID,drop=T]
    freqs <- freqs[freqs > -1 & !is.na(freqs) ] # SNVs with too low coverage are ignored in calc
    if(direction == "low"){
      prop <- sum(freqs <= cutoffsLow[cutoffIndex])/length(freqs)
    }else{
      prop <- sum(freqs >= cutoffsHigh[cutoffIndex])/length(freqs)
    }
    return(prop)
  }

  getPropSuffCoverage <- function(i){
    cutoffIndex <- combos[i,1]
    sampleID <- combos[i,2]
    freqs <- snvFreqs[,sampleID,drop=T]
    prop <- sum(freqs > -1 & !is.na(freqs)) / length(freqs) # SNVs with too low coverage are ignored in calc
    return(prop)
  }

  combos <- expand.grid(cutoffIndex=1:length(cutoffsLow),
                        sampleID=colnames(snvFreqs),
                        stringsAsFactors = F)

  propLow <- sapply(1:nrow(combos),getPropPassingCutoff,direction="low")
  propHigh <- sapply(1:nrow(combos),getPropPassingCutoff,direction="high")
  propSuffCov <- sapply(1:nrow(combos),getPropSuffCoverage)

  cutoffData <- cbind.data.frame(combos,propLow=propLow,propHigh=propHigh,propSuffCov=propSuffCov)
  cutoffData$propPass <- cutoffData$propLow + cutoffData$propHigh

  cutoffData$sampleUsedForMedoidDefn <- ((cutoffData$cutoffIndex) == cutoffHomogIndex & cutoffData$propPass > cutoffPropSNV)
  cutoffData$sampleUsedForMedoidDefn[cutoffData$cutoffIndex != cutoffHomogIndex] <- NA # don't highlight at non cut-off x values
  #current cutoff: 80% of SNVs must be seen in <5% or >95% of reads in the sample

  getPlot<-function(y="propPass"){
    df <- cutoffData %>%
      mutate(x=cutoffIndex-1)
    df %>%
      ggplot(aes_string(x="x", y = y, group="sampleID")) +
      ylim(c(-0.05,1.05))+
      geom_vline(xintercept = cutoffHomog, color = "grey60",linetype="dotted") +
      geom_hline(yintercept = cutoffPropSNV, color = "maroon",linetype="dotted") +
      geom_point(size=0.3)+
      geom_line(alpha=0.5, aes(color= propSuffCov))+
      scale_color_continuous("Proportion of SNV positions \nwith sufficient coverage", low = "grey50")+ #lightgreen
      #geom_point(shape=21,data= filter(df,sampleUsedForMedoidDefn),
      #           aes(fill=sampleUsedForMedoidDefn))+
      geom_point(size=3, alpha=0.6,aes_string(shape="sampleUsedForMedoidDefn",fill=y))+
      scale_shape_manual("Sample used for \nsubspecies discovery",values=c("TRUE"=21,"FALSE"=24) ) +
      #scale_fill_manual("Sample used for subspecies discovery",values=c("TRUE"="maroon","FALSE"="grey50") ) +
      scale_fill_gradientn(values = c(0,cutoffPropSNV,1),
                                    limits=c(0,1),
                                    colors=c("firebrick4","lightgoldenrod1", "darkblue"),
                                    #guide = guide_colourbar(barwidth = 10),
                                    name="Proportion of SNVs \nwhere allele is \napprox. \"homogeneous\" ")+
      theme_bw()+theme(legend.position="bottom",legend.box="horizontal")+
      guides(fill = guide_colourbar(title.position="top", title.hjust = 0),
             size = guide_legend(title.position="top", title.hjust = 0),
             shape = guide_legend(title.position="top", title.hjust = 0),
             colour = guide_colourbar(title.position="top", title.hjust = 0,
                                      override.aes = list(alpha = 0.1)))
  }
  pAll <- getPlot("propPass")+
      xlab("% away from homogenity at SNV (e.g. x = 5: allele is seen in < 5% or > 95% of reads)") +
      #ylab("Proportion of SNVs with the same allele in [< x] or [> 100-x] % of reads (pass the cutoff x)")+
      ylab("Proportion of SNV positions with a (nearly) \"homogeneous\" allele \n (across almost all [ >100-x] or none [ < x ] of reads)")+
      ggtitle(paste0("What proportion of SNVs are 'almost' homogeneous for this species?\n (",species,") Each line is a sample.",
                     " Total # SNVs: ",totalSNVs))

  # pLow <- getPlot("propLow") +
  #   xlab("% away from homogenity at SNV (e.g. x = 5: allele is seen in < 5% of reads)") +
  #   ylab("Proportion of SNVs with the same allele in < x% of reads (pass the cutoff x)")+
  #   ggtitle(paste0("What proportion of SNVs are 'almost' *never* seen for this species?\n (",species,") Each line is a sample.",
  #                  " Total # SNVs: ",totalSNVs))
  #
  # pHigh <- getPlot("propHigh") +
  #   xlab("% away from homogenity at SNV (e.g. x = 5: allele is seen in > 95% of reads)") +
  #   ylab("Proportion of SNVs with the same allele in > 100-x % of reads (pass the cutoff x)")+
  #   ggtitle(paste0("What proportion of SNVs are 'almost' *always* seen for this species?\n (",species,") Each line is a sample.",
  #                  " Total # SNVs: ",totalSNVs))

  ggsave(pAll,filename = paste0(outDir,"/",species,"_snvFreq_HighOrLow.png"),units = "in",width = 7,height = 7,dpi = 200)
  #ggsave(pHigh,filename = paste0(outDir,"/",species,"_snvFreq_High.png"),units = "in",width = 7,height = 7,dpi = 200)
  #ggsave(pLow,filename = paste0(outDir,"/",species,"_snvFreq_Low.png"),units = "in",width = 7,height = 7,dpi = 200)


  pSum <- cutoffData %>% filter((cutoffData$cutoffIndex) == cutoffHomogIndex) %>%
    ggplot(aes(x=propPass))+
    geom_histogram(binwidth = 0.05)+
    xlab("Within a sample: \n% SNV positions where allele is \"fixed\"")+
    ylab("Number of samples")+theme_bw()+
    geom_vline(xintercept = cutoffPropSNV, color = "maroon",linetype="dotted")+
    scale_x_continuous(labels = function(x) paste0(x*100, "%"),limits = c(-0.05,1.05))

  ggsave(pSum,filename = paste0(outDir,"/",species,"_snvFreqFixedHist.png"),width = 4,height = 2,units = "in",dpi = 300 )

}
