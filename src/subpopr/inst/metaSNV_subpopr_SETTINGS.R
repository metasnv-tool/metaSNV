# for phenotype association test (optional)
METADATA.COLS.TO.TEST=c("gender","geographic_location","studies","subject_disease_status")


MAX.PROP.READS.NON.HOMOG = 0.1 # default = 0.1
MIN.PROP.SNV.HOMOG = 0.8 # default = 0.8
CLUSTERING.PS.CUTOFF = 0.8 # default = 0.8
SNV.SUBSPEC.UNIQ.CUTOFF = 0.8 # default = 0.8

LIB.DIR=NULL # set as NULL or leave blank to use default location for R libraries


SAMPLE.ID.SUFFIX=".fr11Repv2UL-min100bp.unique.sorted.bam"


DIST.METH.REPORTS="mann" # what method to use for generating reports: either "mann" or "allele"
BAMS.TO.USE = NULL # default = NULL
ANALYSE.ALLELE.DISTANCES = F
USE.PACKAGE.PREDICTION.STRENGTH = FALSE # default = FALSE

