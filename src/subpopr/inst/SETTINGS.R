N.CORES=12
MAX.PROP.READS.NON.HOMOG = 0.05 # default = 0.1
MIN.PROP.SNV.HOMOG = 0.7 # default = 0.8
CLUSTERING.PS.CUTOFF = 0.8 # default = 0.8
SNV.SUBSPEC.UNIQ.CUTOFF = 0.8 # default = 0.8
ANALYSE.ALLELE.DISTANCES = F # default = T
BAMS.TO.USE = NULL # default = NULL
USE.PACKAGE.PREDICTION.STRENGTH = FALSE # default = FALSE

SUBPOPR.DIR="/g/scb2/bork/rossum/subspecies/testingSubpopr/" # where the subpopr package binary is located
LIB.DIR="/g/scb2/bork/rossum/R/conda/3.5.1" # set as NULL or leave blank to use default location for R libraries

PROJ.DIR="/g/scb2/bork/rossum/subspecies/ocean"

METASNV.TOP.DIR=paste0(PROJ.DIR,"/metaSNV/")
METASNV.PROJ.NAME="outputs"
METASNV.PARAM.DESC="subspecMod"
METASNV.DIR=paste0(METASNV.TOP.DIR,METASNV.PROJ.NAME,"_",METASNV.PARAM.DESC)

#SAMPLE.ID.SUFFIX=".human_gut-v11UL.uniq.coordSorted.bam"
SAMPLE.ID.SUFFIX=".fr11Repv2UL-min100bp.unique.sorted.bam"

SUBPOPR_RESULTS_DIR=paste0(PROJ.DIR,"/subpopr/fr11v2UL-human71/params.",
                           "hr",MAX.PROP.READS.NON.HOMOG*100,
                           ".hs",MIN.PROP.SNV.HOMOG*100,
                           ".ps",CLUSTERING.PS.CUTOFF*100,"/")
OUT.DIR=paste0(SUBPOPR_RESULTS_DIR,"/",METASNV.PARAM.DESC,"/")

#dir.create(path = OUT.DIR,showWarnings = F,recursive = T)

# for abundance and gene content analysis
SPECIES.ABUND.PROFILE.IS.MOTUS=TRUE
SPECIES.ABUNDANCE.PROFILE=paste0(PROJ.DIR,"/motus2/all_samples.motusv2.relabund.rmHeaderRenamed2UL.tsv")
KEGG.PATH=paste0(PROJ.DIR,"/mapped/igc/kegg.counts.all1.scaled.rmHeaderRenamed2UL.txt")
NOG.PATH=paste0(PROJ.DIR,"/mapped/igc/eggnog.counts.all1.NOG.scaled.txt")

# for phenotype association test
METADATA.COLS.TO.TEST=c("gender","geographic_location","studies","subject_disease_status")
METADATA.COL.ID="sampleNames"
METADATA.PATH=paste0(PROJ.DIR,"/metadata_manualEditReduced.csv")


DIST.METH.REPORTS="mann" # what method to use for generating reports: either "mann" or "allele"

