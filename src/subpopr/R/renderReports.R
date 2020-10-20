
renderDetailedSpeciesReport <- function(speciesID, distMethod= "mann", metasnvOutDir, subpopOutDir,
                                        bamSuffix=".bam", rmdDir){

  rmdPath <- paste0(rmdDir, "/detailedSpeciesReport.rmd")

  tf <- paste0(subpopOutDir,"/tmp_",speciesID,"_",format(Sys.time(), "%Y-%m-%d-%H%M%S"))
  dir.create(tf)
  rmarkdown::render(rmdPath,
                    quiet = T,
                    knit_root_dir = getwd(),
                    params = list(
                      speciesID = speciesID,
                      subpopOutDir = subpopOutDir,
                      distMethod = distMethod,
                      metasnvOutDir = metasnvOutDir,
                      bamSuffix = bamSuffix),
                    intermediates_dir=tf, # so that reports can be rendered in parallel
                    output_dir = subpopOutDir,

                    output_file = paste0(speciesID,"_detailedSpeciesReport.html") )
  unlink(tf,recursive = T)
}


renderGeneContentReport <- function(speciesID, subpopOutDir, geneFamilyType,
                                    bamSuffix=".bam", rmdDir) {

  rStatCutoff <- 0.5
  statCutoff <- 0.05
  rmdPath <- paste0(rmdDir, "/geneContentReport.rmd")

  tf <- paste0(subpopOutDir,"/tmp_",speciesID,format(Sys.time(), "%Y-%m-%d-%H%M%S"))
  dir.create(tf)
  rmarkdown::render(rmdPath,quiet = T,
                    knit_root_dir = getwd(),
                    params = list(
                      speciesID = speciesID,
                      subpopOutDir = subpopOutDir,
                      geneFamilyType = geneFamilyType,
                      rStatCutoff = rStatCutoff,
                      statCutoff = statCutoff,
                      bamSuffix = bamSuffix),
                    intermediates_dir=tf, # so that reports can be rendered in parallel,
                    output_dir = subpopOutDir,
                    output_file = paste0(speciesID,"_geneContentReport.html"))
  unlink(x = tf,recursive = T)
}

#renderGeneContentReport(species = "357276")



# resultsDir <-"~/Dropbox/PostDocBork/subspecies/toolDevelopment/testingSubpopr/inSilicoMock/results_metaSNV-minFilter/"
# distMeth <- "mann"
#subpopOutDir="/Users/rossum/Dropbox/PostDocBork/subspecies/toolDevelopment/subpopr/ES_PC/metaSNV/fr11_v1/subspec/fromServer/"

renderTestPhenotypeAssocReport <- function(speciesID, subpopOutDir, metadataFile, categoryColumnNames,
                                           sampleIDColumnName,
                                           sampleExtension = ".bam",
                                           rmdDir) {

  tf <- paste0(subpopOutDir,"/tmp_",speciesID,format(Sys.time(), "%Y-%m-%d-%H%M%S"))
  dir.create(tf)
  #the original subspeceis defined from dist matrix clustering
  rmdPath <- paste0(rmdDir, "/testSubspecMultiPhenoAssoc.rmd")
  rmarkdown::render(rmdPath,quiet = T,
                    knit_root_dir = getwd(),
                    params = list(
                      speciesID = speciesID,
                      subpopOutDir = subpopOutDir,
                      metadataFile = metadataFile,
                      categoryColumnNames = categoryColumnNames,
                      sampleIDColumnName = sampleIDColumnName,
                      sampleExtension = sampleExtension),
                    intermediates_dir=tf, # so that reports can be rendered in parallel,
                    output_dir = subpopOutDir,
                    output_file = paste0(speciesID,paste0("_testSubspecMultiPhenoAssoc.html")) )
  unlink(tf,recursive = T)
}


renderResultsSummaryReport <- function(subpopOutDir,
                                       rmdDir){

  rmdPath <- paste0(rmdDir, "/resultsSummary.Rmd")
  rmarkdown::render(rmdPath,quiet = T,
                    knit_root_dir = getwd(),
                    params = list(
                      resultsDir = subpopOutDir),
                    output_dir = subpopOutDir,
                    output_file = "resultsSummary.html" )
}


# panCanCommResults_subpopr <- "/Users/rossum/Dropbox/PostDocBork/subspecies/pancreaticCancer/2018-06_data/subpopr/stoolSubspecMod/"
# panCanCommResults_metasnv <- "/Users/rossum/Dropbox/PostDocBork/subspecies/pancreaticCancer/2018-06_data/metaSNV/stool/outputs_minFilter/"
# panCanCommResults_motuProfileFilePath = "/Users/rossum/Dropbox/PostDocBork/subspecies/pancreaticCancer/2018-06_data/all_samples.motusv2.relabund.rmHeaderRenamed.tsv"
#
#
#  # panCanCommResults_subpoprMotus <- "/Users/rossum/Dropbox/PostDocBork/subspecies/pancreaticCancer/2018-06_data/subpopr_motus2/stool_NotQC-defaults/noClustering/"
#  # panCanCommResults_metasnvMotus <- "/Users/rossum/Dropbox/PostDocBork/subspecies/pancreaticCancer/2018-06_data/motus2/snv_inputNotQC/stool/metaSNVStyle/"
#  # panCanCommResults_motuProfileFilePath = "/Users/rossum/Dropbox/PostDocBork/subspecies/pancreaticCancer/2018-06_data/all_samples.motusv2.relabund.rmHeaderRenamed.tsv"
#
# # renderDetailedSpeciesReport("meta_mOTU_v2_5386",subpopOutDir=panCanCommResults_subpoprMotus, metasnvOutDir= panCanCommResults_metasnvMotus)
# renderDetailedSpeciesReport("469606",subpopOutDir=panCanCommResults_subpopr, metasnvOutDir= panCanCommResults_metasnvMotus)


# panCanCommResults_subpopr <- "/Users/rossum/Dropbox/PostDocBork/subspecies/pancreaticCancer/2018-06_data/v2UL/oral/subpopr/minFilter/"
# renderTestPhenotypeAssocReport(speciesID = "469606",subpopOutDir = panCanCommResults_subpopr,
#                                categoryColumnName = "status",
#                                sampleIDColumnName = "ID",
#                                sampleExtension = ".ULRepGenomesv11.unique.sorted.bam",
#                                metadataFile = "/Users/rossum/Dropbox/PostDocBork/subspecies/pancreaticCancer/2018-06_data/PC_totalmetadata_07.2018_pg1.csv",
#                                runOnGenotypedSamples = F)

#renderDetailedSpeciesReport("742765",subpopOutDir=panCanCommResults_subpopr, metasnvOutDir= panCanCommResults_metasnv)
# renderGeneContentReport(speciesID = 470146,
#                         subpopOutDir = panCanCommResults_subpopr,
#                         geneFamilyAbundancesFile = "/Users/rossum/Dropbox/PostDocBork/subspecies/pancreaticCancer/2018-06_data/mapped/igc/kegg.counts.all1.scaled.rmHeaderRenamed2UL.txt")

# mockCommResults_subpopr <- "/Users/rossum/Dropbox/PostDocBork/subspecies/toolDevelopment/testingSubpopr/inSilicoMock/subpopr_results/results_metaSNV-minFilter/substructure/"
# mockCommResults_metasnv <- "/Users/rossum/Dropbox/PostDocBork/subspecies/toolDevelopment/testingSubpopr/inSilicoMock/metaSNV_results/outputs_minFilter/"
# renderDetailedSpeciesReport(267377,subpopOutDir=mockCommResults_subpopr, metasnvOutDir= mockCommResults_metasnv)


#mockCommResults_subpopr <- "/Users/rossum/Dropbox/PostDocBork/subspecies/toolDevelopment/testingSubpopr/inSilicoMock/mine/subpopr/results/defaults/"
#mockCommResults_metasnv <- "/Users/rossum/Dropbox/PostDocBork/subspecies/toolDevelopment/testingSubpopr/inSilicoMock/mine/metaSNV/outputs_defaults/"
#renderDetailedSpeciesReport("referenceGenome",subpopOutDir=mockCommResults_subpopr, metasnvOutDir= mockCommResults_metasnv)



