# if loading Cairo fails, can also use conda (e.g. conda install -c anaconda cairo)
if(!requireNamespace("BiocManager", quietly=TRUE)){install.packages("BiocManager", repos="https://cloud.r-project.org/")}
if(!requireNamespace("BiocParallel", quietly=TRUE)){BiocManager::install("BiocParallel"); library(BiocParallel)}
if(!requireNamespace("coin", quietly=TRUE)){install.packages("coin", repos="https://cloud.r-project.org/"); library(coin)}
if(!requireNamespace("fpc", quietly=TRUE)){install.packages("fpc", repos="https://cloud.r-project.org/"); library(fpc)}
if(!requireNamespace("ape", quietly=TRUE)){install.packages("ape", repos="https://cloud.r-project.org/"); library(ape)}
if(!requireNamespace("getopt", quietly=TRUE)){install.packages("getopt", repos="https://cloud.r-project.org/"); library(getopt)}
if(!requireNamespace("optparse", quietly=TRUE)){install.packages("optparse", repos="https://cloud.r-project.org/"); library(optparse)}
if(!requireNamespace("readr", quietly=TRUE)){install.packages("readr", repos="https://cloud.r-project.org/"); library(readr)}
if(!requireNamespace("dplyr", quietly=TRUE)){install.packages("dplyr", repos="https://cloud.r-project.org/"); library(dplyr)}
if(!requireNamespace("ggplot2", quietly=TRUE)){install.packages("ggplot2", repos="https://cloud.r-project.org/"); library(ggplot2)}
if(!requireNamespace("tidyr", quietly=TRUE)){install.packages("tidyr", repos="https://cloud.r-project.org/"); library(tidyr)}
if(!requireNamespace("gridExtra", quietly=TRUE)){install.packages("gridExtra", repos="https://cloud.r-project.org/"); library(gridExtra)}
if(!requireNamespace("lemon", quietly=TRUE)){install.packages("lemon", repos="https://cloud.r-project.org/"); library(lemon)}
if(!requireNamespace("DT", quietly=TRUE)){install.packages("DT", repos="https://cloud.r-project.org/"); library(DT)}
if(!requireNamespace("cluster", quietly=TRUE)){install.packages("cluster", repos="https://cloud.r-project.org/"); library(cluster)}
if(!requireNamespace("ggrepel", quietly=TRUE)){install.packages("ggrepel", repos="https://cloud.r-project.org/"); library(ggrepel)}
if(!requireNamespace("data.table", quietly=TRUE)){install.packages("data.table", repos="https://cloud.r-project.org/"); library(data.table)}
if(!requireNamespace("kableExtra", quietly=TRUE)){install.packages("kableExtra", repos="https://cloud.r-project.org/"); library(kableExtra)}
if(!requireNamespace("rmarkdown", quietly=TRUE)){install.packages("rmarkdown", repos="https://cloud.r-project.org/"); library(rmarkdown)}
if(!requireNamespace("batchtools")){install.packages("batchtools", repos="https://cloud.r-project.org/"); library(batchtools)}
if(!requireNamespace("questionr", quietly=TRUE)){install.packages("questionr", repos="https://cloud.r-project.org/"); library(questionr)}
if(!requireNamespace("futile.logger", quietly=TRUE)){install.packages("futile.logger", repos="https://cloud.r-project.org/"); library(futile.logger)}
if(!requireNamespace("Cairo", quietly=TRUE)){install.packages("Cairo", repos="https://cloud.r-project.org/"); library(Cairo)}
if(!requireNamespace("reshape2", quietly=TRUE)){install.packages("reshape2", repos="https://cloud.r-project.org/"); library(reshape2)}
