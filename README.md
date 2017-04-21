# MetaSNV, a metagenomic SNV calling pipeline


The metaSNV pipeline performs variant calling on aligned metagenomic samples.


Download
========

Via Git:

    git clone git@git.embl.de:costea/metaSNV.git
    
or [download](https://git.embl.de/rmuench/metaSNP/repository/archive.zip?ref=master) a zip file of the repository.

Dependencies
============

* [Boost-1.53.0 or above](http://www.boost.org/users/download/)

* [htslib](http://www.htslib.org/)
 
* Python-2.7 or above

#### Installing dependencies on Ubuntu/debian

On an Ubuntu/debian system, the following sequence of commands will install all
required packages (the first two are only necessary if you have not enabled the
universe repository before):


    sudo add-apt-repository "deb http://archive.ubuntu.com/ubuntu $(lsb_release -sc) universe"
    sudo apt-get update

    sudo apt-get install libhts-dev libboost-dev

### Installing dependencies using anaconda

If you use [anaconda](https://www.continuum.io/downloads), you can create an
environment with all necessary dependencies using the following commands:

    conda create --name metaSNV boost htslib pkg-config
    source activate metaSNV
    export CFLAGS=-I$CONDA_ENV_PATH/include
    export LD_LIBRARY_PATH=$CONDA_ENV_PATH/lib:$LD_LIBRARY_PATH

If you do not have a C++ compiler, anaconda can also install G++:

    conda create --name metaSNV boost htslib pkg-config
    source activate metaSNV

    # Add this command:
    conda install gcc

    export CFLAGS=-I$CONDA_ENV_PATH/include
    export LD_LIBRARY_PATH=$CONDA_ENV_PATH/lib:$LD_LIBRARY_PATH

Setup & Compilation
===================

    make

Workflow:
=========
## Required Files:

* **'all\_samples'**    = a list of all BAM files, one /path/2/sample.bam per line (no duplicates)
* **'ref\_db'**  = the reference database in fasta format (f.i. multi-sequence fasta)

## Optional Files:
* **'db\_ann'** = a gene annotation file for the reference database.

## To use one of the provided reference databases:

    ./getRefDB.sh
    
## 2. Run metaSNV

    metaSNV.py project_dir/ all_samples ref_db ref_fasta [options]

## 3. Part II: Post-Processing (Filtering & Analysis)

### a) Filtering:
Note: requires SNP calling (Part II) to be done!
Caution: Perform this step seperately for individual SNPs and population SNPs.

    usage: metaSNP_filtering.py

        positional arguments:
            perc_FILE               input file with horizontal genome (taxon) coverage (breadth) per sample (percentage covered)
            cov_FILE                input file with average genome (taxon) coverage
                                    (depth) per sample (average number reads per site)
            snp_FILE                input files from SNP calling
            all_samples             list of input BAM files, one per line
            output_dir/             output folder

        optional arguments:
            -h, --help              show this help message and exit
            -p PERC, --perc PERC    Coverage breadth: Horizontal coverage cutoff
                                    (percentage taxon covered) per sample (default: 40.0)
            -c COV, --cov COV       Coverage depth: Average vertical coverage cutoff per
                                    taxon, per sample (default: 5.0)
            -m MINSAMPLES, --minsamples MINSAMPLES
                                    Minimum number of sample required to pass the coverage
                                    cutoffs per genome (default: 2)
            -s SNPC, --snpc SNPC    FILTERING STEP II: SNP coverage cutoff (default: 5.0)
            -i SNPI, --snpi SNPI    FILTERING STEP II: SNP occurence (incidence) cutoff
                                    within samples_of_interest (default: 0.5)

### b) Analysis:
>   TODO: include R scripts for computing pairwise distances and visualization


Example Tutorial
================

## 1. Run the setup & compilation steps and download the provided reference database.

    ./getRefDB.sh

## 2. Go to the EXAMPLE directory and download the samples with the getSamplesScript.sh

    $ cd EXAMPLE
    $ ./getSamplesScript.sh

## 3. Initiate a new project in the parent directory

    $ metaSNV_New tutorial

## 4. Generate the 'all_samples' file

    $ find `pwd`/EXAMPLE/samples -name “*.bam” > tutorial/all_samples

## 5. Prepare and run the coverage estimation

    $ metaSNV_COV tutorial/ tutorial/all_samples > runCoverage
    $ bash runCoverage

## 6. Perform a work load balancing step for run time optimization.

    $ metaSNV_OPT tutorial/ db/Genomev9_definitions 5
    $ bash runCoverage

## 7. Prepare and run the SNV calling step

    $ metaSNV_SNP tutorial/ tutorial/all_samples db/RepGenomesv9.fna -a db/RefOrganismDB_v9_gene.clean -l tutorial/bestsplits/ > runSNPcall
    $ bash runSNPcall

## 8. Run the post processing / filtering steps
### a) Compute allele frequencies for each position that pass the given thresholds.

    $ metaSNV_filtering.py tutorial/tutorial.all_perc.tab tutorial/tutorial.all_cov.tab tutorial/snpCaller/called_SNVs.best_split_* tutorial/all_samples tutorial/filtered/pop/

### b) Compute pair-wise distances between samples on their SNP profiles and create a PCoA plot.
    



Advanced usage (tools and scripts)
==================================

If you are interested in using the pipeline in a more manual way (for example
the metaSNV caller stand alone), you will find the executables for the
individual steps in the `src/` directory.

metaSNV caller
--------------
Calls SNVs from samtools pileup format and generates two outputs.

    usage: ./snpCall [options] < stdin.mpileup > std.out.popSNPs

    Options:
        -f,     faidx indexed reference genome.
        -g,     gene annotation file.
        -i,     individual SNVs.

Note: Expecting samtools mpileup as standard input

### __Output__
1. Population SNVs (pSNVs): 
Population wide variants that occur with a frequency of 1 % at positions with at least 4x coverage.

2. Individual specific SNVs (iSNVs):
Non population variants, that occur with a frequency of 10 % at positions with at least 10x coverage.


[qaCompute](https://github.com/CosteaPaul/qaTools)
-------------------------------------------------
Computes normal and span coverage from a bam/sam file. 
Also counts unmapped and sub-par quality reads.

### __Parameters:__
   m	    -	Compute median coverage for each contig/chromosome. 
   		Will make running a bit slower. Off by default.
   
   q [INT]  -   Quality threshold. Any read with a mapping quality under
                INT will be ignored when computing the coverage.
                
		NOTE: bwa outputs mapping quality 0 for reads that map with
		equal quality in multiple places. If you want to condier this,
		set q to 0.
		
   d        -   Print coverage histrogram over each individual contig/chromosome.
   	        These details will be printed in file <output>.detail
   	        
   p [INT]  -   Print coverage profile to bed file, averaged over given window size.  
   
   i        -   Silent run. Will not print running info to stdout.
    
   s [INT]  -   Compute span coverage. (Use for mate pair libs)
                Instead of actual read coverage, using the options will consider
                the entire span of the insert as a read, if insert size is
		lower than INT. 
 		For an accurate estimation of span coverage, I recommend
		setting an insert size limit INT around 3*std_dev of your lib's 
		insert size distribution.
    
   c [INT]  -   Maximum X coverage to consider in histogram.
    
   h [STR]  -   Use different header. 
                Because mappers sometimes break the headers or simply don't output them, 
		this is provieded as a non-kosher way around it. Use with care!
    
   For more info on the parameteres try ./qaCompute
   

metaSNV_filtering.py
--------------------   
usage: metaSNV filtering [-h] [-p PERC] [-c COV] [-m MINSAMPLES] [-s SNPC]
                         [-i SNPI] 
                         perc_FILE cov_FILE snp_FILE [snp_FILE ...]
                         all_samples output_dir/

metaSNV filtering

positional arguments:
  perc_FILE             input file with horizontal genome (taxon) coverage
                        (breadth) per sample (percentage covered)
  cov_FILE              input file with average genome (taxon) coverage
                        (depth) per sample (average number reads per site)
  snp_FILE              input files from SNP calling
  all_samples           list of input BAM files, one per line
  output_dir/           output folder

optional arguments:
  -h, --help            show this help message and exit
  -p PERC, --perc PERC  Coverage breadth: Horizontal coverage cutoff
                        (percentage taxon covered) per sample (default: 40.0)
  -c COV, --cov COV     Coverage depth: Average vertical coverage cutoff per
                        taxon, per sample (default: 5.0)
  -m MINSAMPLES, --minsamples MINSAMPLES
                        Minimum number of sample that have to pass the
                        filtering criteria in order to write an output for the
                        representative Genome (default: 2)
  -s SNPC, --snpc SNPC  FILTERING STEP II: SNP coverage cutoff (default: 5.0)
  -i SNPI, --snpi SNPI  FILTERING STEP II: SNP occurence (incidence) cutoff
                        within samples_of_interest (default: 0.5)


