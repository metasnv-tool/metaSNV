# MetaSNV, a metagenomic SNV calling pipeline


The metaSNV pipeline performs variant calling on aligned metagenomic samples.


Download
========

Via Git:

    git clone git@git.embl.de:costea/metaSNV.git
    
or [download](https://git.embl.de/costea/metaSNV/repository/archive.zip?ref=master) a zip file of the repository.

Dependencies
============

* [Boost-1.53.0 or above](http://www.boost.org/users/download/)

* [htslib](http://www.htslib.org/)
 
* Python-2.7 or above
    * numpy
    * pandas

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

    conda create --name metaSNV boost htslib pkg-config numpy pandas
    source activate metaSNV
    export CFLAGS=-I$CONDA_ENV_PATH/include
    export LD_LIBRARY_PATH=$CONDA_ENV_PATH/lib:$LD_LIBRARY_PATH

If you do not have a C++ compiler, anaconda can also install G++:

    conda create --name metaSNV boost htslib pkg-config numpy pandas
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

* **'all\_samples'**  = a list of all BAM files, one /path/2/sample.bam per line (no duplicates)
* **'ref\_db'**       = the reference database in fasta format (f.i. multi-sequence fasta)
* **'gen\_pos'**      = a list with start and end positions for each sequence in the reference (format: `sequence\_id  start end`)

## Optional Files:
* **'db\_ann'** = a gene annotation file for the reference database (format: ).

## To use one of the provided reference databases:

    ./getRefDB.sh
    
## 2. Run metaSNV

    metaSNV.py project_dir/ all_samples ref_db [options]

## 3. Part II: Post-Processing (Filtering & Analysis)
Note: requires SNP calling (Part II) to be done!

    metaSNV_post.py project_dir [options]

Example Tutorial
================

## 1. Run the setup & compilation steps and download the provided reference database. 

    $ ./getRefDB.sh
    Select freeze9, as the tutorial files have been mapped against this freeze. 

## 2. Go to the EXAMPLE directory and download the samples with the getSamplesScript.sh

    $ cd EXAMPLE
    $ ./getSamplesScript.sh

## 3. Make sample list

    $ find `pwd`/EXAMPLE/samples -name “*.bam” > sample_list

## 4. Run the SNV calling step

    $ python metaSNV.py tutorial sample_list db/freeze9.genomes.RepGenomesv9.fna --threads 8

## 5. Run filtering and post processing

    $ python metaSNV_post.py tutorial
    
    Voila! Your distances will be in the tutorial/distances folder. Enjoy!

Advanced usage 
==================================

If you want to run a lot of samples and would like to use the power of your cluster, we will print out the commands you need to
run and you can decide on how to schedule and manage them.

## 1. Get the first set of commands
    
    $ python metaSNV.py tutorial sample_list db/freeze9.genomes.RepGenomesv9.fna --n_splits 8 --print-commands
    
    Note the addition of the "--print-commnads". This will print out one-liners that you need to run. When done, run same again.

## 2. Get the second set of commands
 
    $ python metaSNV.py tutorial sample_list db/freeze9.genomes.RepGenomesv9.fna --n_splits 8 --print-commands
    
    This will calculate the "load balancing" and give you the commands for running the SNV calling.
    
## 3. Run post-processing as usual

    $ python metaSNV_post.py tutorial

