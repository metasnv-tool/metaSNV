# MetaSNV, a metagenomic SNV calling pipeline


The metaSNV pipeline performs variant calling on aligned metagenomic samples.


Download
========

Via Git:

    git clone https://github.com/metasnv-tool/metaSNV.git
    
or [download](https://git.embl.de/metasnv-tool/metaSNV/repository/archive.zip?ref=master) a zip file of the repository.

Dependencies
============

**SNV calling**:

* [Boost-1.53.0 or above](http://www.boost.org/users/download/)
* [htslib](http://www.htslib.org/)
* Python-2.7 or above
    * numpy
    * pandas

**Subpopulation calling**:

* R 3.5 or 3.6 [](https://www.r-project.org/)
* Cairo [](http://cairographics.org/)
* Several R libraries, see [here](src/subpopr/R/installOrLoadPackages.R)


#### Installing dependencies on Ubuntu/debian

On an Ubuntu/debian system, the following sequence of commands will install all
required packages (the first two are only necessary if you have not enabled the
universe repository before):


    sudo add-apt-repository "deb http://archive.ubuntu.com/ubuntu $(lsb_release -sc) universe"
    sudo apt-get update
    sudo apt-get install libhts-dev libboost-dev

#### Installing dependencies using anaconda

If you use [anaconda](https://www.anaconda.com/products/individual), you can create an
environment with all necessary dependencies using the following commands:

    conda create --name metaSNV -c bioconda boost htslib pkg-config numpy pandas
    source activate metaSNV
    export CONDA_ENV_PATH=$CONDA_PREFIX
    export CFLAGS=-I$CONDA_ENV_PATH/include
    export LD_LIBRARY_PATH=$CONDA_ENV_PATH/lib:$LD_LIBRARY_PATH

If you do not have a C++ compiler, anaconda can also install G++:

    conda create --name metaSNV -c bioconda boost htslib pkg-config numpy pandas
    source activate metaSNV
    conda install gcc_linux-64 gxx_linux-64 
    export CONDA_ENV_PATH=$CONDA_PREFIX
    export CFLAGS=-I$CONDA_ENV_PATH/include
    export LD_LIBRARY_PATH=$CONDA_ENV_PATH/lib:$LD_LIBRARY_PATH

For subpopulation calling, you will also need

    conda install -c r r-essentials
    conda install -c anaconda cairo

Setup & Compilation
===================

    make
    
To test that all files and dependencies have been properly installed, run the following:

    python metaSNV.py --help
    python metaSNV_Filtering.py --help
    python metaSNV_DistDiv.py --help
    Rscript metaSNV_subpopr.R --help

Workflow:
=========

### Part 0: Input files

* **'all\_samples'**  = a list of all BAM files, one /path/2/sample.bam per line (no duplicates)
* **'ref\_db'**       = the reference database in fasta format (f.i. multi-sequence fasta) and write permission for its directory
* **'db\_ann'**       = [optional] a gene annotation file for the reference database (format: ).

**To use one of the provided reference databases**:

The `ref_db` and `db_ann` files can be downloaded from using:

    ./getRefDB.sh

We recommend using ProGenomes2 (aka Freeze11).


### Part I: Call SNVs

    metaSNV.py project_dir/ all_samples ref_db [options]

### Part II: SNV Post-Processing (Filtering & Analysis)

Note: requires SNV calling (Part I) to be done

    metaSNV_Filtering.py project_dir [options]
    
    metaSNV_DistDiv.py --filt project_dir/filtered/pop [options]

### Part III: Subpopulation detection

Note: requires SNV calling, filtering, and distance calculations to be done (see 'Tutorial 2' below for example)

    metaSNV_subpopr.R -i project_dir [options]



Example Tutorial 1 (no subpopulation calling)
===================

## 1. Run the setup & compilation steps and download the provided reference database. 

    $ ./getRefDB.sh
    Select freeze9, as the tutorial files have been mapped against this freeze. 

## 2. Go to the EXAMPLE directory and download the samples with the getExp.sh

    $ cd EXAMPLE
    $ ./getExp.sh

## 3. Make sample list

    $ find `pwd`/EXAMPLE/samples -name "*.bam" > sample_list

## 4. Run the SNV calling step

    $ python metaSNV.py tutorial sample_list db/freeze9.genomes.RepGenomesv9.fna --threads 8

## 5. Run filtering and post processing

    $ python metaSNV_Filtering.py tutorial 
    $ python metaSNV_DistDiv.py --filt tutorial/filtered/pop --dist
    
    Voila! Your distances will be in the tutorial/distances folder. Enjoy!




Example Tutorial 2 (with subpopulation calling)
===================

## Fetch and unpack test data

    wget http://swifter.embl.de/~ralves/metaSNV_test_data/testdata.tar.xz
    tar xvf testdata.tar.xz && rm -f testdata.tar.gz

# Run all steps of metaSNV2 with the test data

Call SNVs:

    ./metaSNV.py output testdata/all_samples testdata/ref/allReferenceGenomes.fasta

Filter SNVs:

    ./metaSNV_Filtering.py output
    
Calculate distances between samples based on SNV profiles:
    
    ./metaSNV_DistDiv.py --filt output/filtered/pop --dist
    
Detect clusters of samples that correspond to within-species subpopulations:

    ./metaSNV_subpopr.R -i output -g testdata/abunds/geneAbundances.tsv -a testdata/abunds/speciesAbundances.tsv


Advanced usage 
==================================

If you want to run a lot of samples and would like to use the power of your cluster, we will print out the commands you need to
run and you can decide on how to schedule and manage them.

## 1. Get the first set of commands
    
    $ python metaSNV.py tutorial sample_list db/freeze9.genomes.RepGenomesv9.fna --n_splits 8 --print-commands
    
    Note the addition of the "--print-commands". This will print out one-liners that you need to run. When done, run same again.

## 2. Get the second set of commands
 
    $ python metaSNV.py tutorial sample_list db/freeze9.genomes.RepGenomesv9.fna --n_splits 8 --print-commands
    
    This will calculate the "load balancing" and give you the commands for running the SNV calling.
    
## 3. Run post-processing as usual

    $ python metaSNV_Filtering.py tutorial 
    $ python metaSNV_DistDiv.py --filt tutorial/filtered/pop --dist

