# MetaSNV, a metagenomic SNV calling pipeline


metaSNV supports variant (SNV) calling on metagenomic data and population genetic analysis, including subspecies identification and profiling. Input is metagenomic reads mapped against reference genomes (bam files).

- metaSNV v2 paper: coming soon
- metaSNV v1 paper: https://doi.org/10.1371/journal.pone.0182392


See the full documentation for more details and a tutorial on this page.

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

* **'all\_samples'**  = a list of all BAM files, one /path/to/sample.bam per line (no duplicates, no empty mappings)
* **'ref\_db'**       = the reference database in fasta format (i.e. multi-sequence fasta) and write permission for its directory
* **'db\_ann'**       = [optional] a gene annotation file for the reference database (see documentation for custom files).

**To use one of the provided reference databases**:

The `ref_db` and `db_ann` files can be downloaded from using:

    ./getRefDB.sh

We recommend using ProGenomes2 (aka Freeze11). The version provided here has one representative genome per species.


### Part I: Call SNVs

    metaSNV.py output_dir/ all_samples ref_db [options]

### Part II: SNV Post-Processing (Filtering & Analysis)

Note: requires SNV calling (Part I) to be done

    metaSNV_Filtering.py output_dir [options]
    
    metaSNV_DistDiv.py --filt output_dir/filtered/pop [options]

### Part III: Subpopulation detection

Note: requires SNV calling, filtering, and distance calculations to be done (Parts I & II)

    metaSNV_subpopr.R -i output_dir [options]



Example Tutorial 1 (no subspecies identification)
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




Example Tutorial 2 (with subspecies identification)
===================

## Fetch and unpack test data

    wget http://swifter.embl.de/~ralves/metaSNV_test_data/testdata.tar.xz
    tar xvf testdata.tar.xz && rm -f testdata.tar.gz

## Run all steps of metaSNV v2 with the test data

Call SNVs:

    ./metaSNV.py output testdata/all_samples testdata/ref/allReferenceGenomes.fasta

Filter SNVs:

    ./metaSNV_Filtering.py output
    
Calculate distances between samples based on SNV profiles:
    
    ./metaSNV_DistDiv.py --filt output/filtered/pop --dist
    
Detect clusters of samples that correspond to within-species subpopulations:

    ./metaSNV_subpopr.R -i output -g testdata/abunds/geneAbundances.tsv -a testdata/abunds/speciesAbundances.tsv


See the full documentation for details on parameters, inputs, outputs, and the method.
