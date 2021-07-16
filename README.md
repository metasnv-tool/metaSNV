# MetaSNV, a metagenomic SNV calling pipeline


metaSNV supports variant (SNV) calling on metagenomic data and population genetic analysis, including subspecies identification and profiling. Input is metagenomic reads mapped against reference genomes (bam files).

- metaSNV v2 paper: coming soon
- metaSNV v1 paper: https://doi.org/10.1371/journal.pone.0182392


See the full [documentation](https://github.com/metasnv-tool/metaSNV/tree/master/documentation) for more details and an test example on this page below.

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
    
**To use one of the provided reference databases**:

We recommend using genomes from ProGenomes2. The version provided here is a subset with one representative genome per species.

To download the species genome reference fasta file (`ref_db`) and the gene annotation file (`db_ann`), run the following:

***** TO ADD ********

These files will take approx. 25 GB of space.

Workflow:
=========

This summarises the general workflow. See the full documentation for details.

### Part 0: Input files

* **'all\_samples'**  = a list of all BAM files, one /path/to/sample.bam per line (no duplicates, no empty mappings)
* **'ref\_db'**       = the reference database in fasta format (i.e. multi-sequence fasta) and write permission for its directory
* **'db\_ann'**       = [optional] a gene annotation file for the reference database (see documentation for custom files).

### Part I: Call SNVs

    metaSNV.py output_dir/ all_samples ref_db [options]

### Part II: SNV Post-Processing: Filtering & Analysis

Note: requires SNV calling (Part I) to be done

    metaSNV_Filtering.py output_dir [options]
    
    metaSNV_DistDiv.py --filt output_dir/filtered/pop [options]

### Part III: Subpopulation detection

Note: requires SNV calling, filtering, and distance calculations to be done (Parts I & II)

    metaSNV_subpopr.R -i output_dir [options]

To determine abundances of subspecies relative to the whole community, you will also need to provide species abundance profiles. To determine subspecies-associated gene content, you will need to provide per-metagenome gene abundance profiles. See the full documentation for details.

Example Tutorial 
===================

This test example uses in silico generated data and will take more space and time to complete than the previous one due to the larger number of samples required for subspecies identification.

1. Fetch and unpack test data

```
wget http://swifter.embl.de/~ralves/metaSNV_test_data/testdata.tar.xz
tar xvf testdata.tar.xz && rm -f testdata.tar.gz
```

This will need 400M of space and should take less than 10 minutes to run in total.

**Run all steps of metaSNV v2 with the test data:**

2. Call SNVs:

```
python metaSNV.py output testdata/all_samples testdata/ref/allReferenceGenomes.fasta
```

3. Filter SNVs:

```
python metaSNV_Filtering.py output
```

4. Calculate distances between samples based on SNV profiles:

```
python metaSNV_DistDiv.py --filt output/filtered/pop --dist
```

5. Detect clusters of samples that correspond to within-species subpopulations:

```
Rscript metaSNV_subpopr.R -i output -g testdata/abunds/geneAbundances.tsv -a testdata/abunds/speciesAbundances.tsv
```

See the full [documentation](https://github.com/metasnv-tool/metaSNV/tree/master/documentation) for details on parameters, inputs, outputs, and the method.
