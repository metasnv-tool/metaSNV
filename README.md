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

```
wget http://swifter.embl.de/~ralves/metaSNV_reference_data/progenomes2_speciesReps_genomes.fna
wget http://swifter.embl.de/~ralves/metaSNV_reference_data/progenomes2_speciesReps_annotations.txt
```

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

This test example uses in silico generated data so that results compute quickly. There are 160 metagenomic samples with abundance of 3 species. The 'species' are called "refGenome1clus", "refGenome2clus", and "refGenome3clus", with 1, 2, and 3 subspecies, respectively. See the manual for details on the formats of the output files.

This will need 450M of disk space and, after the downloads are complete, should take approximately 5 minutes to run with 3 processors/threads or 15 minutes unparallelised.

1. Fetch and unpack test data

```
wget http://swifter.embl.de/~ralves/metaSNV_test_data/testdata.tar.xz
tar xvf testdata.tar.xz && rm -f testdata.tar.gz
```

**Run all steps of metaSNV v2 with the test data:**

2. Call SNVs:

```
python metaSNV.py --threads 3 output testdata/all_samples testdata/ref/allReferenceGenomes.fasta
```

Your SNVs are now in `output/snpCaller/called_SNPs`. You should have 6238 SNVs in this file, one per line.

3. Filter SNVs:

```
python metaSNV_Filtering.py --n_threads 3 output
```

This command filtered your SNVs (with default paramters) and calculated allele frequencies. Your filtered SNV allele frequencies are now in the `output/filtered/pop/` folder. Each species has its own file.

The number of SNVs expected per file are:
|# SNVs|File|
|-|-|
|1023|output/filtered/pop/refGenome1clus.filtered.freq|
|2075|output/filtered/pop/refGenome2clus.filtered.freq|
|3061|output/filtered/pop/refGenome3clus.filtered.freq|


4. Calculate distances between samples based on SNV profiles:

```
python metaSNV_DistDiv.py --n_threads 3 --filt output/filtered/pop --dist
```

This command calculated pairwise dissimilarities between samples based on filtered SNV allele frequencies. Your filtered SNV allele frequencies are now in the `output/distances` folder. Each species has its own file with 160 samples (161 lines with the header).

5. Detect clusters of samples that correspond to within-species subpopulations:

```
Rscript metaSNV_subpopr.R --procs 3 -i output -g testdata/abunds/geneAbundances.tsv -a testdata/abunds/speciesAbundances.tsv
```

This command detected subspecies. The results are in the `results/params.hr10.hs80.ps80.gs80/output/` folder. See the manual for a description of all the files produced. A summary of the results is provided in `results/params.hr10.hs80.ps80.gs80/output/resultsSummary.html` and `results/params.hr10.hs80.ps80.gs80/output/summary_allResults.csv` an overview of the results per species is provided in the files named with pattern: `refGenome*clus_detailedSpeciesReport.html`. 

This will include PCoA plots illustrating the clustering and tables with the number of distinctive genes. 
<img src="https://user-images.githubusercontent.com/6667809/126517961-d86d2c60-809c-4815-821a-8a76b1f4f3d5.png" alt="plot.png" width="50%"/>


See the full [documentation](https://github.com/metasnv-tool/metaSNV/tree/master/documentation) for details on parameters, inputs, outputs, and the method.
