# Installing from source

## Source code

The latest source code can be downloaded via Git:

    git clone https://github.com/metasnv-tool/metaSNV.git
    
or [as a zip file](https://git.embl.de/metasnv-tool/metaSNV/repository/archive.zip?ref=master).

### Dependencies

**SNV calling**:

* [Boost-1.53.0 or above](http://www.boost.org/users/download/)
* [htslib](http://www.htslib.org/)
* Python-2.7 or above
    * numpy
    * pandas

**Subpopulation calling**:

* [R 3.5 or 3.6](https://www.r-project.org/)
* [Cairo](http://cairographics.org/)
* Several R libraries, see [here](src/subpopr/R/installOrLoadPackages.R)


#### Installing dependencies on Ubuntu/debian

On an Ubuntu/debian system, the following sequence of commands will install all
required packages (the first two are only necessary if you have not enabled the
universe repository before):

```
sudo add-apt-repository "deb http://archive.ubuntu.com/ubuntu $(lsb_release -sc) universe"
sudo apt-get update
sudo apt-get install libhts-dev libboost-dev
```

#### Installing dependencies using anaconda

If you use [anaconda](https://www.anaconda.com/products/individual), you can create an
environment with all necessary dependencies using the following commands:

```
conda create --name metaSNV -c bioconda boost htslib pkg-config numpy pandas
source activate metaSNV
export CONDA_ENV_PATH=$CONDA_PREFIX
export CFLAGS=-I$CONDA_ENV_PATH/include
export LD_LIBRARY_PATH=$CONDA_ENV_PATH/lib:$LD_LIBRARY_PATH
```

If you do not have a C++ compiler, anaconda can also install G++:

```
conda create --name metaSNV -c bioconda boost htslib pkg-config numpy pandas
source activate metaSNV
conda install gcc_linux-64 gxx_linux-64 
export CONDA_ENV_PATH=$CONDA_PREFIX
export CFLAGS=-I$CONDA_ENV_PATH/include
export LD_LIBRARY_PATH=$CONDA_ENV_PATH/lib:$LD_LIBRARY_PATH
```

For subpopulation calling, you will also need

```
conda install -c r r-essentials
conda install -c anaconda cairo
```

And several R libraries, see [here](src/subpopr/R/installOrLoadPackages.R) and call function installOrLoadPackagesfunction() from within R. e.g.

```
source("installOrLoadPackages.R")
installOrLoadPackagesfunction(doInstall=TRUE, doSuppressPackageStartupMessages=FALSE)
```

### Setup & Compilation

```
make
```
    
To test that all files and dependencies have been properly installed, run the following:

```
python metaSNV.py --help
python metaSNV_Filtering.py --help
python metaSNV_DistDiv.py --help
Rscript metaSNV_subpopr.R --help
```
    
**To use one of the provided reference databases**:

We recommend using genomes from ProGenomes2. The version provided here is a subset with one representative genome per species.

To download the species genome reference fasta file (`ref_db`) and the gene annotation file (`db_ann`), run the following:

```
wget http://swifter.embl.de/~ralves/metaSNV_reference_data/progenomes2_speciesReps_genomes.fna
wget http://swifter.embl.de/~ralves/metaSNV_reference_data/progenomes2_speciesReps_annotations.txt
```

These files will take approx. 25 GB of space.

