# MetaSNV, a metagenomic SNV calling pipeline


metaSNV supports variant (SNV) calling on metagenomic data and population genetic analysis, including subspecies identification and profiling. Input is metagenomic reads mapped against reference genomes (bam files).

- metaSNV v2 paper: coming soon
- [metaSNV v1 paper](https://doi.org/10.1371/journal.pone.0182392): Costea, et al. 2017. _metaSNV: A tool for metagenomic strain level analysis._ PLOS One.


See the full [documentation](https://github.com/metasnv-tool/metaSNV/tree/master/documentation) for more details and an test example on this page below.

## Installing

### Installing via conda

Create a new conda environment with metaSNV installed:

```
conda create --name metaSNV -c bioconda -c conda-forge 'metasnv>=2.0.1'
```

Install metaSNV in an existing conda environment:

```
conda install -c bioconda -c conda-forge 'metasnv>=2.0.1'
```

To test that all files and dependencies have been properly installed, run the following:

```
metaSNV.py --help
metaSNV_Filtering.py --help
metaSNV_DistDiv.py --help
metaSNV_subpopr.R --help
```

### Installing from source

Refer to the [developer guide](DEVELOPER.md).

## Preparing to analyse your data

We recommend using genomes from ProGenomes2 as your species references. The version provided here is a subset with one representative genome per species (the longest of the represenatatives).

To use the provided reference database, download the species genome reference fasta file (`ref_db`) and the gene annotation file (`db_ann`) with the following. 
These files will take approx. 25 GB of space.

```
wget http://swifter.embl.de/~ralves/metaSNV_reference_data/progenomes2_speciesReps_genomes.fna
wget http://swifter.embl.de/~ralves/metaSNV_reference_data/progenomes2_speciesReps_annotations.txt
```

## Workflow

This summarises the general workflow. See the full documentation for details.

### Part 0: Input files

* **'all\_samples'**  = a list of all BAM files, one /path/to/sample.bam per line (no duplicates, no empty mappings)
* **'ref\_db'**       = the reference database in fasta format (i.e. multi-sequence fasta) and write permission for its directory
* **'db\_ann'**       = [optional] a gene annotation file for the reference database (see documentation for custom files).

### Part I: Call SNVs

```
metaSNV.py output_dir/ all_samples ref_db [options]
```

### Part II: SNV Post-Processing: Filtering & Analysis

Note: requires SNV calling (Part I) to be done

```
metaSNV_Filtering.py output_dir [options]

metaSNV_DistDiv.py --filt output_dir/filtered/pop [options]
```

### Part III: Subpopulation detection

Note: requires SNV calling, filtering, and distance calculations to be done (Parts I & II)

```
metaSNV_subpopr.R -i output_dir [options]
```

To determine abundances of subspecies relative to the whole community, you will also need to provide species abundance profiles. To determine subspecies-associated gene content, you will need to provide per-metagenome gene abundance profiles. See the full documentation for details.

## Example Tutorial

This test example uses in silico generated data so that results compute quickly. There are 160 metagenomic samples with abundance of 3 species. The 'species' are called "refGenome1clus", "refGenome2clus", and "refGenome3clus", with 1, 2, and 3 subspecies, respectively. See the manual for details on the formats of the output files.

This will need 450M of disk space and, after the downloads are complete, should take approximately 5 minutes to run with 3 processors/threads or 15 minutes unparallelised.

### 1. Fetch and unpack test data

```
wget http://swifter.embl.de/~ralves/metaSNV_test_data/testdata.tar.xz
tar xvf testdata.tar.xz && rm -f testdata.tar.xz
```

**Run all steps of metaSNV v2 with the test data:**

### 2. Call SNVs:

```
metaSNV.py --threads 3 output testdata/all_samples testdata/ref/allReferenceGenomes.fasta
```

or if running from source

```
python metaSNV.py --threads 3 output testdata/all_samples testdata/ref/allReferenceGenomes.fasta
```

To also detect whether SNVs result in codon changes, add `--db_ann testdata/ref/metaSNV_anntotationsAll.txt` which will add to the SNV output: the gene within which the SNV was detected and the original and resultant codon. The file `testdata/ref/metaSNV_anntotationsAll.txt` contains gene names and locations within the reference genomes. See the documentation for details.

Your SNVs are now in `output/snpCaller/`.

- If you ran with >1 thread: Your SNVs are now in files named with the pattern: `output/snpCaller/called_SNPs.best_split_*`. If you ran with 3 threads, then each species will have it's own file and you should have 1080, 2094, and 3064 SNVs per file, one per line.
- If you ran with 1 thread: Your SNVs are now in `output/snpCaller/called_SNPs`. You should have 6238 SNVs in this file, one per line.

### 3. Filter SNVs:

```
metaSNV_Filtering.py --n_threads 3 output
```

or if running from source

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


### 4. Calculate distances between samples based on SNV profiles:

```
metaSNV_DistDiv.py --n_threads 3 --filt output/filtered/pop --dist
```

or if running from source

```
python metaSNV_DistDiv.py --n_threads 3 --filt output/filtered/pop --dist
```

This command calculated pairwise dissimilarities between samples based on filtered SNV allele frequencies. Your filtered SNV allele frequencies are now in the `output/distances` folder. Each species has its own file with 160 samples (161 lines with the header).

### 5. Detect clusters of samples that correspond to within-species subpopulations:

```
metaSNV_subpopr.R --procs 3 -i output -g testdata/abunds/geneAbundances.tsv -a testdata/abunds/speciesAbundances.tsv
```

or if running from source

```
Rscript metaSNV_subpopr.R --procs 3 -i output -g testdata/abunds/geneAbundances.tsv -a testdata/abunds/speciesAbundances.tsv
```

This command detected subspecies. The results are in the `results/params.hr10.hs80.ps80.gs80/output/` folder. See the manual for a description of all the files produced.

A summary of the results is provided in `results/params.hr10.hs80.ps80.gs80/output/resultsSummary.html` and `results/params.hr10.hs80.ps80.gs80/output/summary_allResults.csv`.

![image](https://user-images.githubusercontent.com/6667809/126518980-b8521bbd-0cb8-4433-b16d-1abb4e5bf1ed.png)


An overview of the results per species is provided in the files named with pattern: `refGenome*clus_detailedSpeciesReport.html`.

This will include PCoA plots illustrating the clustering and tables with the number of distinctive genes. For example, for in silico "species" "refGenome3clus":

<img src="https://user-images.githubusercontent.com/6667809/126517961-d86d2c60-809c-4815-821a-8a76b1f4f3d5.png" alt="plot.png" width="50%"/>

See the full [documentation](https://github.com/metasnv-tool/metaSNV/tree/master/documentation) for details on parameters, inputs, outputs, and the method.
