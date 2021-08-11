# Changelog
All notable changes to this project will be documented in this file.


## Release - 2021-08-10

### Project changes:

- metaSNV code and maintenance moved to [GitHub](https://github.com/metasnv-tool/metaSNV).
- List of all the resources before the migration:
    - [EMBL GitLab](https://git.embl.de/costea/metaSNV)
    - [Website](http://metasnv.embl.de/)
    - [Initial paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0182392)
    
### Code changes:

**SNV-calling**
-	made compatible with Python 3
-   handle non-numeric species IDs
-   refined coverage counting to account for hard and soft clipping
-	support for ProGenomes2 as reference genomic database (aka "Freeze 11")
-	addition of "mOTUs mode" to support using mOTUs marker sequences as species references

**SNV post-processing**
-	added calculation of piN and piS (`metaSNV_distDiv.py` parameter `--divNS`)
-	added option to only compute SNV statistics on positions present in 90% of the samples (`metaSNV_distDiv.py` parameter `--matched`)
-	multi-threaded option for filtering & post processing
-	filtering output written to directory coded with filtering thresholds, facilitating multiple filtering runs
-   added subspecies calling module

## [v1.0.3] - 2018-04-25
Last release issued from the [EMBL GitLab](https://git.embl.de/costea/metaSNV).


