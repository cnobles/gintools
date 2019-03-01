# gintools
Tools used for analysis and interpretation of foriegn DNA integration elements, such as retroviruses, retroviral-based gene therapy vectors, retrotransposons, or DNA oligos incorporated into a host genome by a Non-Homologous End Joining pathway.

Functions contained in the `gintools` package can be separated into different catagories:

* **Nucleotide sequence processing** - managing nucleotide sequence data from machine output, including demultiplexing, trimming, filtering, and consolidating. Alignment can be handled by the user's preferred sequence aligner (BLAT, BWA, ...).

* **Alignment interpretation** - interpret alignment information from PSL or SAM/BAM inputs. Identify important locations, standarize across samples, or resolve edges of observed nucleotide alignments to help reduce noise inherent in the data from PCR amplification and sequencing.

* **Analytics** - functions designed to assist in analysis of processed data. Track observations across specimens, determine abundances of clones, or consolidate data into a summary of observations. 

* **Utilities** - functions designed to assist with the other catagories and serve to make seemingly simple operations just that. 

## Install
To install `gintools`, simply run the following command within an R session:

```
devtools::install_github("https://github.com/cnobles/gintools.git")

# Or

devtools::install_github("cnobles/gintools")
```

## Functions by catagory

#### Nucleotide sequence processing
* **banmat** : A binary ambiguous nucleotide matrix based on NUC4.4.


#### Alignment interpretation
* **standardize_sites** : Returns a GRanges object where the site positions have been standardized with all other sites in the dataset which are within the gap distance.
* **refine_breakpoints** : Returns a GRanges object where the break point positions have been adjused based on positional clusterting and read counts within the dataset and specified distances.


#### Analytics
* **track_clones** : Returns a GRangesList of integration sites shared between multiple GRanges objects.
* **condense_intsites** : Returns a GRanges object containing a single integration site in each row, removing all breakpoint information.
* **determine_abundance** : Returns a data.frame with position ids and abundances, calculated by the number of unique fragment lengths or utilizing the sonicLength package.
* **cluster_kv** : Returns the cluster or group membership based on connections between key nodes based on supplied values.


#### Utilities
* **generate_posid** : Creates a character vector of position IDs in the format of [chromosome(+/-/*)position] given a GRanges object or ID information.
* **db_to_granges** : Converts an INSPIIRED database query into a GRanges object.
* **unique_granges** : Considers and keeps metadata columns when identifying unique ranges within a GRanges object.
* **pop_calcs** : Calculations for describing features of populations, i.e. Shannon Diversity, Clonality, Gini Index ...
* **vzip** : Combines two or more vectors in a "zipping" fashion and returns a single vector.
* **vintersect** : Identify intersecting values in two or more vectors.
* **vcollapse** : Collapse row contents into a single vector from a data.frame or matrix.

## Dependencies
The `gintools` package depends on `R 3.2` or higher and will `Import` the following packages during installation unless they are already present:

* dplyr (>= 0.7)
* GenomicRanges (>= 1.26)
* igraph (>= 1.0.1)
* IRanges (>= 2.10)
* magrittr(>= 1.0)
* Matrix (>= 1.2)
* parallel(>= 3.2)
* S4Vectors (>= 0.12)

The following packages are suggested for increased utility:

* Biostrings
* digest
* geneRxCluster
* plyr
* ShortRead
* stringr
* sonicLength
* testthat
