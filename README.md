# gintools
Tools used for analysis and interpretation of genomic DNA integration elements, such as retroviruses, retroviral-based gene therapy vectors, retrotransposons, and more. Built to work with GenomicRanges and the INSPIIRED pipeline. 

## List of functions:
**generate_posid** : Creates a character vector of position IDs in the format of [chromosome(+/-/*)position] given a GRanges object or ID information.

**db_to_granges** : Converts an INSPIIRED database query into a GRanges object.

**unique_granges** : Considers and keeps metadata columns when identifying unique ranges within a GRanges object.

**track_clones** : Returns a GRangesList of integration sites shared between multiple GRanges objects.

**condense_intsites** : Returns a GRanges object containing a single integration site in each row, removing all breakpoint information.

**determine_abundance** : Returns a data.frame with position ids and abundances, calculated by the number of unique fragment lengths or utilizing the sonicLength package.

**pop_calcs** : Calculations for describing features of populations, i.e. Shannon Diversity, Clonality, Gini Index ...

**standardize_sites** : Returns a GRanges object where the site positions have been standardized with all other sites in the dataset which are within the gap distance.

**refine_breakpoints** : Returns a GRanges object where the break point positions have been adjused based on positional clusterting and read counts within the dataset and specified distances.

**vzip** : Combines two or more vectors in a "zipping" fashion and returns a single vector.
