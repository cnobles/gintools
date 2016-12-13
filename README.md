# gintools
Tools used for analysis and interpretation of genomic DNA integration elements, such as retroviruses, retroviral-based gene therapy vectors, retrotransposons, and more. Built to work with GenomicRanges and the INSPIIRED pipeline. 

## List of functions:
generate_posid() : Creates a character vector of position IDs in the format of [chromosome(+/-/*)position] given a GRanges object or ID information.

db_to_granges() : Converts an INSPIIRED database query into a GRanges object.

graph_clusters() : Generates a partial undirected graph connecting integration sites within a specified genomic gap distance.

serial_cluster() : Returns a GRanges object of the same length and order as input with additional metadata columns specifying the group or clusterID of given nt windows.

track_clones() : Returns a GRangesList of integration sites shared between multiple GRanges objects.

condense_intsites() : Returns a GRanges object containing a single integration site in each row, removing all breakpoint information.

determine_abundance() : Returns a data.frame with position ids and abundances, calculated by the number of unique fragment lengths or utilizing the sonicLength package.

normalize_multihit_clusters() : Annotates multihit GRanges objects with new multihit.ID's which can be used to identify the same multihit in another sample.

population_calcs() : Calculations for describing features of populations, i.e. Shannon Diversity, Clonality, Gini Index ...

scan_format() : Format a set of GRanges objects for geneRxCluster scan statistics package.

standardize_intsites() : Returns a GRanges object where the integration site positions have been standardized with all other sites in the dataset which are within the gap distance.

refine_breakpoints() : Returns a GRanges object where the break point positions have been adjused based on positional clusterting and read counts within the dataset and specified distances.

connect_satalite_verticies() : Used by standardize_intsites() and refine_breakpoints(), this function returns a graph based on input GRanges positions where nodes within a 'gap' distance from clusters are now connected to the boundries of the clusters.

break_connecting_source_paths() : Used by standadize_intsites() and refine_breakpoints(), this function returns a graph where only one source is present per cluster.

connect_adjacent_clusters() : Used by standardize_intsites(), this function returns a graph where adjacent clusters having sources within the gap distance are joined.

resolve_cluster_sources() : Used by standardize_intsites(), this function returns a graph where each cluster only has a single primary source node, but does not reduce the size of the cluster.

sinks() and sources() : Given a directed graph, returns a numerical vector of sink / source nodes.
