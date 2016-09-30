# gintools
Tools used for analysis and interpretation of genomic DNA integration elements, such as retroviruses, retroviral-based gene therapy vectors, retrotransposons, and more. Built to work with GenomicRanges and the INSPIIRED pipeline. 

## List of functions:
generate_posid() : Creates a character vector of position IDs in the format of [chromosome(+/-/*)position] given a GRanges object or ID information.

db_to_granges() : Converts an INSPIIRED database query into a GRanges object.

graph_clusters() : Generates a partial undirected graph connecting integration sites within a specified genomic gap distance.

serial_cluster() : Returns a GRanges object of the same length and order as input with additional metadata columns specifying the group or clusterID of given nt windows.
