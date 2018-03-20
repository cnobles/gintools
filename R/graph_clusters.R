#' Generate a graph from integrations sites within a specified genomic gap
#' distance
#'
#' \code{graph_clusters} generates a partial undirected graph connecting each
#' integration sites within a specified genomic gap distance.
#'
#' @description Given a set of integration sites, this function generates a
#' partial undirected graph of integration sites where nodes are the specific
#' sites (one row from the GRanges object input) and edges mark which sites are
#' within the gap distance from eachother. To generate the partial graph, nodes
#' are broken into two groups, axil nodes and nodes. Axil nodes are the first
#' occurances (or left most occurance) of the unique integration sites
#' (chromsome + strand + position) which cluster together. Nodes are all other
#' occurances, which are always connected to one axil node by definition.
#'
#' @usage
#' graph_clusters(sites, gap)
#'
#' @param sites a GRanges object where each row is a single integration site or
#' range.
#'
#' @param gap an integer specifying the distance to consider between sites to
#' call an edge.
#'
#' @examples
#' sites <- GRanges(
#'   seqnames = rep("chr1", 7),
#'   ranges = IRanges(
#'     start = c(50, 49, 44, 50, 50, 50, 60),
#'     width = seq(20, 26, 1)),
#'   strand = rep("+", 7))
#'
#' graph <- graph_clusters(sites, gap = 5L)
#'
#' plot.igraph(graph)
#'
#' @author Christopher Nobles, Ph.D.
#'

graph_clusters <- function(sites, gap){
  sites$clus.key <- 1:length(sites)
  fl.sites <- GenomicRanges::flank(sites, width = -1, start = TRUE)
  rd.sites <- GenomicRanges::reduce(
    fl.sites, min.gapwidth = gap, with.revmap = TRUE)
  revmap <- rd.sites$revmap
  axil_nodes <- as.numeric(S4Vectors::Rle(
    values = sites[sapply(revmap, "[[", 1)]$clus.key,
    lengths = sapply(revmap, length)
  ))
  nodes <- sites[unlist(IRanges::as.list(revmap))]$clus.key
  edgelist <- unique(matrix( c(axil_nodes, nodes), ncol = 2 ))
  igraph::graph.edgelist(edgelist, directed = FALSE)
}
