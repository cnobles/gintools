#' Identified clusters for multihit integration sites.
#'
#' \code{cluster_multihits} returns a GRanges object with cluster information.
#'
#' @description Given a GRanges object of integration sites which align to
#' multiple location on the reference genome (multihit sites), this function
#' will append a column with cluster designation. Multihits are considered to be
#' part of the same cluster if they share any alignments with other multihits.
#' This function tackles this memory intensive process through iterive measures,
#' involving randomly sampling from the given set of data to start cluster
#' formation, refinement, and then resamples from the remaining data till all
#' information has been used.
#'
#' @usage
#' cluster_multihits(multihits, read_col = NULL)
#'
#' cluster_multihits(multihits, read_col = NULL, max_gap = 5L, iterations = 5L)
#'
#' @param multihits a GRanges object containing sets of integration sites
#' or ranges (one alignment per row) and containing several columns, inlcluding
#' "ID" (read identifier) and "key_pair" (unique identifier R2 and R1
#' sequences).
#'
#' @param read_col character string matching the name of the column of the
#' GRanges object given in 'multihits' which designates which ranges are
#' associated together. For example, the read name can be used here to show
#' which alignments were associated with the same read.
#'
#' @param max_gap an integer designating the nucleotide distance or window for 
#' which to group integration sites.
#'
#' @param iterations integer The number of interations of cluster generation to
#' perform and the number of random subsets that will be made from the data.
#' More iterations leads to less overhead memory and more time required.
#'
#' @examples
#'
#' @author Christopher Nobles, Ph.D.
#'

cluster_multihits <- function(multihits, max_gap = 5L, iterations = 5L){
  stopifnot(class(multihits) == "GRanges")
  isThere <- names(GenomicRanges::mcols(multihits))
  if(!all(c("ID", "readPairKey") %in% isThere)){
    stop("ID and readPairKey columns not found in multihits input.")}

  read_key <- data.frame(
    "readPairKey" = sample(
      x = unique(GenomicRanges::mcols(multihits)$readPairKey),
      size = length(unique(GenomicRanges::mcols(multihits)$readPairKey))))

  if(iterations > 1){
    read_key$proc_group <- ceiling(
        seq_along(read_key$readPairKey)/(nrow(read_key)/iterations))
  }else{
    read_key$proc_group <- rep(1, nrow(read_key))
  }

  multihit_list <- split(
    multihits,
    read_key$proc_group[
      match(multihits$readPairKey, read_key$readPairKey)])

  axil_gr <- GRanges()
  edgelist <- matrix(ncol = 2)

  lapply(multihit_list, function(mhits){
    fl_mhits <- c(
      axil_gr, GenomicRanges::flank(mhits, width = -1, start = TRUE))
    red_mhits <- GenomicRanges::reduce(
      fl_mhits, min.gapwidth = max_gap, with.revmap = TRUE)
    revmap <- IRanges::as.list(red_mhits$revmap)

    axil_nodes <- as.character(S4Vectors::Rle(
      values = fl_mhits$readPairKey[sapply(revmap, "[", 1)],
      length = sapply(revmap, length)
    ))
    nodes <- fl_mhits$readPairKey[unlist(revmap)]
    el <- unique(matrix( c(axil_nodes, nodes), ncol = 2 ))
    clus <- igraph::clusters(igraph::graph.edgelist(el, directed = FALSE))
    clus_key <- data.frame(
      row.names = unique(as.character(t(el))),
      "clusID" = clus$membership)

    axils <- fl_mhits[fl_mhits$readPairKey %in% axil_nodes]
    axil_gr <<- unique(GenomicRanges::flank(axils, width = -1, start = TRUE))
    edgelist <<- rbind(edgelist, el)
  })

  edgelist <- S4Vectors::na.exclude(edgelist)
  cluster_data <- igraph::clusters(igraph::graph.edgelist(
    edgelist, directed = FALSE))
  key <- data.frame(row.names = unique(as.character(t(edgelist))),
                    "multihitID" = cluster_data$membership)
  multihits$multihitID <- key[multihits$readPairKey, "multihitID"]

  clusteredMultihitPositions <- split(
    GenomicRanges::flank(multihits, width = -1, start = TRUE),
    multihits$multihitID)
  clusteredMultihitNames <- lapply(
    clusteredMultihitPositions,
    function(x) unique(x$readPairKey) )
  clusteredMultihitPositions <- GenomicRanges::GRangesList(lapply(
    clusteredMultihitPositions,
    function(x) unname(unique(GenomicRanges::granges(x))) ))
  multihits_medians <- round(median(width(split(multihits, multihits$ID))))
  clusteredMultihitLengths <- lapply(
    clusteredMultihitNames,
    function(x){
      readIDs <- unique(multihits[multihits$readPairKey %in% x]$ID)
      df <- data.frame(table(multihits_medians[readIDs]))
      names(df) <- c("length", "count")
      df
  })

  list(
    "unclusteredMultihits" = multihits,
    "clusteredMultihitPositions" = clusteredMultihitPositions,
    "clusteredMultihitLengths" = clusteredMultihitLengths)
}

