#' Scan subject and query ranges for windows of enrichment using the gRxCluster
#' function
#'
#' @usage gintools:::.scan_clusters(subject, query, grouping, kvals, nperm, ...)
#'
#' @param subject GRanges object (group = 0 / FALSE)
#' @param query GRanges object (group = 1 / TRUE)
#' @param grouping character Name of metadata column to use for grouping, if
#' NULL, sites will all be treated as the same group.
#' @param kvals integer vector of window widths
#' @param nperm number of permutations for FDR calculation
#' @param ... other args to pass to gRxCluster
#' @details Wrapper around gRxCluster from the geneRxCluster package.
#' @examples
#' gr1 <- .generate_test_granges(n_sites = 50)
#' gr2 <- .generate_test_granges(n_sites = 50)
#' .scan_clusters(gr1, gr2, kvals = 5L:10L, nperm = 20)
#' @author Christopher Nobles, Ph.D.
#'

.scan_clusters <- function(subject, query, grouping = NULL, kvals, nperm, ...){
  require(GenomicRanges)
  require(dplyr)
  require(geneRxCluster)

  subject <- .quick_clus_format(subject, grouping)
  query <- .quick_clus_format(query, grouping)
  subject$grp <- 0
  query$grp <- 1
  data <- rbind(subject, query)
  data <- arrange(data, chr, pos)
  gRxCluster(
    object = data$chr,
    starts = data$pos,
    group = data$grp,
    kvals = kvals,
    nperm = nperm,
    ... = ...
  )
}
