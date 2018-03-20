#' Generate clusterIDs for integration sites given multiple gap distances.
#'
#' \code{serial_cluster} returns a GRanges object of the same length and order
#' as input with additional metadata columns specifying the group or clusterID
#' of given nt windows.
#'
#' @description Given a GRanges object, this function identifies groups or
#' clusters of integration sites within specified distances. The output is of
#' the same length and order of the input, and clusterID information is appended
#' to the metadata columns.
#'
#' @usage
#' serial_cluster(sites, gaps = c(0L, 1L, 2L))
#'
#' @param sites a GRanges object where each row is a different integration site
#' or range.
#'
#' @param gaps an integer vector of any length specifying the nucleotide
#' distance allowable between sites to call in the same group or cluster.
#'
#' @examples
#' sites <- GRanges(
#'   seqnames = rep("chr1", 7),
#'   ranges = IRanges(
#'     start = c(50, 49, 44, 50, 50, 50, 60),
#'     width = seq(20, 26, 1)),
#'   strand = rep("+", 7))
#'
#' serial_cluster(sites, gaps = c(1L, 5L, 10L))
#'
#' @author Christopher Nobles, Ph.D.
#'

serial_cluster <- function(sites, gaps = c(0L, 1L, 2L)){
  if(is.null(names(sites))){
    wasNull <- TRUE
    names(sites) <- seq(1:length(sites))
  }else{
    wasNull <- FALSE
  }
  origin_order <- names(sites)

  lapply(gaps, function(gap){
    clusID <- paste0("clusID.", gap, "nt")
    fl.sites <- GenomicRanges::flank(sites, -1, start = TRUE)
    red.sites <- GenomicRanges::reduce(
      fl.sites,
      min.gapwidth = gap,
      with.revmap = TRUE)

    revmap <- IRanges::as.list(red.sites$revmap)
    groups <- S4Vectors::Rle(
      values = 1:length(revmap),
      lengths = sapply(revmap, length)
    )
    mod.sites <- sites[unlist(revmap)]
    GenomicRanges::mcols(mod.sites)[, clusID] <- groups
    sites <<- mod.sites[origin_order]
  })
  if(wasNull) names(sites) <- NULL
  sites
}
