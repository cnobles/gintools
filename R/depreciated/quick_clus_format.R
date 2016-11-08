#' Quickly format a GRanges object for geneRxCluster
#'
#' @usage .quick_clus_format(gr)
#'
#' @param gr GRanges object consistent with integration sites
#' @param grouping character name of metadata column in GRanges object to be
#' used for defining groups. Groups will be made unique independently, meaning
#' if the same event is observed in two different groups, they are both kept,
#' but if they occur in the same group, they will be made one.
#' @details Quickly formats a GRanges object for use with the geneRxCluster
#' functions, specifically the scan statistics of gRxCluster.
#' @examples
#' gr <- .generate_test_granges()
#' .quick_clus_format(gr, grouping = NULL)
#' @author Christopher Nobles, Ph.D.
#'

.quick_clus_format <- function(gr, grouping = NULL){
  if(is.null(grouping)) gr$grp <- grouping <- "grp"
  group_col <- grep(grouping, names(mcols(gr)))
  gr <- sort(gr)
  df <- data.frame(
    "chr" = seqnames(gr),
    "pos" = ifelse(strand(gr) == "+", start(gr), end(gr)),
    "grp" = mcols(gr)[, group_col]
  )
  df <- distinct(df)
  df[, c(1,2)]
}
