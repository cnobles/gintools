#' Format a set of GRanges objects for geneRxCluster
#'
#' @usage scan_format(query, subject, grouping = NULL)
#'
#' @param query GRanges object consistent with integration sites, will be
#' associated with grp = 0.
#' @param subject GRanges object consistent with integration sites, will be
#' associated with grp = 1.
#' @param grouping character name of metadata column in GRanges object to be
#' used for defining groups. Groups will be made unique independently, meaning
#' if the same event is observed in two different groups, they are both kept,
#' but if they occur in the same group, they will be made one.
#'
#' @details Formats a GRanges object for use with the geneRxCluster
#' functions, specifically the scan statistics of gRxCluster.
#'
#' @examples
#' gr1 <- .generate_test_granges()
#' gr2 <- .generate_test_granges()
#' scan_format(gr1, gr2)
#'
#' @author Christopher Nobles, Ph.D.
#' @export

scan_format <- function(query, subject, grouping = NULL){
  format <- function(gr, grouping = grouping){
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

  subject <- format(subject, grouping)
  query <- format(query, grouping)
  subject$grp <- 0
  query$grp <- 1
  data <- rbind(subject, query) %>% arrange(chr, pos) %>% as.data.frame(.)
  data
}

