#' Determine sources for clusters where they are absent.
#'
#' \code{resolve_source_deficiency} returns a graph where each cluster has at
#' least one source.
#'
#' @description Given a list of unique integration site positions (reduced
#' GRanges object) and a directed graph of connected components, this function
#' identifies clusters missing source nodes, the edges to remove based on the
#' bias, and then returns the edited graph.
#'
#' @usage
#' resolve_source_deficiency(red.sites, graph, bias = "upstream")
#'
#' @param red.sites GRanges object which has been reduced to single nt positions
#' and contains the revmap from the original GRanges object. The object must
#' also contain a column for cluster membership (clusID) and a column for
#' abundance (fragLengths).
#'
#' @param graph a directed graph built from the red.sites object. Each node
#' corresponds to a row in the red.sites object.
#'
#' @param bias character either "upstream" or "downstream" dictates a bias
#' imposed on the analysis. Selecting one will preference the direction when
#' breaking redundant edges.
#'
#' @examples
#' gr <- .generate_test_granges(stdev = 3)
#' red.sites <- reduce(
#'   flank(gr, -1, start = TRUE),
#'   min.gapwidth = 0L,
#'   with.revmap = TRUE)
#' red.sites$siteID <- seq(1:length(red.sites))
#' revmap <- as.list(red.sites$revmap)
#' red.sites$fragLengths <- sapply(revmap, length)
#' red.hits <- GenomicRanges::as.data.frame(
#'   findOverlaps(red.sites, maxgap = 1L, ignoreSelf = TRUE))
#' g <- make_empty_graph(n = length(red.sites), directed = TRUE) %>%
#'   add_edges(unlist(mapply(
#'     c, red.hits$queryHits, red.hits$subjectHits, SIMPLIFY = FALSE)))
#' red.sites$clusID <- clusters(g)$membership
#'
#' resolve_source_deficiency(red.sites, g)
#'
#' @author Christopher Nobles, Ph.D.


resolve_source_deficiency <- function(red.sites, graph, bias = "upstream"){
  clus <- clusters(graph)$membership
  clus.w.srcs <- clus[sources(graph)]
  clus.wo.srcs <- unique(clus[!clus %in% clus.w.srcs])

  if(length(clus.wo.srcs) > 0){
    gr <- red.sites[red.sites$clusID %in% clus.wo.srcs]
    red.sites$pos <- ifelse(
      strand(red.sites) == "+",
      start(red.sites),
      end(red.sites))

    edges <- igraph::as_data_frame(graph, "edges")
    df <- edges %>%
      filter(from %in% gr$siteID) %>%
      filter(to %in% gr$siteID) %>%
      dplyr::mutate(from.pos = red.sites$pos[from]) %>%
      dplyr::mutate(to.pos = red.sites$pos[to]) %>%
      dplyr::mutate(strand = unique(strand(c(red.sites[from], red.sites[to])))) %>%
      dplyr::mutate(is.upstream = ifelse(
        strand == "+",
        from.pos < to.pos,
        from.pos > to.pos))

    if(bias == "upstream"){
      df <- df[df$is.upstream == FALSE,]
    }else{
      df <- df[df$is.upstream == TRUE,]
    }

    edges.to.remove <- unlist(mapply(c, df$from, df$to, SIMPLIFY = FALSE))

    edge.ids.to.remove <- get.edge.ids(graph, edges.to.remove)
  }else{
    edge.ids.to.remove <- c()
  }

  delete_edges(graph, edge.ids.to.remove)
}
