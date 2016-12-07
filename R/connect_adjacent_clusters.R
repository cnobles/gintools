#' Connect adjacent clusters when their sources are within a specified distance.
#'
#' \code{connect_adjacent_clusters} returns a graph where adjacent clusters with
#' sources within the gap distance are joined.
#'
#' @description Given a list of unique integration site positions (reduced
#' GRanges object) and a directed graph of connected components, this function
#' identifies clusters where source nodes are within a gap distance away from
#' each other, and connects the source nodes with a directional edge, based
#' first on abundance, and secondly on an upstream bias for tie breaking.
#'
#' @usage
#' connect_adjacent_clusters(red.sites, graph, gap)
#'
#' @param red.sites GRanges object which has been reduced to single nt positions
#' and contains the revmap from the original GRanges object. The object must
#' also contain a column for cluster membership (clusID) and a column for
#' abundance (fragLengths).
#'
#' @param graph a directed graph built from the red.sites object. Each node
#' corresponds to a row in the red.sites object.
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
#' red.hits <- red.hits %>%
#'   mutate(q_pos = start(red.sites[queryHits])) %>%
#'   mutate(s_pos = start(red.sites[subjectHits])) %>%
#'   mutate(q_fragLengths = red.sites[queryHits]$fragLengths) %>%
#'   mutate(s_fragLengths = red.sites[subjectHits]$fragLengths) %>%
#'   mutate(strand = unique(strand(
#'     c(red.sites[queryHits], red.sites[subjectHits])))) %>%
#'   mutate(is.upstream = ifelse(
#'     strand == "+",
#'     q_pos < s_pos,
#'     q_pos > s_pos)) %>%
#'   mutate(keep = q_fragLengths > s_fragLengths) %>%
#'   mutate(keep = ifelse(
#'     q_fragLengths == s_fragLengths,
#'     is.upstream,
#'     keep)) %>%
#'   filter(keep)
#' g <- make_empty_graph(n = length(red.sites), directed = TRUE) %>%
#'   add_edges(unlist(mapply(
#'     c, red.hits$queryHits, red.hits$subjectHits, SIMPLIFY = FALSE)))
#' red.sites$clusID <- clusters(g)$membership
#' g <- connect_satalite_vertices(red.sites, g, gap = 2L)
#' red.sites$clusID <- clusters(g)$membership
#' g <- break_connecting_source_paths(red.sites, g)
#' red.sites$clusID <- clusters(g)$membership
#'
#' connect_adjacent_clusters(red.sites, g, gap = 5L)
#'
#' @author Christopher Nobles, Ph.D.
#' @export

connect_adjacent_clusters <- function(red.sites, graph, gap){
  src.nodes <- sources(graph)
  near.sources <- findOverlaps(
    red.sites[src.nodes],
    maxgap = gap,
    ignoreSelf = TRUE
  )

  if(length(near.sources) > 0){
    # Identify sources of clusters within the largets satalite gap distance
    # and identify the directionality of the edge to create based first on
    # abundance (source will likely have greater abundance) and then by
    # upstream bias (more likely the origin site is upstream of drifting mapped
    # reads).
    near.src.df <- data.frame(
      node.i = src.nodes[queryHits(near.sources)],
      node.j = src.nodes[subjectHits(near.sources)]
    ) %>%
      dplyr::mutate(abund.i = red.sites[node.i]$fragLengths) %>%
      dplyr::mutate(abund.j = red.sites[node.j]$fragLengths) %>%
      dplyr::mutate(pos.i = start(red.sites[node.i])) %>%
      dplyr::mutate(pos.j = start(red.sites[node.j])) %>%
      dplyr::mutate(strand = as.character(strand(red.sites[node.i]))) %>%
      dplyr::mutate(is.upstream = ifelse(
        strand == "+",
        pos.i < pos.j,
        pos.i > pos.j))

    redundant_graph <- make_graph(
      edges = unlist(with(
        near.src.df,
        mapply(
          c,
          1:nrow(near.src.df),
          match(paste(node.i, node.j), paste(node.j, node.i))
    ))))
    redundant_groups <- clusters(redundant_graph)$membership

    near.src.df <- mutate(near.src.df, redundant_grp = redundant_groups) %>%
      filter(abund.i >= abund.j) %>%
      group_by(redundant_grp) %>%
      dplyr::mutate(group_size = n()) %>%
      dplyr::mutate(keep = ifelse(
        group_size == 1,
        TRUE,
        is.upstream)) %>%
      filter(keep)

    edges.to.connect.near.srcs <- unlist(with(
      near.src.df,
      mapply(c, node.i, node.j, SIMPLIFY = FALSE)
    ))
  }else{
    edges.to.connect.near.srcs <- c()
  }

  add.edges(graph, edges.to.connect.near.srcs)
}
