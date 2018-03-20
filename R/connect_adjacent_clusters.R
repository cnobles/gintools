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
#' connect_adjacent_clusters(red.sites, graph, gap, bais)
#'
#' @param red.sites GRanges object which has been reduced to single nt positions
#' and contains the revmap from the original GRanges object. The object must
#' also contain a column for cluster membership (clusID) and a column for
#' abundance (fragLengths).
#'
#' @param graph a directed graph built from the red.sites object. Each node
#' corresponds to a row in the red.sites object.
#'
#' @param bias either "upsteam" or "downstream", designating which position to
#' choose if other decision metrics are tied.
#'
#' @examples
#' gr <- gintools:::generate_test_granges(stdev = 3)
#' red.sites <- reduce(
#'   flank(gr, -1, start = TRUE),
#'   min.gapwidth = 0L,
#'   with.revmap = TRUE)
#' red.sites$siteID <- seq(1:length(red.sites))
#' revmap <- as.list(red.sites$revmap)
#' red.sites$fragLengths <- sapply(revmap, length)
#' red.hits <- GenomicRanges::as.data.frame(
#'   findOverlaps(red.sites, maxgap = 1L, drop.self = TRUE))
#' red.hits <- red.hits %>%
#'   mutate(q_pos = start(red.sites[queryHits])) %>%
#'   mutate(s_pos = start(red.sites[subjectHits])) %>%
#'   mutate(q_abund = red.sites[queryHits]$abundance) %>%
#'   mutate(s_abund = red.sites[subjectHits]$abundance) %>%
#'   mutate(strand = unique(strand(
#'     c(red.sites[queryHits], red.sites[subjectHits])))) %>%
#'   mutate(is.upstream = ifelse(
#'     strand == "+",
#'     q_pos < s_pos,
#'     q_pos > s_pos)) %>%
#'   mutate(keep = q_abund > s_abund) %>%
#'   mutate(keep = ifelse(
#'     q_abund == s_abund,
#'     is.upstream,
#'     keep)) %>%
#'   filter(keep)
#' g <- make_empty_graph(n = length(red.sites), directed = TRUE) %>%
#'   add_edges(unlist(mapply(
#'     c, red.hits$queryHits, red.hits$subjectHits, SIMPLIFY = FALSE)))
#' red.sites$clusID <- clusters(g)$membership
#' g <- connect_satalite_vertices(red.sites, g, gap = 2L, "upstream")
#' red.sites$clusID <- clusters(g)$membership
#' g <- break_connecting_source_paths(red.sites, g, "upstream")
#' red.sites$clusID <- clusters(g)$membership
#'
#' connect_adjacent_clusters(red.sites, g, gap = 5L, "upstream")
#'
#' @author Christopher Nobles, Ph.D.
#'
#' @importFrom magrittr %>%
#'

connect_adjacent_clusters <- function(red.sites, graph, gap, bias){
  src.nodes <- sources(graph)
  near.sources <- GenomicRanges::findOverlaps(
    red.sites[src.nodes],
    maxgap = gap,
    drop.self = TRUE
  )

  if(length(near.sources) > 0){
    # Identify sources of clusters within the largets satalite gap distance
    # and identify the directionality of the edge to create based first on
    # abundance (source will likely have greater abundance) and then by
    # upstream bias (more likely the origin site is upstream of drifting mapped
    # reads).
    if(bias == "upstream"){
      near.src.df <- data.frame(
        node.i = src.nodes[queryHits(near.sources)],
        node.j = src.nodes[subjectHits(near.sources)]
      ) %>%
        dplyr::mutate(abund.i = red.sites[node.i]$abundance) %>%
        dplyr::mutate(abund.j = red.sites[node.j]$abundance) %>%
        dplyr::mutate(pos.i = GenomicRanges::start(red.sites[node.i])) %>%
        dplyr::mutate(pos.j = GenomicRanges::start(red.sites[node.j])) %>%
        dplyr::mutate(strand = as.character(
          GenomicRanges::strand(red.sites[node.i]))) %>%
        dplyr::mutate(is.upstream = ifelse(
          strand == "+",
          pos.i < pos.j,
          pos.i > pos.j))
    }else if(bias == "downstream"){
      near.src.df <- data.frame(
        node.i = src.nodes[queryHits(near.sources)],
        node.j = src.nodes[subjectHits(near.sources)]
      ) %>%
        dplyr::mutate(abund.i = red.sites[node.i]$abundance) %>%
        dplyr::mutate(abund.j = red.sites[node.j]$abundance) %>%
        dplyr::mutate(pos.i = GenomicRanges::start(red.sites[node.i])) %>%
        dplyr::mutate(pos.j = GenomicRanges::start(red.sites[node.j])) %>%
        dplyr::mutate(strand = as.character(
          GenomicRanges::strand(red.sites[node.i]))) %>%
        dplyr::mutate(is.downstream = ifelse(
          strand == "+",
          pos.i > pos.j,
          pos.i < pos.j))
    }else{
      stop("No bias specified. Please choose either 'upstream' or 'downstream'.")
    }

    redundant_graph <- igraph::make_graph(
      edges = unlist(with(
        near.src.df,
        mapply(
          c,
          1:nrow(near.src.df),
          match(paste(node.i, node.j), paste(node.j, node.i))
    ))))
    redundant_groups <- igraph::clusters(redundant_graph)$membership

    if(bias == "upstream"){
      near.src.df <- dplyr::mutate(
          near.src.df, redundant_grp = redundant_groups) %>%
        dplyr::filter(abund.i >= abund.j) %>%
        dplyr::group_by(redundant_grp) %>%
        dplyr::mutate(group_size = n()) %>%
        dplyr::mutate(keep = ifelse(
          group_size == 1,
          TRUE,
          is.upstream)) %>%
        dplyr::filter(keep)
    }else if(bias == "downstream"){
      near.src.df <- dplyr::mutate(
          near.src.df, redundant_grp = redundant_groups) %>%
        dplyr::filter(abund.i >= abund.j) %>%
        dplyr::group_by(redundant_grp) %>%
        dplyr::mutate(group_size = n()) %>%
        dplyr::mutate(keep = ifelse(
          group_size == 1,
          TRUE,
          is.downstream)) %>%
        dplyr::filter(keep)
    }

    edges.to.connect.near.srcs <- unlist(with(
      near.src.df,
      mapply(c, node.i, node.j, SIMPLIFY = FALSE)
    ))
  }else{
    edges.to.connect.near.srcs <- c()
  }

  igraph::add.edges(graph, edges.to.connect.near.srcs)
}
