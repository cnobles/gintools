#' Resolve primary sources from clusters with multiple souce nodes.
#'
#' \code{resolve_cluster_sources} returns a graph where each cluster only
#' has a single primary source node.
#'
#' @description Given a list of unique integration site positions (reduced
#' GRanges object) and a directed graph of connected components, this function
#' identifies clusters where multiple source nodes exist and then identifies
#' which source should be considered the primary source node, first based on
#' abundance and then
#'
#' @usage
#' resolve_cluster_sources(red.sites, graph)
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
#' red.sites$siteID <- seq_along(red.sites)
#' revmap <- as.list(red.sites$revmap)
#' red.sites$abundance <- lengths(revmap)
#' red.hits <- GenomicRanges::as.data.frame(
#'   findOverlaps(red.sites, maxgap = 0L, drop.self = TRUE))
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
#' g <- connect_adjacent_clusters(red.sites, g, gap = 5L, "upstream")
#' red.sites$clusID <- clusters(g)$membership
#'
#' resolve_cluster_sources(red.sites, g, "upstream")
#'
#' @author Christopher Nobles, Ph.D.
#'
#' @importFrom magrittr %>%

resolve_cluster_sources <- function(red.sites, graph, bias){
  src_nodes <- sources(graph)
  sources_p_clus <- IRanges::IntegerList(split(
    src_nodes, igraph::clusters(graph)$membership[src_nodes]))
  clus_w_multi_sources <- sources_p_clus[S4Vectors::lengths(sources_p_clus) > 1]

  if(length(clus_w_multi_sources) > 0){
    if(bias == "upstream"){
      resolve_df <- data.frame(
          node = unlist(clus_w_multi_sources),
          clus = as.numeric(S4Vectors::Rle(
            values = seq_along(clus_w_multi_sources),
            lengths = S4Vectors::lengths(clus_w_multi_sources)))) %>%
        dplyr::mutate(abund = red.sites[node]$abund) %>%
        dplyr::group_by(clus) %>%
        dplyr::mutate(
          top.abund = abund == max(abund),
          strand = as.character(GenomicRanges::strand(red.sites[node])),
          pos = GenomicRanges::start(red.sites[node])) %>%
        dplyr::group_by(clus, top.abund) %>%
        dplyr::mutate(
          grp.size = n(),
          is.upstream = ifelse(strand == "+", pos == min(pos), pos == max(pos)),
          src = ifelse(
            top.abund == TRUE & grp.size > 1, is.upstream, top.abund)) %>%
        dplyr::ungroup() %>%
        as.data.frame()
    }else if(bias == "downstream"){
      resolve_df <- data.frame(
        node = unlist(clus_w_multi_sources),
        clus = as.numeric(S4Vectors::Rle(
          values = seq_along(clus_w_multi_sources),
          lengths = S4Vectors::lengths(clus_w_multi_sources)))) %>%
        dplyr::mutate(
          abund = red.sites[node]$abund) %>%
        dplyr::group_by(clus) %>%
        dplyr::mutate(
          top.abund = abund == max(abund),
          strand = as.character(GenomicRanges::strand(red.sites[node])),
          pos = GenomicRanges::start(red.sites[node])) %>%
        dplyr::group_by(clus, top.abund) %>%
        dplyr::mutate(
          grp.size = n(),
          is.downstream = ifelse(
            strand == "+", pos == max(pos), pos == min(pos)),
          src = ifelse(
            top.abund == TRUE & grp.size > 1, is.downstream, top.abund)) %>%
        dplyr::ungroup() %>%
        as.data.frame()
    }else{
      stop("No bias specified. Please choose either 'upstream' or 'downstream'.")
    }

    src_nodes <- resolve_df[resolve_df$src == TRUE,]$node
    snk_nodes <- lapply(seq_along(clus_w_multi_sources), function(i){
      resolve_df[resolve_df$clus == i & resolve_df$src == FALSE,]$node
    })

    # Accomidates multiple sinks for singular sources in a cluster.
    resolve_edges <- do.call(c, lapply(seq_along(src_nodes), function(i){
      src <- src_nodes[i]
      snk <- snk_nodes[[i]]
      vzip(rep(src, length(snk)), snk) }))
  }else{
    resolve_edges <- c()
  }

  igraph::add.edges(graph, resolve_edges)
}
