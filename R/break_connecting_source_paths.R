#' Break graph paths which connect sources.
#'
#' \code{break_connecting_source_paths} returns a graph where only one source
#' is present per cluster.
#'
#' @description Given a list of unique integration site positions (reduced
#' GRanges object) and a directed graph of connected components, this function
#' identifies clusters with multiple sources, the paths between those sources,
#' and removes edges along the path so that each cluster only has one source
#' node. Edge removal is first based on nucleotide distance (greater distance
#' prefered), then based on abundance (lowest abundance prefered), then on an
#' upstream bias (downstream connection will be removed when everything ties).
#'
#' @usage
#' break_connecting_source_paths(red.sites, graph, bias)
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
#' gr <- .generate_test_granges(stdev = 3)
#' red.sites <- reduce(
#'   flank(gr, -1, start = TRUE),
#'   min.gapwidth = 0L,
#'   with.revmap = TRUE)
#' red.sites$siteID <- seq(1:length(red.sites))
#' revmap <- as.list(red.sites$revmap)
#' red.sites$abundance <- sapply(revmap, length)
#' red.hits <- GenomicRanges::as.data.frame(
#'   findOverlaps(red.sites, maxgap = 1L, ignoreSelf = TRUE))
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
#' g <- connect_satalite_vertices(red.sites, g, gap = 2L, bias = "upstream")
#' red.sites$clusID <- clusters(g)$membership
#' break_connecting_source_paths(red.sites, g, "upstream")
#'
#' @author Christopher Nobles, Ph.D.
#' @export

break_connecting_source_paths <- function(red.sites, graph, bias){
  src.nodes <- sources(graph)
  sources.p.clus <- split(src.nodes, clusters(graph)$membership[src.nodes])
  clus.w.multi.sources <- sources.p.clus[sapply(sources.p.clus, length) > 1]

  if(length(clus.w.multi.sources) > 0){
    adj.pairs <- do.call(c, lapply(clus.w.multi.sources, function(x){
      lapply(1:(length(x)-1), function(i) c(x[i], x[i+1]))
    }))

    snk.nodes <- sinks(graph)

    edges.to.edit <- data.frame(
      "src_node_i" = sapply(adj.pairs, "[[", 1),
      "src_node_j" = sapply(adj.pairs, "[[", 2)
    ) %>%
      dplyr::mutate("src_node_i_abund" = as.numeric(red.sites[src_node_i]$abundance)) %>%
      dplyr::mutate("src_node_j_abund" = as.numeric(red.sites[src_node_j]$abundance))

    source.paths <- mapply(function(x,y){
      all_simple_paths(as.undirected(graph), x, y)},
      edges.to.edit$src_node_i,
      edges.to.edit$src_node_j)

    edges.to.edit <- mutate(edges.to.edit, sink_node = snk.nodes[
      sapply(source.paths, function(x){
        which(snk.nodes %in% x)
      })])

    # Identify the nodes adjacent to sinks between connected sources
    # then filter adjacent pairs to identify which edge should be 'clipped'.
    # Filtering based first on adjacent node distance (edges with greater
    # distance get clipped), then abundance (lower abund gets clipped), then
    # biasing on upstream edges over downstream (downstream is clipped for
    # tie breaking).

    if(bias == "upstream"){
      target.edges <- bind_rows(lapply(1:nrow(edges.to.edit), function(i){
        sink <- edges.to.edit[i, "sink_node"]
        path <- as.numeric(source.paths[[i]])
        pos <- grep(sink, path)
        data.frame(
          "sink" = rep(sink, 2),
          "adj_node" = c(path[pos-1], path[pos+1])
        )
      })) %>%
        dplyr::mutate(sink_pos = start(red.sites[sink])) %>%
        dplyr::mutate(adj_pos = start(red.sites[adj_node])) %>%
        dplyr::mutate(adj_abund = red.sites[adj_node]$abundance) %>%
        dplyr::mutate(nt_dist = abs(sink_pos - adj_pos)) %>%
        dplyr::mutate(strand = as.character(strand(red.sites[sink]))) %>%
        dplyr::mutate(is.upstream = ifelse(
          strand == "+",
          sink_pos < adj_pos,
          sink_pos > adj_pos)) %>%
        group_by(sink) %>%
        filter(nt_dist == max(nt_dist)) %>%
        filter(adj_abund == min(adj_abund)) %>%
        dplyr::mutate(group_size = n()) %>%
        dplyr::mutate(keep = ifelse(
            group_size == 1,
            TRUE,
            !is.upstream)) %>%
        filter(keep)
    }else if(bias == "downstream"){
      target.edges <- bind_rows(lapply(1:nrow(edges.to.edit), function(i){
        sink <- edges.to.edit[i, "sink_node"]
        path <- as.numeric(source.paths[[i]])
        pos <- grep(sink, path)
        data.frame(
          "sink" = rep(sink, 2),
          "adj_node" = c(path[pos-1], path[pos+1])
        )
      })) %>%
        dplyr::mutate(sink_pos = start(red.sites[sink])) %>%
        dplyr::mutate(adj_pos = start(red.sites[adj_node])) %>%
        dplyr::mutate(adj_abund = red.sites[adj_node]$abundance) %>%
        dplyr::mutate(nt_dist = abs(sink_pos - adj_pos)) %>%
        dplyr::mutate(strand = as.character(strand(red.sites[sink]))) %>%
        dplyr::mutate(is.downstream = ifelse(
          strand == "+",
          sink_pos > adj_pos,
          sink_pos < adj_pos)) %>%
        group_by(sink) %>%
        filter(nt_dist == max(nt_dist)) %>%
        filter(adj_abund == min(adj_abund)) %>%
        dplyr::mutate(group_size = n()) %>%
        dplyr::mutate(keep = ifelse(
          group_size == 1,
          TRUE,
          !is.downstream)) %>%
        filter(keep)
    }else{
      stop("No bias specified. Please choose either 'upstream' or 'downstream'.")
    }

    break.edges <- unlist(with(
      target.edges,
      mapply(c, sink, adj_node, SIMPLIFY = FALSE)
    ))

    edge.ids.to.break <- get.edge.ids(graph, break.edges, directed = FALSE)
  }else{
    edge.ids.to.break <- c()
  }

  delete_edges(graph, edge.ids.to.break)
}
