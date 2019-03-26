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
#' gr <- generate_test_granges(stdev = 3)
#' red.sites <- reduce(
#'   flank(gr, -1, start = TRUE),
#'   min.gapwidth = 0L,
#'   with.revmap = TRUE)
#' red.sites$siteID <- seq(1:length(red.sites))
#' revmap <- as.list(red.sites$revmap)
#' red.sites$abundance <- sapply(revmap, length)
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
#' g <- connect_satalite_vertices(red.sites, g, gap = 2L, bias = "upstream")
#' red.sites$clusID <- clusters(g)$membership
#' break_connecting_source_paths(red.sites, g, "upstream")
#'
#' @author Christopher Nobles, Ph.D.
#'
#' @importFrom magrittr %>%
#'

break_connecting_source_paths <- function(red.sites, graph, bias){
  src_nodes <- sources(graph)
  sources_p_clus <- IRanges::IntegerList(split(
    src_nodes, igraph::clusters(graph)$membership[src_nodes]))
  clus_w_multi_sources <- sources_p_clus[S4Vectors::lengths(sources_p_clus) > 1]

  if(length(clus_w_multi_sources) > 0){
    adj_pairs <- do.call(c, lapply(clus_w_multi_sources, function(x){
      lapply(1:(length(x)-1), function(i) c(x[i], x[i+1]))
    }))

    snk_nodes <- sinks(graph)

    edges_to_edit <- data.frame(
        src.node.i = unlist(adj_pairs)[
          IRanges::start(IRanges::IntegerList(adj_pairs)@partitioning)],
        src.node.j = unlist(adj_pairs)[
          IRanges::end(IRanges::IntegerList(adj_pairs)@partitioning)]) %>%
      dplyr::mutate(
        src.node.i.abund = as.numeric(red.sites[src.node.i]$abund),
        src.node.j.abund = as.numeric(red.sites[src.node.j]$abund),
        sink.node = IRanges::start(
          IRanges::findOverlapPairs(
            IRanges::IRanges(src.node.i, src.node.j),
            IRanges::IRanges(snk_nodes, width = 1))@second))

    # Identify the nodes adjacent to sinks between connected sources
    # then filter adjacent pairs to identify which edge should be 'clipped'.
    # Filtering based first on adjacent node distance (edges with greater
    # distance get clipped), then abundance (lower abund gets clipped), then
    # biasing on upstream edges over downstream (downstream is clipped for
    # tie breaking).

    if(bias == "upstream"){
      target_edges <- dplyr::bind_rows(lapply(
          seq_len(nrow(edges_to_edit)), function(i){
            sink <- edges_to_edit[i, "sink.node"]
            path <- unlist(igraph::all_simple_paths(
              igraph::as.undirected(graph),
              edges_to_edit[i, "src.node.i"],
              edges_to_edit[i, "src.node.j"]))
            pos <- which(path == sink)
            data.frame(
              sink = rep(sink, 2),
              adj.node = c(path[pos-1], path[pos+1])) })) %>%
        dplyr::mutate(
          sink.pos = GenomicRanges::start(red.sites[sink]),
          adj.pos = GenomicRanges::start(red.sites[adj.node]),
          adj.abund = red.sites[adj.node]$abund,
          nt.dist = abs(sink.pos - adj.pos),
          strand = as.character(GenomicRanges::strand(red.sites[sink])),
          is.upstream = ifelse(
            strand == "+", sink.pos < adj.pos, sink.pos > adj.pos)) %>%
        dplyr::group_by(sink) %>%
        dplyr::filter(nt.dist == max(nt.dist)) %>%
        dplyr::filter(adj.abund == min(adj.abund)) %>%
        dplyr::mutate(
          group_size = n(),
          keep = ifelse(group_size == 1, TRUE, !is.upstream)) %>%
        dplyr::filter(keep) %>%
        dplyr::ungroup() %>%
        as.data.frame()
    }else if(bias == "downstream"){
      target_edges <- dplyr::bind_rows(lapply(
          seq_len(nrow(edges_to_edit)), function(i){
            sink <- edges_to_edit[i, "sink.node"]
            path <- unlist(igraph::all_simple_paths(
              igraph::as.undirected(graph),
              edges_to_edit[i, "src.node.i"],
              edges_to_edit[i, "src.node.j"]))
            pos <- which(path == sink)
            data.frame(
              sink = rep(sink, 2),
              adj.node = c(path[pos-1], path[pos+1])) })) %>%
        dplyr::mutate(
          sink.pos = GenomicRanges::start(red.sites[sink]),
          adj.pos = GenomicRanges::start(red.sites[adj.node]),
          adj.abund = red.sites[adj.node]$abund,
          nt.dist = abs(sink.pos - adj.pos),
          strand = as.character(
            GenomicRanges::strand(red.sites[sink])),
          is.downstream = ifelse(
            strand == "+", sink.pos > adj.pos, sink.pos < adj.pos)) %>%
        dplyr::group_by(sink) %>%
        dplyr::filter(nt.dist == max(nt.dist)) %>%
        dplyr::filter(adj.abund == min(adj.abund)) %>%
        dplyr::mutate(
          group_size = n(),
          keep = ifelse(group_size == 1, TRUE, !is.downstream)) %>%
        dplyr::filter(keep) %>%
        dplyr::ungroup() %>%
        as.data.frame()
    }else{
      stop("No bias specified. Please choose either 'upstream' or 'downstream'.")
    }

    break_edges <- with(target_edges, vzip(sink, adj.node))

    edge_ids_to_break <- igraph::get.edge.ids(
      graph, break_edges, directed = FALSE)
  }else{
    edge_ids_to_break <- c()
  }

  igraph::delete_edges(graph, edge_ids_to_break)
}
