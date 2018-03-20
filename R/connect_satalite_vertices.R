#' Connect integration sites which lie a specific gap distance away from a
#' cluster of integration sites.
#'
#' \code{connect_satalite_vertices} returns a graph where nodes within 'gap'
#' distance from clusters are now connected to the edge of the clusters.
#'
#' @description Given a list of unique integration site positions (reduced
#' GRanges object) and a directed graph of connected components, this function
#' identifies each node within gap range of clusters and creates an edge to
#' connect the cluster to the 'satalite' node. Edges are drawn from the last
#' node in the cluster to the 'satalite' node, but directionality is determined
#' first by abundance and secondly by an upstream bias.
#'
#' @usage
#' connect_satalite_vertices(red.sites, graph, gap, bais)
#'
#' @param red.sites GRanges object which has been reduced to single nt positions
#' and contains the revmap from the original GRanges object. The object must
#' also contain a column for cluster membership (clusID) and a column for
#' abundance (fragLengths).
#'
#' @param graph a directed graph built from the red.sites object. Each node
#' corresponds to a row in the red.sites object.
#'
#' @param gap integer nucleotide (nt) gap distance to consider joining to
#' clusters.
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
#'
#' connect_satalite_vertices(red.sites, g, gap = 2L, "upstream")
#'
#' @author Christopher Nobles, Ph.D.
#'
#' @importFrom magrittr %>%

connect_satalite_vertices <- function(red.sites, graph, gap, bias){
  clus_mem <- igraph::clusters(graph)$membership
  clus.ranges <- unlist(GenomicRanges::reduce(
    GenomicRanges::split(red.sites, clus_mem),
    min.gapwidth = (gap-1)))
  sata.hits <- as.data.frame(
    GenomicRanges::findOverlaps(clus.ranges, maxgap = gap, drop.self = TRUE))
  names(sata.hits) <- c("source_clus", "sata_clus")

  red.df <- GenomicRanges::as.data.frame(red.sites)

  if(nrow(sata.hits) > 0){
    clus.data <- red.df %>%
      dplyr::group_by(clusID) %>%
      dplyr::summarize(
        clus_pos_mean = as.integer(mean(start)),
        min_abund = min(abundance),
        sum_abund = sum(abundance))

    if(bias == "upstream"){
      sata.hits <- sata.hits %>%
        dplyr::mutate(source_pos = clus.data[source_clus,]$clus_pos_mean) %>%
        dplyr::mutate(sata_pos = clus.data[sata_clus,]$clus_pos_mean) %>%
        dplyr::mutate(min_src_abund = clus.data[.$source_clus,]$min_abund) %>%
        dplyr::mutate(min_sat_abund = clus.data[.$sata_clus,]$min_abund) %>%
        dplyr::mutate(sum_src_abund = clus.data[.$source_clus,]$sum_abund) %>%
        dplyr::mutate(sum_sat_abund = clus.data[.$sata_clus,]$sum_abund) %>%
        dplyr::mutate(is_upstream = source_pos < sata_pos) %>%
        dplyr::filter(
          as.integer(min_src_abund) >= as.integer(min_sat_abund)) %>%
        dplyr::filter(sum_src_abund > sum_sat_abund) %>%
        dplry::filter(abs(source_clus - sata_clus) == 1)
    }else if(bias == "downstream"){
      sata.hits <- sata.hits %>%
        dplyr::mutate(source_pos = clus.data[source_clus,]$clus_pos_mean) %>%
        dplyr::mutate(sata_pos = clus.data[sata_clus,]$clus_pos_mean) %>%
        dplyr::mutate(min_src_abund = clus.data[.$source_clus,]$min_abund) %>%
        dplyr::mutate(min_sat_abund = clus.data[.$sata_clus,]$min_abund) %>%
        dplyr::mutate(sum_src_abund = clus.data[.$source_clus,]$sum_abund) %>%
        dplyr::mutate(sum_sat_abund = clus.data[.$sata_clus,]$sum_abund) %>%
        dplyr::mutate(is_downstream = source_pos > sata_pos) %>%
        dplyr::filter(
          as.integer(min_src_abund) >= as.integer(min_sat_abund)) %>%
        dplyr::filter(sum_src_abund > sum_sat_abund) %>%
        dplyr::filter(abs(source_clus - sata_clus) == 1)
    }else{
      stop("No bias specified. Please choose either 'upstream' or 'downstream'.")
    }

    if(nrow(sata.hits) > 0){
      clus.map <- GenomicRanges::findOverlaps(clus.ranges, red.sites)
      clus.list <- split(subjectHits(clus.map), queryHits(clus.map))

      if(bias == "upstream"){
        sata.hits <- sata.hits %>%
          dplyr::mutate(source_node = ifelse(
            sata.hits$is_upstream,
            sapply(clus.list[sata.hits$source_clus], dplyr::last),
            sapply(clus.list[sata.hits$source_clus], dplyr::first)
          )) %>%
          dplyr::mutate(sata_node = ifelse(
            is_upstream,
            sapply(clus.list[sata_clus], dplyr::first),
            sapply(clus.list[sata_clus], dplyr::last)
          ))
      }else if(bias == "downstream"){
        sata.hits <- sata.hits %>%
          dplyr::mutate(source_node = ifelse(
            sata.hits$is_downstream,
            sapply(clus.list[sata.hits$source_clus], dplyr::first),
            sapply(clus.list[sata.hits$source_clus], dplyr::last)
          )) %>%
          dplyr::mutate(sata_node = ifelse(
            is_downstream,
            sapply(clus.list[sata_clus], dplyr::last),
            sapply(clus.list[sata_clus], dplyr::first)
          ))
      }

      sata.edges <- unlist(with(
        sata.hits,
        mapply(c, source_node, sata_node, SIMPLIFY = FALSE)
      ))
    }else{
      sata.edges <- c()
    }
  }else{
    sata.edges <- c()
  }
  igraph::add_edges(graph, sata.edges)
}
