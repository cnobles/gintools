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
#' also contain a column for cluster membership (clus.id) and a column for
#' abundance (fragLengths).
#' @param graph a directed graph built from the red.sites object. Each node
#' corresponds to a row in the red.sites object.
#' @param gap integer nucleotide (nt) gap distance to consider joining to
#' clusters.
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
#' red.sites$fragLengths <- lengths(revmap)
#' red.hits <- GenomicRanges::as.data.frame(
#'   findOverlaps(red.sites, maxgap = 0L, drop.self = TRUE))
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
#' red.sites$clus.id <- clusters(g)$membership
#'
#' connect_satalite_vertices(red.sites, g, gap = 2L, "upstream")
#'
#' @author Christopher Nobles, Ph.D.
#'
#' @importFrom magrittr %>%

connect_satalite_vertices <- function(red.sites, graph, gap, bias){
  clus_mem <- igraph::clusters(graph)$membership
  
  clus_ranges <- unlist(GenomicRanges::reduce(
    GenomicRanges::split(red.sites, clus_mem),
    min.gapwidth = (gap-1)
  ))
  
  sata_hits <- as.data.frame(
    GenomicRanges::findOverlaps(
      clus_ranges, maxgap = gap - 1L, drop.self = TRUE
    )
  )
  
  names(sata_hits) <- c("source.clus", "sata.clus")

  red_df <- GenomicRanges::as.data.frame(red.sites)

  if(nrow(sata_hits) > 0){
    clus_data <- red_df %>%
      dplyr::group_by(clus.id) %>%
      dplyr::summarize(
        clus.pos.mean = as.integer(mean(start)),
        min.abund = min(abund),
        sum.abund = sum(abund))

    if(bias == "upstream"){
      sata_hits <- sata_hits %>%
        dplyr::mutate(
          source.pos = clus_data[source.clus,]$clus.pos.mean,
          sata.pos = clus_data[sata.clus,]$clus.pos.mean,
          min.src.abund = clus_data[.$source.clus,]$min.abund,
          min.sat.abund = clus_data[.$sata.clus,]$min.abund,
          sum.src.abund = clus_data[.$source.clus,]$sum.abund,
          sum.sat.abund = clus_data[.$sata.clus,]$sum.abund,
          is.upstream = source.pos < sata.pos) %>%
        dplyr::filter(
          as.integer(min.src.abund) >= as.integer(min.sat.abund),
          sum.src.abund > sum.sat.abund,
          abs(source.clus - sata.clus) == 1)
    }else if(bias == "downstream"){
      sata_hits <- sata_hits %>%
        dplyr::mutate(
          source.pos = clus_data[source.clus,]$clus.pos.mean,
          sata.pos = clus_data[sata.clus,]$clus.pos.mean,
          min.src.abund = clus_data[.$source.clus,]$min.abund,
          min.sat.abund = clus_data[.$sata.clus,]$min.abund,
          sum.src.abund = clus_data[.$source.clus,]$sum.abund,
          sum.sat.abund = clus_data[.$sata.clus,]$sum.abund,
          is.downstream = source.pos > sata.pos) %>%
        dplyr::filter(
          as.integer(min.src.abund) >= as.integer(min.sat.abund),
          sum.src.abund > sum.sat.abund,
          abs(source.clus - sata.clus) == 1)
    }else{
      stop("No bias specified. Please choose either 'upstream' or 'downstream'.")
    }

    if(nrow(sata_hits) > 0){
      clus.map <- GenomicRanges::findOverlaps(clus_ranges, red.sites)
      clus.list <- split(
        S4Vectors::subjectHits(clus.map), S4Vectors::queryHits(clus.map))

      if(bias == "upstream"){
        sata_hits <- sata_hits %>%
          dplyr::mutate(
            source.node = ifelse(
              sata_hits$is.upstream,
              sapply(clus.list[sata_hits$source.clus], dplyr::last),
              sapply(clus.list[sata_hits$source.clus], dplyr::first)),
            sata.node = ifelse(
              is.upstream,
              sapply(clus.list[sata.clus], dplyr::first),
              sapply(clus.list[sata.clus], dplyr::last)))
      }else if(bias == "downstream"){
        sata_hits <- sata_hits %>%
          dplyr::mutate(
            source.node = ifelse(
              sata_hits$is.downstream,
              sapply(clus.list[sata_hits$source.clus], dplyr::first),
              sapply(clus.list[sata_hits$source.clus], dplyr::last)),
            sata.node = ifelse(
              is.downstream,
              sapply(clus.list[sata.clus], dplyr::last),
              sapply(clus.list[sata.clus], dplyr::first)))
      }

      sata.edges <- with(sata_hits, vzip(source.node, sata.node))
    }else{
      sata.edges <- c()
    }
  }else{
    sata.edges <- c()
  }
  igraph::add_edges(graph, sata.edges)
}
