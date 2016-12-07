#' Standardize genomic integration site positions within a dataset
#'
#' \code{standardize_intsites} returns a GRanges object where the integration
#' site positions have been standardized with sites within the gap distance.
#'
#' @description Given a GRanges object ...
#'
#' @usage
#' standardize_intsites(sites)
#'
#' standardize_intsites(sites, min.gap = 1L, sata.gap = 5L)
#'
#' @param sites GRanges Object integration sites to standardize.
#' @param min.gap integer minimum gap to consider combine integration sites.
#' @param sata.gap integer maximum distance to consider combining integration
#' sites.
#'
#' @examples
#' gr <- .generate_test_granges()
#'
#' standardize_intsites(gr)
#'
#' @author Christopher Nobles, Ph.D.
#' @export

standardize_intsites <- function(sites, min.gap = 1L, sata.gap = 5L){
  require(gintools)

  # Start by reducing the sites object down to only the "intSite" positions
  # and storing the revmap for going back to the original sites.

  message(paste0("Standardizing ", length(sites), " positions."))
  message("Generating initial graph by connecting positions with 1 nt difference.")

  red.sites <- reduce(
    flank(sites, -1, start = TRUE),
    min.gapwidth = 0L,
    with.revmap = TRUE)
  red.sites$siteID <- seq(1:length(red.sites))
  revmap <- as.list(red.sites$revmap)
  red.sites$fragLengths <- sapply(revmap, length)
  # Not true if original sites are not dereplicated or unique
  red.hits <- GenomicRanges::as.data.frame(
    findOverlaps(red.sites, maxgap = min.gap, ignoreSelf = TRUE))

  red.hits <- red.hits %>%
    dplyr::mutate(q_pos = start(red.sites[queryHits])) %>%
    dplyr::mutate(s_pos = start(red.sites[subjectHits])) %>%
    dplyr::mutate(q_fragLengths = red.sites[queryHits]$fragLengths) %>%
    dplyr::mutate(s_fragLengths = red.sites[subjectHits]$fragLengths) %>%
    dplyr::mutate(strand = as.vector(strand(red.sites[queryHits]))) %>%
    dplyr::mutate(is.upstream = ifelse(
      strand == "+",
      q_pos < s_pos,
      q_pos > s_pos)) %>%
    dplyr::mutate(keep = q_fragLengths > s_fragLengths) %>%
    dplyr::mutate(keep = ifelse(
      q_fragLengths == s_fragLengths,
      is.upstream,
      keep)) %>%
    filter(keep)

  # Construct initial graph, the graph will be modified during other parts of
  # the function and can be used to identify clusters as well as 'sources'
  # within those clusters as it's a directed graph. Sources are identified and
  # used to identify the original position given a distribution.
  g <- make_empty_graph(n = length(red.sites), directed = TRUE) %>%
    add_edges(unlist(mapply(
      c, red.hits$queryHits, red.hits$subjectHits, SIMPLIFY = FALSE)))

  red.sites$clusID <- clusters(g)$membership

  message(paste0("Initial cluster count: ", clusters(g)$no))

  # When clusters are formed of positions with equivalent abundance, both
  # directional edges are created. This type of cluster does not have sources
  # or sinks. The following set resolves which directed edges to remove, but
  # imposes an upstream bias.

#  g <- resolve_source_deficiency(red.sites, g, bias = "upstream")
#  red.sites$clusID <- clusters(g)$membership

  # Identify satalite positions that should be included in clusters up to 5nt
  # away. This portion of the function tries to reach out to up to 5nt from the
  # boundry of a cluster to see if there are any further "satalite" positions
  # that have been annotated. It does this by iterively increasing the size
  # from 2nt to 5nt by 1nt increments.

  message(paste0("Connecting satalite positions up to ", sata.gap, " nt apart."))

  lapply(2:sata.gap, function(gap){
    g <<- connect_satalite_vertices(red.sites, g, gap)
    red.sites$clusID <<- clusters(g)$membership
  })

  message("Clusters after satalite connecting: ", clusters(g)$no)

  # During these steps of expanding the graph, it's possible that clusters grew
  # larger than they needed to be. This could be identified by seeing multiple
  # sources within a single cluster. This may arrize by reaching out and
  # connecting the edge of one smaller cluster to a larger one. If this occures,
  # the path between the sources will be 'snipped' between the nodes that share
  # the largest distance appart from one another.

  g <- break_connecting_source_paths(red.sites, g)
  red.sites$clusID <- clusters(g)$membership

  message(paste0("Clusters after clipping: ", clusters(g)$no))

  # In the end, sources that are within the range of satalites (5L default),
  # should be grouped together. These sources should be connected by an edge,
  # pointing toward the larger source. The chosen source will have the highest
  # abundance of widths / sonic breaks / fragment lengths.

  message("Connecting clusters with source nodes within ", sata.gap, " nt.")

  g <- connect_adjacent_clusters(red.sites, g, gap = sata.gap)
  red.sites$clusID <- clusters(g)$membership

  message(paste0("Final cluster count: ", clusters(g)$no))

  # If wide clusters are generated, these can confound the cluster source. For
  # this reason, the "true" source for the cluster will be resolved by
  # picking the source with the greatest number of fragment lengths, where ties
  # are decided by randon picking.

  g <- resolve_cluster_sources(red.sites, g)
  red.sites$clusID <- clusters(g)$membership

  # As cluster membership has been determined, the remaining section of the
  # function serves to consolidate the information and correct the original
  # sites object with the standardized positions.

  src.nodes <- sources(g)

  clus.data <- data.frame(
    "clusID" = 1:clusters(g)$no,
    "chr" = seqnames(red.sites[src.nodes]),
    "strand" = strand(red.sites[src.nodes]),
    "position" = start(red.sites[src.nodes]),
    "width" = width(unlist(range(
      GenomicRanges::split(red.sites, red.sites$clusID))))
  )

  sites <- sites[unlist(as.list(red.sites$revmap))]
  sites$clusID <- as.numeric(Rle(
    clusters(g)$membership,
    sapply(red.sites$revmap, length)))
  sites$called.pos <- ifelse(strand(sites) == "+", start(sites), end(sites))
  sites$adj.pos <- clus.data[match(sites$clusID, clus.data$clusID), "position"]

  message(paste0("Cumulative displacement per range: ",
                 sum(abs(sites$called.pos - sites$adj.pos))/length(sites)))

  ranges(sites) <- IRanges(
    start = ifelse(strand(sites) == "+", sites$adj.pos, start(sites)),
    end = ifelse(strand(sites) == "+", end(sites), sites$adj.pos)
  )

  sites
}
