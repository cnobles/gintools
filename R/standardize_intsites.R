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
#' standardize_intsites(sites, counts = "counts", min.gap = 1L, sata.gap = 5L)
#'
#' @param sites GRanges Object integration sites to standardize.
#' @param counts character string. Name of column holding range count
#' information (ie. read counts). If not supplied, assume count of 1 for each
#' row in sites.
#' @param min.gap integer minimum gap to consider combine integration sites.
#' @param sata.gap integer maximum distance to consider combining integration
#' sites.
#'
#' @examples
#' gr <- gintools:::generate_test_granges()
#'
#' standardize_intsites(gr)
#'
#' @author Christopher Nobles, Ph.D.
#' @export

standardize_intsites <- function(sites, counts = NULL, min.gap = 1L, sata.gap = 5L){

  # Retain original order
  sites$ori.order <- seq_along(sites)

  # Identify counts or abundance info, assume 1 if not given, but check for
  # a column named counts first and use, otherwise error out if column is not
  # found.

  if(!is.null(counts)){
    if(counts %in% names(GenomicRanges::mcols(sites))){
      counts_pos <- grep(counts, names(GenomicRanges::mcols(sites)))
      sites$func.counts <- GenomicRanges::mcols(sites)[,counts_pos]
    }else{
      stop("Could not identify 'counts' column.")
    }
  }else{
    message("Assuming abundance of 1 for each row of sites object.")
    sites$func.counts <- rep(1, length(sites))
  }

  # Start by reducing the sites object down to only the "intSite" positions
  # and storing the revmap for going back to the original sites.

  message(paste0("Standardizing ", length(sites), " positions."))
  message("Generating initial graph by connecting positions with 1 nt difference.")

  red.sites <- GenomicRanges::reduce(
    GenomicRanges::flank(sites, -1, start = TRUE),
    min.gapwidth = 0L,
    with.revmap = TRUE)
  red.sites$siteID <- seq_along(red.sites)
  df.counts <- GenomicRanges::mcols(sites) %>%
    as.data.frame(row.names = NULL)
  df.counts <- df.counts[unlist(red.sites$revmap),] %>%
    dplyr::mutate(grp = as.numeric(S4Vectors::Rle(
      values = seq_along(red.sites),
      lengths = S4Vectors::lengths(red.sites$revmap)))) %>%
    dplyr::group_by(grp) %>%
    dplyr::summarise(abund = sum(func.counts)) %>%
    dplyr::ungroup()

  red.sites$abund <- df.counts$abund

  red.hits <- GenomicRanges::as.data.frame(
    GenomicRanges::findOverlaps(red.sites, maxgap = min.gap, drop.self = TRUE))

  red.hits <- red.hits %>%
    dplyr::mutate(
      q_pos = GenomicRanges::start(red.sites[queryHits]),
      s_pos = GenomicRanges::start(red.sites[subjectHits]),
      q_abund = red.sites[queryHits]$abund,
      s_abund = red.sites[subjectHits]$abund,
      strand = as.vector(GenomicRanges::strand(red.sites[queryHits])),
      is.upstream = ifelse(strand == "+", q_pos < s_pos, q_pos > s_pos),
      keep = q_abund > s_abund,
      keep = ifelse(q_abund == s_abund, is.upstream, keep)) %>%
    dplyr::filter(keep)

  # Construct initial graph, the graph will be modified during other parts of
  # the function and can be used to identify clusters as well as 'sources'
  # within those clusters as it's a directed graph. Sources are identified and
  # used to identify the original position given a distribution.
  g <- igraph::make_empty_graph(n = length(red.sites), directed = TRUE) %>%
    igraph::add_edges(vzip(red.hits$queryHits, red.hits$subjectHits))

  red.sites$clusID <- igraph::clusters(g)$membership

  message(paste0("Initial cluster count: ", igraph::clusters(g)$no))

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

  null <- lapply(2:sata.gap, function(gap){
    g <<- connect_satalite_vertices(red.sites, g, gap, "upstream")
    red.sites$clusID <<- igraph::clusters(g)$membership
  })

  message("Clusters after satalite connecting: ", igraph::clusters(g)$no)

  # During these steps of expanding the graph, it's possible that clusters grew
  # larger than they needed to be. This could be identified by seeing multiple
  # sources within a single cluster. This may arrize by reaching out and
  # connecting the edge of one smaller cluster to a larger one. If this occures,
  # the path between the sources will be 'snipped' between the nodes that share
  # the largest distance appart from one another.

  g <- break_connecting_source_paths(red.sites, g, "upstream")
  red.sites$clusID <- igraph::clusters(g)$membership

  message(paste0("Clusters after clipping: ", igraph::clusters(g)$no))

  # In the end, sources that are within the range of satalites (5L default),
  # should be grouped together. These sources should be connected by an edge,
  # pointing toward the larger source. The chosen source will have the highest
  # abundance of widths / sonic breaks / fragment lengths.

  message("Connecting clusters with source nodes within ", sata.gap, " nt.")

  g <- connect_adjacent_clusters(red.sites, g, gap = sata.gap, "upstream")
  red.sites$clusID <- igraph::clusters(g)$membership

  message(paste0("Final cluster count: ", igraph::clusters(g)$no))

  # If wide clusters are generated, these can confound the cluster source. For
  # this reason, the "true" source for the cluster will be resolved by
  # picking the source with the greatest number of fragment lengths, where ties
  # are decided by randon picking.

  g <- resolve_cluster_sources(red.sites, g, "upstream")
  red.sites$clusID <- igraph::clusters(g)$membership

  # As cluster membership has been determined, the remaining section of the
  # function serves to consolidate the information and correct the original
  # sites object with the standardized positions.

  src.nodes <- sources(g)

  clus.data <- data.frame(
    "clusID" = seq_along(igraph::clusters(g)$csize),
    "chr" = GenomicRanges::seqnames(red.sites[src.nodes]),
    "strand" = GenomicRanges::strand(red.sites[src.nodes]),
    "position" = GenomicRanges::start(red.sites[src.nodes]),
    "width" = GenomicRanges::width(unlist(range(
      GenomicRanges::split(red.sites, red.sites$clusID))))
  )

  sites <- sites[unlist(IRanges::as.list(red.sites$revmap))]
  sites$clusID <- as.numeric(S4Vectors::Rle(
    igraph::clusters(g)$membership,
    S4Vectors::lengths(red.sites$revmap)))
  sites$called.pos <- ifelse(
    GenomicRanges::strand(sites) == "+",
    GenomicRanges::start(sites),
    GenomicRanges::end(sites))
  sites$adj.pos <- clus.data[match(sites$clusID, clus.data$clusID), "position"]

  message(paste0("Cumulative displacement per range: ",
                 sum(abs(sites$called.pos - sites$adj.pos))/length(sites)))

  ranges(sites) <- IRanges::IRanges(
    start = ifelse(
      GenomicRanges::strand(sites) == "+",
      sites$adj.pos,
      GenomicRanges::start(sites)),
    end = ifelse(
      GenomicRanges::strand(sites) == "+",
      GenomicRanges::end(sites),
      sites$adj.pos)
  )

  sites <- sites[order(sites$ori.order)]
  sites$ori.order <- sites$func.counts <- NULL
  sites
}
