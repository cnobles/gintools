#' Standardize genomic site positions within a dataset
#'
#' \code{standardize_sites} returns a GRanges object where the site positions
#' have been standardized with sites within the gap distance.
#'
#' @description Given a GRanges object ...
#'
#' @usage
#' standardize_sites(input.sites)
#'
#' standardize_sites(input.sites, counts.col = "counts", min.gap = 1L, sata.gap = 5L)
#'
#' @param input.sites GRanges object of sites to standardize.
#' @param counts.col character string. Name of column holding range count
#' information (ie. read counts). If not supplied, assume count of 1 for each
#' row in sites.
#' @param min.gap integer minimum gap to consider combine sites.
#' @param sata.gap integer maximum distance to consider combining sites.
#' @param details logical If TRUE, data columns will be appended to the metadata
#' of the GRanges object noting the original and adjusted position for the site.
#' FALSE by default.
#' @param quiet logical during processing, should messages be suppressed which
#' report findings? TRUE by default.
#'
#' @examples
#' gr <- gintools:::generate_test_granges()
#'
#' standardize_sites(gr)
#'
#' @author Christopher Nobles, Ph.D.
#' @export
#'

standardize_sites <- function(input.sites, counts.col = NULL, min.gap = 1L,
                              sata.gap = 5L, details = FALSE, quiet = TRUE){
  stopifnot(class(input.sites) == "GRanges")
  sites <- GenomicRanges::granges(input.sites)

  # Retain original order
  sites$ori.order <- seq_along(sites)

  # Identify counts or abundance info, assume 1 if not given, but check for
  # a column named counts first and use, otherwise error out if column is not
  # found.

  if(!is.null(counts.col)){
    if(counts.col %in% names(GenomicRanges::mcols(input.sites))){
      counts_pos <- grep(counts.col, names(GenomicRanges::mcols(input.sites)))
      sites$counts <- GenomicRanges::mcols(input.sites)[,counts_pos]
    }else{
      stop("Could not identify 'counts' column.")
    }
  }else{
    if(!quiet) message("Assuming abundance of 1 for each row of sites object.")
    sites$counts <- rep(1, length(input.sites))
  }

  # Start by reducing the sites object down to only the site's positions
  # and store the revmap for going back to the original sites object.
  if(!quiet){
    message(paste0("Standardizing ", length(sites), " positions."))
    message(
      "Generating initial graph by connecting positions with 1 nt difference.")
  }

  # Construct reduced sites object to initialize the processing
  red_sites <- GenomicRanges::reduce(
    GenomicRanges::flank(sites, -1, start = TRUE),
    min.gapwidth = 0L,
    with.revmap = TRUE)
  red_sites$site.id <- seq_along(red_sites)

  # Summarise count data for reduced site object
  red_counts <- GenomicRanges::mcols(sites) %>%
    as.data.frame(row.names = NULL)
  red_counts <- red_counts[unlist(red_sites$revmap),] %>%
    dplyr::mutate(grp = as.numeric(S4Vectors::Rle(
      values = seq_along(red_sites),
      lengths = S4Vectors::lengths(red_sites$revmap)))) %>%
    dplyr::group_by(grp) %>%
    dplyr::summarise(abund = sum(counts)) %>%
    dplyr::ungroup()

  red_sites$abund <- red_counts$abund

  # Identify which sites are adjacent to each other
  red_hits <- GenomicRanges::as.data.frame(
    GenomicRanges::findOverlaps(
      red_sites, maxgap = min.gap - 1L, drop.self = TRUE
    )
  )

  # Organize data and filter for constructing directed graph
  red_hits <- red_hits %>%
    dplyr::mutate(
      hits.q = queryHits,
      hits.s = subjectHits,
      pos.q = GenomicRanges::start(red_sites[hits.q]),
      pos.s = GenomicRanges::start(red_sites[hits.s]),
      abund.q = red_sites[hits.q]$abund,
      abund.s = red_sites[hits.s]$abund,
      strand = as.vector(GenomicRanges::strand(red_sites[hits.q])),
      is.upstream = ifelse(strand == "+", pos.q < pos.s, pos.q > pos.s),
      keep = abund.q > abund.s,
      keep = ifelse(abund.q == abund.s, is.upstream, keep)) %>%
    dplyr::filter(keep)

  # Construct initial graph, the graph will be modified during other parts of
  # the function and can be used to identify clusters as well as 'sources'
  # within those clusters as it's a directed graph. Sources are identified and
  # used to identify the original position given a distribution.
  g <- igraph::make_empty_graph(n = length(red_sites), directed = TRUE) %>%
    igraph::add_edges(with(red_hits, vzip(hits.q, hits.s)))

  red_sites$clus.id <- igraph::clusters(g)$membership

  if(!quiet) message(paste0("Initial cluster count: ", igraph::clusters(g)$no))

  # When clusters are formed of positions with equivalent abundance, both
  # directional edges are created. This type of cluster does not have sources
  # or sinks. The following set resolves which directed edges to remove, but
  # imposes an upstream bias.

#  g <- resolve_source_deficiency(red_sites, g, bias = "upstream")
#  red_sites$clus.id <- clusters(g)$membership

  # Identify satalite positions that should be included in clusters up to 5nt
  # away. This portion of the function tries to reach out to up to 5nt from the
  # boundry of a cluster to see if there are any further "satalite" positions
  # that have been annotated. It does this by iterively increasing the size
  # from 2nt to 5nt by 1nt increments.

  if(!quiet){
    message(paste0(
      "Connecting satalite positions up to ", sata.gap, " nt apart."))
  }

  null <- lapply(2:sata.gap, function(gap){
    g <<- connect_satalite_vertices(red_sites, g, gap, "upstream")
    red_sites$clus.id <<- igraph::clusters(g)$membership
  })

  if(!quiet){
    message("Clusters after satalite connecting: ", igraph::clusters(g)$no)
  }

  # During these steps of expanding the graph, it's possible that clusters grew
  # larger than they needed to be. This could be identified by seeing multiple
  # sources within a single cluster. This may arrize by reaching out and
  # connecting the edge of one smaller cluster to a larger one. If this occures,
  # the path between the sources will be 'snipped' between the nodes that share
  # the largest distance appart from one another.

  g <- break_connecting_source_paths(red_sites, g, "upstream")
  red_sites$clus.id <- igraph::clusters(g)$membership

  if(!quiet){
    message(paste0("Clusters after clipping: ", igraph::clusters(g)$no))
  }

  # In the end, sources that are within the range of satalites (5L default),
  # should be grouped together. These sources should be connected by an edge,
  # pointing toward the larger source. The chosen source will have the highest
  # abundance of widths / sonic breaks / fragment lengths.

  if(!quiet){
    message("Connecting clusters with source nodes within ", sata.gap, " nt.")
  }

  g <- connect_adjacent_clusters(red_sites, g, gap = sata.gap, "upstream")
  red_sites$clus.id <- igraph::clusters(g)$membership

  if(!quiet) message(paste0("Final cluster count: ", igraph::clusters(g)$no))

  # If wide clusters are generated, these can confound the cluster source. For
  # this reason, the "true" source for the cluster will be resolved by
  # picking the source with the greatest number of fragment lengths, where ties
  # are decided by randon picking.

  g <- resolve_cluster_sources(red_sites, g, "upstream")
  red_sites$clus.id <- igraph::clusters(g)$membership

  # As cluster membership has been determined, the remaining section of the
  # function serves to consolidate the information and adjust the original
  # sites object with the standardized positions.

  src.nodes <- sources(g)

  clus.data <- data.frame(
    "clus.id" = seq_along(igraph::clusters(g)$csize),
    "chr" = GenomicRanges::seqnames(red_sites[src.nodes]),
    "strand" = GenomicRanges::strand(red_sites[src.nodes]),
    "position" = GenomicRanges::start(red_sites[src.nodes]),
    "width" = GenomicRanges::width(unlist(range(
      GenomicRanges::split(red_sites, red_sites$clus.id))))
  )

  # Index the sites object and add adjusted data to metadata
  sites <- sites[unlist(IRanges::as.list(red_sites$revmap))]
  sites$clus.id <- as.numeric(S4Vectors::Rle(
    igraph::clusters(g)$membership,
    S4Vectors::lengths(red_sites$revmap)))
  sites$called.pos <- ifelse(
    GenomicRanges::strand(sites) == "+",
    GenomicRanges::start(sites),
    GenomicRanges::end(sites))
  sites$adj.pos <- clus.data[
    match(sites$clus.id, clus.data$clus.id), "position"]

  if(!quiet){
    message(paste0(
      "Cumulative displacement per range: ",
      sum(abs(sites$called.pos - sites$adj.pos))/length(sites)))
  }

  # Adjust range information within sites object
  GenomicRanges::ranges(sites) <- IRanges::IRanges(
    start = ifelse(
      GenomicRanges::strand(sites) == "+",
      sites$adj.pos,
      GenomicRanges::start(sites)),
    end = ifelse(
      GenomicRanges::strand(sites) == "+",
      GenomicRanges::end(sites),
      sites$adj.pos)
  )

  # Transfer range information from sites to output_sites object
  sites <- sites[order(sites$ori.order)]
  output_sites <- input.sites
  GenomicRanges::ranges(output_sites) <- GenomicRanges::ranges(sites)

  if(details){
    GenomicRanges::mcols(output_sites)$called.pos <- GenomicRanges::mcols(
      sites)$called.pos
    GenomicRanges::mcols(output_sites)$adj.pos <- GenomicRanges::mcols(
      sites)$adj.pos
  }

  output_sites
}
