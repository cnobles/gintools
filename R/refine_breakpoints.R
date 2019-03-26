#' Refine or resolve sonic break points within a dataset
#'
#' \code{refine_breakpoints} returns a GRanges object where the sonic break
#' point positions have been adjusted based on the dataset and counts.
#'
#' @description Given a GRanges object ....
#'
#' @usage
#' refine_breakpoints(input.sites)
#'
#' refine_breakpoints(input.sites, count = "counts", min.gap = 1L, sata.gap = 3L)
#'
#' @param input.sites GRanges object. Ranges of alignments or sites to adjust.
#' @param counts.col character string. Name of column holding range count
#' information (ie. read counts). If not supplied, assume count of 1 for each
#' row in sites.
#' @param min.gap integer minimum gap to consider combining break points.
#' @param sata.gap integer maximum distance to consider combining break points.
#' @param details logical If TRUE, data columns will be appended to the metadata
#' of the GRanges object noting the original and adjusted position for the site.
#' FALSE by default.
#' @param quiet logical during processing, should messages be suppressed which
#' report findings? TRUE by default.
#'
#' @examples
#' gr <- gintools:::generate_test_granges()
#'
#' refine_breakpoints(gr)
#'
#' @author Christopher L. Nobles, Ph.D.
#'
#' @importFrom magrittr %>%
#'
#' @export
#'

refine_breakpoints <- function(input.sites, counts.col = NULL, min.gap = 1L,
                               sata.gap = 3L, details = FALSE, quiet = TRUE){

  stopifnot(class(input.sites) == "GRanges")
  sites <- GenomicRanges::granges(input.sites)

  if(!quiet){
    message(paste0("Refining ", length(sites), " break points."))
    message(
      "Generating initial graph by connecting positions with 1 nt difference.")
  }

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
    sites$counts <- rep(1, length(sites))
  }

  # Reduce the genomic locations of break points down to only unique positions,
  # and identify the abundance of the positions

  red_sites <- GenomicRanges::reduce(
    GenomicRanges::flank(sites, -1, start = FALSE),
    min.gapwidth = 0L,
    with.revmap = TRUE)
  red_sites$breakpoint.id <- seq_along(red_sites)

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
      is.downstream = ifelse(strand == "+", pos.q > pos.s, pos.q < pos.s),
      keep = abund.q > abund.s,
      keep = ifelse(abund.q == abund.s, is.downstream, keep)) %>%
    dplyr::filter(keep)

  # Construct initial graph, the graph will be modified during other parts of
  # the function and can be used to identify clusters as well as 'sources'
  # within those clusters as it's a directed graph. Sources are identified and
  # used to identify the original position given a distribution.

  g <- igraph::make_empty_graph(n = length(red_sites), directed = TRUE) %>%
    igraph::add_edges(with(red_hits, vzip(hits.q, hits.s)))

  red_sites$clus.id <- igraph::clusters(g)$membership

  if(!quiet) message(paste0("Initial cluster count: ", igraph::clusters(g)$no))

  # Identify satalite positions that should be included in clusters up to the
  # sata.gap max. This portion of the function tries to reach out to up to the
  # sata.gap distance from the boundry of a cluster to see if there are any
  # further "satalite" positions that have been annotated. It does this by
  # iterively increasing the size from 2nt to the sata.gap by 1nt increments.

  if(!quiet){
    message(paste0(
      "Connecting satalite positions up to ", sata.gap, " nt apart."))
  }

  null <- lapply(2:sata.gap, function(gap){
    g <<- connect_satalite_vertices(red_sites, g, gap, "downstream")
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

  g <- break_connecting_source_paths(red_sites, g, "downstream")
  red_sites$clus.id <- igraph::clusters(g)$membership

  if(!quiet){
    message(paste0("Final break point cluster count: ", igraph::clusters(g)$no))
  }

  # As cluster membership has been determined, the remaining section of the
  # function serves to consolidate the information and correct the original
  # sites object with the adjusted breakpoints.

  src_nodes <- sources(g)

  clus_data <- data.frame(
    "clus.id" = seq_along(igraph::clusters(g)$csize),
    "chr" = GenomicRanges::seqnames(red_sites[src_nodes]),
    "strand" = GenomicRanges::strand(red_sites[src_nodes]),
    "breakpoint" = GenomicRanges::start(red_sites[src_nodes]),
    "width" = GenomicRanges::width(unlist(range(
      GenomicRanges::split(red_sites, red_sites$clus.id))))
  )

  sites <- sites[unlist(IRanges::as.list(red_sites$revmap))]
  sites$clus.id <- as.numeric(S4Vectors::Rle(
    igraph::clusters(g)$membership,
    S4Vectors::lengths(red_sites$revmap)))
  sites$called.bp <- ifelse(
    GenomicRanges::strand(sites) == "+",
    GenomicRanges::end(sites),
    GenomicRanges::start(sites))
  sites$adj.bp <- clus_data[
    match(sites$clus.id, clus_data$clus.id), "breakpoint"]

  if(!quiet){
    message(paste0(
      "Cumulative displacement per range: ",
      round(
        sum( abs(sites$called.bp - sites$adj.bp) ) / length(sites),
        digits = 1)))
  }

  GenomicRanges::ranges(sites) <- IRanges::IRanges(
    start = ifelse(
      GenomicRanges::strand(sites) == "+",
      GenomicRanges::start(sites),
      sites$adj.bp),
    end = ifelse(
      GenomicRanges::strand(sites) == "+",
      sites$adj.bp,
      GenomicRanges::end(sites))
  )

  # Transfer range information from sites to output_sites object
  sites <- sites[order(sites$ori.order)]
  output_sites <- input.sites
  GenomicRanges::ranges(output_sites) <- GenomicRanges::ranges(sites)

  if(details){
    GenomicRanges::mcols(output_sites)$called.pos <- GenomicRanges::mcols(
      sites)$called.bp
    GenomicRanges::mcols(output_sites)$adj.pos <- GenomicRanges::mcols(
      sites)$adj.bp
  }

  output_sites
}
