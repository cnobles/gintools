#' Refine or resolve sonic break points within a dataset
#'
#' \code{refine_breakpoints} returns a GRanges object where the sonic break
#' point positions have been adjusted based on the dataset and counts.
#'
#' @description Given a GRanges object ....
#'
#' @usage
#' refine_breakpoints(sites)
#'
#' refine_breakpoints(sites, count = "counts", min.gap = 1L, sata.gap = 3L)
#'
#' @param sites GRanges object. Integration site ranges to adjust.
#' @param counts character string. Name of column holding range count
#' information (ie. read counts). If not supplied, assume count of 1 for each
#' row in sites.
#' @param min.gap integer minimum gap to consider combining break points.
#' @param sata.gap integer maximum distance to consider combining break points.
#'
#' @examples
#' gr <- .generate_test_granges()
#'
#' refine_breakpoints(gr)
#'
#' @author Christopher L. Nobles, Ph.D.
#' @export

refine_breakpoints <- function(sites, counts = NULL, min.gap = 1L, sata.gap = 3L){
  require(gintools)

  message(paste0("Refining ", length(sites), " break points."))
  message("Generating initial graph by connecting positions with 1 nt difference.")

  # Retain original order
  sites$ori.order <- 1:length(sites)

  # Identify counts or abundance info, assume 1 if not given, but check for
  # a column named counts first and use, otherwise error out if column is not
  # found.

  if(!is.null(counts)){
    if(counts %in% names(mcols(sites))){
      counts_pos <- grep(counts, names(mcols(sites)))
      sites$func.counts <- mcols(sites)[,counts_pos]
    }else{
      stop("Could not identify 'counts' column.")
    }
  }else{
    message("Assuming abundance of 1 for each row of sites object.")
    sites$func.counts <- rep(1, length(sites))
  }

  # Reduce the genomic locations of break points down to only unique positions,
  # and identify the abundance of the positions

  red.sites <- reduce(
    flank(sites, -1, start = FALSE),
    min.gapwidth = 0L,
    with.revmap = TRUE)
  red.sites$breakpointID <- seq(1:length(red.sites))
  revmap <- as.list(red.sites$revmap)
  red.sites$abundance <- sapply(revmap, function(x){
    sum(sites[x]$func.counts)
  })
  red.hits <- GenomicRanges::as.data.frame(
    findOverlaps(red.sites, maxgap = min.gap, ignoreSelf = TRUE))

  red.hits <- red.hits %>%
    dplyr::mutate(q_pos = start(red.sites[queryHits])) %>%
    dplyr::mutate(s_pos = start(red.sites[subjectHits])) %>%
    dplyr::mutate(q_abund = red.sites[queryHits]$abundance) %>%
    dplyr::mutate(s_abund = red.sites[subjectHits]$abundance) %>%
    dplyr::mutate(strand = as.vector(strand(red.sites[queryHits]))) %>%
    dplyr::mutate(is.downstream = ifelse(
      strand == "+",
      q_pos > s_pos,
      q_pos < s_pos)) %>%
    dplyr::mutate(keep = q_abund > s_abund) %>%
    dplyr::mutate(keep = ifelse(
      q_abund == s_abund,
      is.downstream,
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

  # Identify satalite positions that should be included in clusters up to the
  # sata.gap max. This portion of the function tries to reach out to up to the
  # sata.gap distance from the boundry of a cluster to see if there are any
  # further "satalite" positions that have been annotated. It does this by
  # iterively increasing the size from 2nt to the sata.gap by 1nt increments.

  message(paste0("Connecting satalite positions up to ", sata.gap, " nt apart."))

  lapply(2:sata.gap, function(gap){
    g <<- connect_satalite_vertices(red.sites, g, gap, "downstream")
    red.sites$clusID <<- clusters(g)$membership
  })

  message("Clusters after satalite connecting: ", clusters(g)$no)

  # During these steps of expanding the graph, it's possible that clusters grew
  # larger than they needed to be. This could be identified by seeing multiple
  # sources within a single cluster. This may arrize by reaching out and
  # connecting the edge of one smaller cluster to a larger one. If this occures,
  # the path between the sources will be 'snipped' between the nodes that share
  # the largest distance appart from one another.

  g <- break_connecting_source_paths(red.sites, g, "downstream")
  red.sites$clusID <- clusters(g)$membership

  message(paste0("Final break point cluster count: ", clusters(g)$no))

  # As cluster membership has been determined, the remaining section of the
  # function serves to consolidate the information and correct the original
  # sites object with the adjusted breakpoints.

  src.nodes <- sources(g)

  clus.data <- data.frame(
    "clusID" = 1:clusters(g)$no,
    "chr" = seqnames(red.sites[src.nodes]),
    "strand" = strand(red.sites[src.nodes]),
    "breakpoint" = start(red.sites[src.nodes]),
    "width" = width(unlist(range(
      GenomicRanges::split(red.sites, red.sites$clusID))))
  )

  sites <- sites[unlist(as.list(red.sites$revmap))]
  sites$clusID <- as.numeric(Rle(
    clusters(g)$membership,
    sapply(red.sites$revmap, length)))
  sites$called.bp <- ifelse(strand(sites) == "+", end(sites), start(sites))
  sites$adj.bp <- clus.data[match(sites$clusID, clus.data$clusID), "breakpoint"]

  message(paste0("Cumulative displacement per range: ",
                 round(
                   sum(abs(sites$called.bp - sites$adj.bp))/length(sites),
                   digits = 1)))

  ranges(sites) <- IRanges(
    start = ifelse(strand(sites) == "+", start(sites), sites$adj.bp),
    end = ifelse(strand(sites) == "+", sites$adj.bp, end(sites))
  )

  sites <- sites[order(sites$ori.order)]
  sites$ori.order <- sites$func.counts <- sites$clusID <- NULL
  sites
}
