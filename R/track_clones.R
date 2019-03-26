#' Identifies or tracks integration sites through multiple GRanges objects.
#'
#' \code{track_clones} returns a GRangesList of integration sites shared between
#' multiple GRanges objects.
#'
#' @description Given a GRangesList of integration site sets, this function will
#' return a GRangesList (separated by position ID) of all sites shared within
#' the gap distance between the sets. Using the track.origin option
#' (default = TRUE) will append a column to the metadata of the GRanges objects
#' designating which set they originated in. If the GRangesList has names, these
#' designations will be the names of the list items. Note: as the output is a
#' GRangesList separated by position ID, some list objects may appear to come
#' only from one origin set, but they will overlap within the gap distance to
#' another list object from a different origin.
#'
#' @usage
#' track_clones(sites.list)
#'
#' track_clones(sites.list, gap = 5L, track.origin = TRUE)
#'
#' @param sites.list a GRangesList object containing sets of integration sites
#' or ranges (one per row).
#'
#' @param gap an integer designating the nucleotide distance or window for which
#' to group integration sites.
#'
#' @param track.origin logical for whether to append information regarding the
#' origin of the integration site. Use if there is no sample information in the
#' metadata columns of the input GRanges objects.
#'
#' @param quiet logical for whether or not to message to the terminal processing
#' findings. True by default.
#'
#' @examples
#' gr <- GRanges(
#'   seqnames = rep("chr1", 100),
#'   ranges = IRanges(
#'     start = 1:100,
#'     width = sample(30:100, 100, replace = TRUE)),
#'   strand = rep(c("+", "-"), 50))
#'
#' gr1 <- sample(gr, 50)
#' gr2 <- sample(gr, 50)
#'
#' grl <- GRangesList(gr1, gr2)
#'
#' track_clones(grl, gap = 0L, track.origin = TRUE)
#'
#' @author Christopher Nobles, Ph.D.
#' @export

track_clones <- function(sites.list, gap = 5L, track.origin = TRUE,
                         quiet = TRUE){
  if(class(sites.list) == "list"){
    grl_sites <- GenomicRanges::GRangesList(sites.list)
  }else if(class(sites.list) == "GRangesList"){
    grl_sites <- sites.list
  }else{
    stop(
      "Input sites.list must be either a ",
      "list of GRanges or GRangesList object.")
  }

  # Track origin of matched site
  if(track.origin){
    if(is.null(names(grl_sites))) names(grl_sites) <- seq_along(grl_sites)
    grl_sites <- GenomicRanges::GRangesList(lapply(
      seq_along(grl_sites), function(i){
        sites <- sites.list[[i]]
        sites$origin <- rep(names(grl_sites[i]), length(sites))
        sites
    }))
  }

  # Identify shared sites across listed objects
  ovlp_grps <- GenomicRanges::findOverlaps(
    GenomicRanges::flank(grl_sites, width = -1, start = TRUE),
    maxgap = gap - 1L,
    drop.self = TRUE,
    drop.redundant = TRUE
  )

  # Identify which specific sites are found in multiple objects of list
  if(length(ovlp_grps) > 0){
    ovlp_sites <- unlist(GenomicRanges::GRangesList(lapply(
      seq_along(ovlp_grps),
      function(i){
        query <- grl_sites[[S4Vectors::queryHits(ovlp_grps[i])]]
        subject <- grl_sites[[S4Vectors::subjectHits(ovlp_grps[i])]]
        hits <- GenomicRanges::findOverlaps(
          GenomicRanges::flank(query, -1, start = TRUE),
          GenomicRanges::flank(subject, -1, start = TRUE),
          maxgap = gap - 1L)
        if(length(hits) > 0){
          sites <- c(
            query[S4Vectors::queryHits(hits)],
            subject[S4Vectors::subjectHits(hits)])
        }else{
          sites <- GenomicRanges::GRanges()
        }
        sites
    })))
  }else{
    message("No overlaping sites found between groups.")
  }

  # Output grl of sites tracked across samples
  if(length(ovlp_sites) > 0){
    sites <- unique_granges(ovlp_sites)
    ovlp_list <- GenomicRanges::split(sites, generate_posid(sites))
  }else{
    ovlp_list <- GenomicRanges::GRangesList()
  }
  if(!quiet) message(paste("Number of overlaping sites:", length(ovlp_list)))
  ovlp_list
}
