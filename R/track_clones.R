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

track_clones <- function(sites.list, gap = 5L, track.origin = TRUE){
  if(class(sites.list) == "list"){grl.sites <- GRangesList(sites.list)}
  grl.sites <- sites.list

  if(track.origin){
    if(is.null(names(grl.sites))) names(grl.sites) <- 1:length(grl.sites)
    grl.sites <- GRangesList(lapply(1:length(grl.sites), function(i){
      sites <- sites.list[[i]]
      sites$origin <- rep(names(grl.sites[i]), length(sites))
      sites
    }))
  }

  ovlp.grps <- findOverlaps(
    flank(grl.sites, width = -1, start = TRUE),
    maxgap = gap,
    ignoreSelf = TRUE,
    ignoreRedundant = TRUE
  )

  if(length(ovlp.grps) > 0){
    ovlp.sites <- GenomicRanges::unlist(GRangesList(lapply(
      1:length(ovlp.grps),
      function(i){
        query <- grl.sites[[queryHits(ovlp.grps[i])]]
        subject <- grl.sites[[subjectHits(ovlp.grps[i])]]
        hits <- findOverlaps(
          flank(query, -1, start = TRUE),
          flank(subject, -1, start = TRUE),
          maxgap = gap)
        if(length(hits) > 0){
          sites <- c(query[queryHits(hits)], subject[subjectHits(hits)])
        }else{
          sites <- GRanges()
        }
        sites
    })))
  }else{
    message("No overlaping sites found between groups.")
  }

  if(length(ovlp.grps) > 0){
    ovlp.sites$posid <- generate_posid(ovlp.sites)
    sites.dfr <- distinct(
      GenomicRanges::as.data.frame(ovlp.sites, row.names = NULL))
    ranges <- IRanges(start = sites.dfr$start, end = sites.dfr$end)
    sites.gr <- GRanges(
      seqnames = sites.dfr$seqnames,
      ranges = ranges,
      strand = sites.dfr$strand
    )
    mcols(sites.gr) <- sites.dfr[,c(6:length(sites.dfr))]
    ovlp.list <- split(sites.gr, sites.gr$posid)
  }else{
    ovlp.list <- GRangesList()
  }
  message(paste("Number of overlaping sites:", length(ovlp.list)))
  ovlp.list
}
