#' Condense integration site GRanges to distinct genomic positions.
#'
#' \code{condense_intsites} condenses or collapses GRanges data to distinct
#' integration positions in the reference genome.
#'
#' @description Given a GRanges object representing the genomic ranges of
#' sequenced DNA flanking the viral/vector integrating site, this function
#' returns only the distinct integration site positions. Further, using the
#' options for return.abundance and method, the function can return abundance
#' or estimated abundance information based on fragment length (see sonicLength
#' package). Initial GRanges object can carry multiple samples which can be
#' separated and analyzed independently by using the grouping option.
#'
#' @usage
#' condense_intsites(sites_to_condense)
#'
#' condense_intsites(sites_to_condense, grouping = NULL,
#' return.abundance = FALSE, method = "fragLen", replicates = "replicates")
#'
#' @param sites_to_condense a GRanges object where each row represents a genomic
#' range of reference genome DNA.
#'
#' @param grouping metadata column name (input as character) which designates
#' which genomic ranges belong to which sample.
#'
#' @param return.abundance logical that will subject the data to abundance
#' calculations, default is FALSE.
#'
#' @param method either "fragLen" (default) or "estAbund" to be passed to
#' determine_abundance function. method = "fragLen" determines the abundance by
#' the number of unique genomic ranges widths while the method = "estAbund"
#' utilizes the sonicLength package.
#'
#' @param replicates character name of the column carrying the replicate
#' information, relevant for abundance calculations and irrelent after using the
#' this function.
#'
#' @examples
#' gr <- gintools:::generate_test_granges()
#' std.gr <- standardize_intsites(gr)
#' condense_intsites(std.gr)
#'
#' condense_intsites(std.gr, return.abundance = TRUE)
#'
#' @author Christopher Nobles, Ph.D.
#' @export

condense_intsites <- function(sites_to_condense, grouping = NULL,
                              return.abundance = FALSE, method = "fragLen",
                              replicates = "replicates"){
  mcols <- GenomicRanges::mcols(sites_to_condense)

  grp <- which(grouping == names(mcols))
  if(length(grp) != 0){
    group <- as.character(mcols[,grp])
  }else{
    group <- rep("group1", length(sites_to_condense))
  }

  sites_to_condense$posID <- generate_posid(sites_to_condense)
  sites_to_condense$groupID <- paste0(group, "^", sites_to_condense$posID)

  groupIDs <- S4Vectors::unique(sites_to_condense$groupID)
  first.hits <- match(groupIDs, sites_to_condense$groupID)

  condensed.gr <- GenomicRanges::flank(sites_to_condense, -1, start = TRUE)
  condensed.gr <- condensed.gr[first.hits]

  if(return.abundance){
    abund.dfr <- determine_abundance(
      sites_to_condense, grouping = grouping,
      replicates = replicates, method = method)
    abund.dfr$groupID <- paste0(abund.dfr$group, "^", abund.dfr$posID)
    mcols <- GenomicRanges::mcols(condensed.gr)
    all.cols <- merge(mcols, abund.dfr, by = "groupID")
    order <- match(condensed.gr$groupID, all.cols$groupID)
    all.cols <- all.cols[order,]
    GenomicRanges::mcols(condensed.gr) <- all.cols
    condensed.gr$posID.x <- NULL
    condensed.gr$posID.y <- NULL
    condensed.gr$posID <- generate_posid(condensed.gr)
  }
  condensed.gr$groupID <- NULL
  condensed.gr$group <- NULL
  names(condensed.gr) <- NULL
  message("Check to make sure all metadata columns are still relavent.")
  condensed.gr
}
