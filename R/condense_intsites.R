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
#' condense_intsites(sites.to.condense)
#'
#' condense_intsites(sites.to.condense, grouping = NULL,
#' return.abundance = FALSE, method = "fragLen", replicates = "replicates")
#'
#' @param sites.to.condense a GRanges object where each row represents a genomic
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
#' @param quiet logical indicating if output message should be displayed.
#'
#' @examples
#' gr <- gintools:::generate_test_granges()
#' std.gr <- standardize_sites(gr)
#' condense_intsites(std.gr)
#'
#' condense_intsites(std.gr, return.abundance = TRUE)
#'
#' @author Christopher Nobles, Ph.D.
#' @export

condense_intsites <- function(sites.to.condense, grouping = NULL,
                              return.abundance = FALSE, method = "fragLen",
                              replicates = "replicates", quiet = TRUE){
  
  mcols <- GenomicRanges::mcols(sites.to.condense)

  grp <- which(grouping == names(mcols))
  
  if( length(grp) != 0 ){
    group <- as.character(mcols[,grp])
  }else{
    group <- rep("group1", length(sites.to.condense))
  }

  sites.to.condense$posid <- generate_posid(sites.to.condense)
  sites.to.condense$group.id <- paste0(group, sites.to.condense$posid)

  group_ids <- S4Vectors::unique(sites.to.condense$group.id)
  first_hits <- match(group_ids, sites.to.condense$group.id)

  condensed_gr <- GenomicRanges::flank(sites.to.condense, -1, start = TRUE)
  condensed_gr <- condensed_gr[first_hits]

  if( return.abundance ){
    
    abund_df <- determine_abundance(
      sites.to.condense, grouping = grouping,
      replicates = replicates, method = method
    )
    
    abund_df$group.id <- paste0(abund_df$group, abund_df$posid)
    mcols <- GenomicRanges::mcols(condensed_gr)
    all_cols <- merge(mcols, abund_df, by = "group.id")
    order <- match(condensed_gr$group.id, all_cols$group.id)
    all_cols <- all_cols[order,]
    GenomicRanges::mcols(condensed_gr) <- all_cols
    condensed_gr$posid.x <- NULL
    condensed_gr$posid.y <- NULL
    condensed_gr$posid <- generate_posid(condensed_gr)
    
  }
  
  condensed_gr$group.id <- NULL
  condensed_gr$group <- NULL
  names(condensed_gr) <- NULL
  
  if( !quiet ){
    message("Check to make sure all metadata columns are still relavent.")
  }
  
  condensed_gr
  
}
