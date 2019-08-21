#' Determine the abundance of integration sites within a GRanges object.
#'
#' \code{determine_abundance} returns the abundance of all distinct integration
#' sites within the given GRanges object based on the unique fragment length
#' (width) or utilizes the sonicLength package.
#'
#' @description Given a GRanges object containing the genomic ranges of
#' flanking genomic DNA from integration sites, the function returns abundance
#' information based on unique widths for each site (method = "fragLen",
#' default) or applies a more powerful tool developed in sonicLength package
#' (method = "estAbund"). If multiple samples are given, they can be grouped
#' separately by giving the metadata column to the grouping option.
#'
#' @usage
#' determine_abundance(sites)
#'
#' determine_abundance(sites, grouping = NULL, replicates = NULL,
#' method = "fragLen")
#'
#' @param sites GRanges object with flanking genomic DNA ranges.
#'
#' @param grouping Character name of the metadata column with grouping
#' information.
#'
#' @param replicates Character name of the metadata column which dictates
#' replicate information.
#'
#' @param method Character either "fragLen" for determining abundance by unique
#' fragment length or "estAbund" to use the sonicLength package.
#'
#' @examples
#' gr <- gintools:::generate_test_granges()
#' std.gr <- standardize_sites(gr)
#'
#' determine_abundance(std.gr)
#'
#' @author Christopher Nobles, Ph.D.
#' @export

determine_abundance <- function(sites, grouping = NULL, replicates = NULL,
                                method = "fragLen"){
  
  sites$posid <- generate_posid(sites)
  in_mcols <- GenomicRanges::mcols(sites)
  GenomicRanges::mcols(sites) <- NULL

  grp <- which(grouping == names(in_mcols))
  
  if( length(grp) != 0 ){
    group <- as.character(in_mcols[,grp])
  }else{
    group <- rep("group1", length(sites))
  }

  reps <- which(replicates == names(in_mcols))
  
  if( length(reps) != 0 ){
    replicates <- as.character(in_mcols[,reps])
  }else{
    replicates <- rep("1", length(sites))
  }

  sites_df <- data.frame(
    "posid"= generate_posid(sites),
    "fragLen" = GenomicRanges::width(sites),
    "replicates" = replicates,
    "group" = group
  )

  if( method == "fragLen" ){
    
    abundCalc <- function(locations, fragLen, replicates, group){
      
      group <- unique(group)
      dfr <- data.frame(
        "loci.id" = locations,
        "fragLen" = fragLen,
        "replicates" = replicates)
      dfr_dist <- dplyr::distinct(dfr)
      sites_list <- split(dfr_dist, dfr_dist$loci.id)
      abundances <- sapply(sites_list, function(x) nrow(x) )
      abundances <- abundances[unique(dfr$loci.id)]

      abund_df <- data.frame(
        "posid" = names(abundances),
        "estAbund" = abundances,
        "group" = rep(group, length(abundances)))
      abund_df
    
    }
    
  }else if( method == "estAbund" ){
    
    if( !requireNamespace("sonicLength", quietly = TRUE) ){
      stop(
        "sonicLength package not installed. ",
        "Please install or change method for ",
        "abundance calculation to 'fragLen'.",
        call. = FALSE
      )
    }

    abundCalc <- function(locations, fragLen, replicates, group){
      
      group <- unique(group)
      
      if( length(unique(sites_df$replicates)) == 1 ){
        theta_list <- sonicLength::estAbund(
          locations=locations, lengths=fragLen)
      }else{
        theta_list <- sonicLength::estAbund(
          locations=locations, lengths=fragLen, replicates=replicates)
      }
      
      posid <- names(theta_list$theta)
      abundances <- theta_list$theta

      if( length(unique(sites_df$replicates)) == 1 ){
        abund_df <- data.frame("posid" = posid,"estAbund" = abundances)
      }else{
        abund_df <- data.frame(
          "posid" = posid,
          "estAbund" = abundances,
          "replicates" = length(unique(theta_list$data$replicates)))
      }
      abund_df$group <- group
      abund_df
      
    }
    
  }else{
    
    stop("Must choose either fragLen or estAbund for method.")
    
  }

  sites_list <- split(sites_df, sites_df$group)

  abund_list <- lapply(sites_list, function(x){
    abundCalc(
      locations = x$posid,
      fragLen = x$fragLen,
      replicates = x$replicates,
      group = x$group
    )
  })

  abund_list <- lapply(abund_list, function(x){
    
    x %>%
      dplyr::mutate(
        totalAbund = sum(x$estAbund),
        estAbund = round(estAbund),
        relAbund = estAbund/totalAbund,
        relRank = rank(-1*relAbund, ties.method = "max")
      )

  })

  abund_df <- dplyr::bind_rows(abund_list)
  abund_df$totalAbund <- NULL
  
  if( nrow(abund_df) > 0 ){
    abund_df <- abund_df[
      , c("posid", "group", "estAbund", "relAbund", "relRank")
    ]
  }else{
    abund_df <- data.frame(
      "posid" = character(), "group" = character(), "estAbund" = numeric(),
      "relAbund" = numeric(), "relRank" = integer(), stringsAsFactors = FALSE
    )
  }
  
  abund_df
  
}
