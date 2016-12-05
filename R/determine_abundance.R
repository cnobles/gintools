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
#' gr <- .generate_test_granges()
#' std.gr <- standardize_intsites(gr)
#'
#' determine_abundance(std.gr)
#'
#' @author Christopher Nobles, Ph.D.
#' @export

determine_abundance <- function(sites, grouping = NULL, replicates = NULL,
                                method = "fragLen"){
  sites$posID <- generate_posid(sites)
  mcols <- mcols(sites)
  mcols(sites) <- NULL

  grp <- which(grouping == names(mcols))
  if(length(grp) != 0){
    group <- as.character(mcols[,grp])
  }else{
    group <- rep("group1", length(sites))
  }

  reps <- which(replicates == names(mcols))
  if(length(reps) != 0){
    replicates <- as.character(mcols[,reps])
  }else{
    replicates <- rep("1", length(sites))
  }

  sites.dfr <- data.frame("posID"= generate_posid(sites),
                          "fragLen"= width(sites),
                          "replicates"=replicates,
                          "group"=group)

  if(method == "fragLen"){
    abundCalc <- function(locations, fragLen, replicates, group){
      group <- unique(group)
      dfr <- data.frame("locationID" = locations,
                        "fragLen" = fragLen,
                        "replicates" = replicates)
      dfr_dist <- distinct(dfr)
      sites_list <- split(dfr_dist, dfr_dist$locationID)
      abundances <- sapply(sites_list, function(x){nrow(x)})
      abundances <- abundances[unique(dfr$locationID)]

      abund.dfr <- data.frame("posID" = names(abundances),
                              "estAbund" = abundances,
                              "group" = rep(group, length(abundances)))
      abund.dfr
    }
  }else if(method == "estAbund"){
    abundCalc <- function(locations, fragLen, replicates, group){
      group <- unique(group)
      if(length(unique(sites.dfr$replicates)) == 1){
        theta_list <- estAbund(locations=locations, lengths=fragLen)
      }else{
        theta_list <- estAbund(locations=locations, lengths=fragLen,
                               replicates=replicates)
      }
      posID <- names(theta_list$theta)
      abundances <- theta_list$theta
      if(length(unique(sites.dfr$replicates)) == 1){
        abund.dfr <- data.frame("posID" = posID,
                                "estAbund" = abundances)
      }else{
        abund.dfr <- data.frame("posID" = posID,
                                "estAbund" = abundances,
                                "replicates" = length(unique(
                                  theta_list$data$replicates)))
      }
      abund.dfr$group <- group
      abund.dfr
    }
  }else{
    stop("Must choose either fragLen or estAbund for method.")
  }

  sites.list <- split(sites.dfr, sites.dfr$group)

  abund.list <- lapply(sites.list, function(x){
    abundCalc(locations = x$posID,
              fragLen = x$fragLen,
              replicates = x$replicates,
              group = x$group)
  })

  abund.list <- lapply(abund.list, function(x){
    x <- x %>%
      dplyr::mutate(totalAbund = sum(x$estAbund)) %>%
      dplyr::mutate(estAbund = round(estAbund)) %>%
      dplyr::mutate(relAbund = estAbund/totalAbund) %>%
      dplyr::mutate(relRank = rank(-1*relAbund, ties.method="max"))
    x
  })

  abund.dfr <- bind_rows(abund.list)

  abund.dfr$totalAbund <- NULL

  abund.dfr <- abund.dfr[, c("posID", "group", "estAbund", "relAbund", "relRank")]
  abund.dfr
}
