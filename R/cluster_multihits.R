#' Identified clusters for multihit integration sites.
#'
#' \code{cluster_multihits} returns a GRanges object with cluster information.
#'
#' @description Given a GRanges object of integration sites which align to
#' multiple location on the reference genome (multihit sites), this function
#' will append a column with cluster designation. Multihits are considered to be
#' part of the same cluster if they share any alignments with other multihits.
#' This function tackles this memory intensive process through iterive measures,
#' involving randomly sampling from the given set of data to start cluster
#' formation, refinement, and then resamples from the remaining data till all
#' information has been used.
#'
#' @usage
#' cluster_multihits(multihits, read_col = NULL)
#'
#' cluster_multihits(multihits, read_col = NULL, max_gap = 5L, iterations = 5L)
#'
#' @param multihits a GRanges object containing sets of integration sites
#' or ranges (one alignment per row) and containing the column from 'read_col'.
#'
#' @param read_col character string matching the name of the column of the
#' GRanges object given in 'multihits' which designates which ranges are
#' associated together. For example, the read name can be used here to show
#' which alignments were associated with the same read.
#'
#' @param gap an integer designating the nucleotide distance or window for which
#' to group integration sites.
#'
#' @param iterations integer The number of interations of cluster generation to
#' perform and the number of random subsets that will be made from the data.
#' More iterations leads to less overhead memory and more time required.
#'
#' @examples
#'
#' @author Christopher Nobles, Ph.D.
#' @export

cluster_multihits <- function(multihits, read_col = NULL,
                              max_gap = 5L, iterations = 5L){

}

