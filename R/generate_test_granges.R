#' Generate a test GRanges object
#'
#' @usage .generate_test_granges(n_sites = 5, n_reads_p_site = 20, site_range = 1:1000, read_width_range = 30:100)
#'
#' @param n_sites number of sites to generate
#' @param n_reads_p_site number of reads per site
#' @param site_range numeric vector to sample from for postion of the integration
#' sites
#' @param read_width_range range of widths to sample from for each read
#' @param stdev numeric, standard deviation of the integration site distribution
#' around the selected positions
#' @param positions a numeric vector of preselected positions.
#'
#' @details Varies the integration start position in with a normal distribution
#' around the random integration site and then assignes a random uniform
#' distribution of widths for the sites.
#'
#' @examples
#' .generate_test_granges()
#'
#' @author Christopher Nobles, Ph.D.
#'

.generate_test_granges <- function(n_sites = 5, n_reads_p_site = 20,
                                   site_range = 1:1000,
                                   read_width_range = 30:100,
                                   stdev = 1,
                                   positions = NULL){
  if(is.null(positions)){
    positions <- sort(sample(site_range, n_sites, replace = TRUE))
  }

  message("True positions: ", paste(sort(positions), collapse = ", "))
  sort(GRanges(
    seqnames = rep("chr1", n_sites*n_reads_p_site),
    ranges = IRanges(
      start = sapply(positions, function(x){
        x + sample(round(rnorm(n_reads_p_site, mean = 0, sd = stdev)))
      }),
      width = sample(read_width_range, n_reads_p_site, replace = TRUE)),
    strand = rep("+", n_sites*n_reads_p_site)
  ))
}
