#' Determine all unique rows in GRanges object
#'
#' \code{unique_granges} returns a GRanges object of only unique observations (
#' removing all duplicated rows), yet considers the meta data information.
#'
#' @description Given a GRanges object, this function returns all unique rows,
#' or observations, including the meta data information. Calling the function
#' \code{unique()} on a GRanges object only returns unique ranges, and does not
#' account for meta data information, unlike to a data.frame. Using the option
#' sum.counts = TRUE and specifying the counts.col = "" name, will sum the
#' numerical values within the column for all combined rows.
#'
#' @usage
#' unique_granges(sites)
#'
#' unique_granges(sites, sum.counts = TRUE, counts.col = "counts")
#'
#' @param sites A GRanges object with or without metadata columns.
#' @param sum.counts logical, if TRUE, the column specified by counts.col will
#' be summed for all rows combined. If FALSE, all columns are considered as
#' unique values.
#' @param counts.col character string specifying the name of the column to use
#' as counts.
#'
#' @examples
#' gr <- .generate_test_granges(
#'   n_sites = 1,
#'   n_reads_p_site = 12,
#'   site_range = 1:20,
#'   read_width_range = 20:30
#' )
#' gr <- refine_breakpoints(gr)
#' gr <- standardize_intsites(gr)
#' gr <- granges(gr)
#' gr$sample <- rep(c("A","B"), 6)
#' gr$counts <- rep(1:4, c(3,3,3,3))
#'
#' # Calling unique() on gr returns a miss interpreted data set
#' unique(gr)
#'
#' # Using unique_granges() without options returns all distinct rows
#' unique_granges(gr)
#'
#' # Using the options for sum.counts and counts.col, rows counts are added when
#' # combined together.
#' unique_granges(gr, sum.counts = TRUE, counts.col = "counts")
#'
#' @author Christopher L. Nobles, Ph.D.
#' @export

unique_granges <- function(sites, sum.counts = FALSE, counts.col = NULL){
  require(dplyr)
  # Checks and balance
  if(!class(sites) == "GRanges"){
    stop("Sites object is not a GRanges class.")}
  if(sum.counts & is.null(counts.col)){
    stop("Please specify the names of the column with count information.")}
  if(!is.null(counts.col)){
    if(!counts.col %in% names(mcols(sites))){
      stop("Could not find counts column name in sites object.")}}

  # Convert sites to a data.frame and remove duplicates
  df <- GenomicRanges::as.data.frame(sites)
  cols <- names(df)

  if(sum.counts){
    counts.pos <- match(counts.col, cols)}

  # Sum counts if needed
  if(!sum.counts){
    df <- dplyr::distinct(df)
  }else{
    df$counts <- df[,cols[counts.pos]]
    groups <- lapply(cols[-counts.pos], as.symbol)
    df <- dplyr::group_by_(df, .dots = groups) %>%
      dplyr::summarise(counts = sum(counts))
    names(df) <- c(cols[-counts.pos], cols[counts.pos])
  }

  # Rebuild GRanges object
  gr <- GRanges(
    seqnames = df$seqnames,
    ranges = IRanges(start = df$start, end = df$end),
    strand = df$strand,
    seqinfo = seqinfo(sites)
  )

  mcols(gr) <- df[,6:length(df)]
  gr
}