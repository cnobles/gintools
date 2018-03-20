#' Convert an INSPIIRED database query into a GRanges object.
#'
#' \code{db_to_granges} specifically converts a query from an INSPIIRED database
#' to a GRanges object.
#'
#' @description After querying data from an INSPIIRED integration site database,
#' this function can convert the returned data.frame into a GRanges object given
#' only the data.frame. The function generates a GRanges object containing the
#' seqnames, starts, ends, and strands information as well as the sampleName and
#' the specimen (parsed from the sampleName by the "-" delimiter). Additionally,
#' if there are other data in columns of the data.frame, they can be transfered
#' to the GRanges object by using the "keep.additional.columns = TRUE" option,
#' but will be dropped if FALSE. Additional columns will be have there column
#' names switched to lower case, but this feature can be turned off by using the
#' option 'lower_case = FALSE'.
#'
#' @usage
#' db_to_granges(dfr_from_db)
#'
#' db_to_granges(dfr_from_db, keep.additional.columns = TRUE, lower_case = TRUE)
#'
#' @param dfr_from_db a data.frame containing queried information from an
#' INSPIIRED integration site database. Any data.frame can be used, but requires
#' the columns: chr, position, breakpoint, strand, and sampleName.
#'
#' @param keep_additional_columns logical to specify whether to keep any columns
#' besides those required in the data.frame.
#'
#' @param lower_case logical, default is TRUE. Changes additional column names
#' to all lower case characters. FALSE leaves the additional column names as
#' original format.
#'
#' @param split_sampleName logical, default is TRUE. Changes behavior of how to
#' treat the sampleName column. If replicate information is not denoted by a
#' "-#" at the end of the sampleName, then creating the specimen list will fail.
#' Set to FALSE if replicate information is not included.
#'
#' @examples
#' dfr <- data.frame(
#'   "chr" = c("chr1", "chr2", "chr2", "chr3"),
#'   "position" = c(5379927, 92775920, 2719573, 7195924),
#'   "breakpoint" = c(5380070, 92775995, 2719450, 7195890),
#'   "strand" = c("+", "+", "-", "-"),
#'   "sampleName" = c("GTSP0001-3", "GTSP0003-1", "GTSP0002-4", "GTSP0003-1"),
#'   "is.multihit" = c(FALSE, TRUE, FALSE, TRUE),
#'   stringsAsFactors = FALSE)
#'
#' db_to_granges(dfr_from_db = dfr)
#'
#' db_to_granges(dfr_from_db = dfr, keep_additional_columns = FALSE)
#'
#' @author Christopher Nobles, Ph.D.
#' @export

db_to_granges <- function(dfr_from_db, keep_additional_columns = TRUE,
                          lower_case = TRUE, split_sampleName = TRUE){
  dfr <- dfr_from_db
  dfr <- dfr[, !duplicated(names(dfr))]

  if(any(names(dfr) == "length")){
    breakpoints <- ifelse(dfr$strand == "+", dfr$length, -dfr$length)
    dfr$breakpoint <- dfr$position + breakpoints
  }

  ranges <- IRanges::IRanges(
    start = ifelse(dfr$strand == "+", dfr$position, dfr$breakpoint),
    end = ifelse(dfr$strand == "+", dfr$breakpoint, dfr$position))
  gr <- GenomicRanges::GRanges(seqnames = dfr$chr,
                ranges = ranges,
                strand = dfr$strand)

  if(split_sampleName){
    specimen.list <- strsplit(dfr$sampleName, split="-")
    mcols <- data.frame("sampleName" = dfr$sampleName,
                        "specimen" = sapply(specimen.list, "[[", 1),
                        stringsAsFactors = FALSE)
  }else{
    mcols <- data.frame("sampleName" = dfr$sampleName, stringsAsFactors = FALSE)
  }

  if(keep_additional_columns){
    std.columns <- c("sampleName", "position", "chr", "strand", "breakpoint")
    are.there <- match(std.columns, colnames(dfr))
    add.cols <- grep(TRUE, is.na(match(names(dfr), names(dfr[,are.there]))))
    cols <- data.frame(dfr[, add.cols])
    colnames(cols) <- colnames(dfr[add.cols])
    mcols <- cbind(mcols, cols)
  }

  GenomicRanges::mcols(gr) <- mcols
  if(lower_case) names(GenomicRanges::mcols(gr)) <- tolower(
    names(GenomicRanges::mcols(gr)))
  gr
}
