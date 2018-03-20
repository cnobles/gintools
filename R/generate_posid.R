#' Generate position IDs.
#'
#' \code{generate_posid} generates a character vector of position IDs
#'
#' @description Given a set of genomic integration positions, the
#' function generates a character vector of integrated positions in the
#' following format: chr(+/-)pos, where chr is the chromosome, (+/-) designates
#' the orientation of the integration, and pos is the position of the integrated
#' element.
#'
#' @usage
#' generate_posid(sites = NULL)
#'
#' generate_posid(seqnames = NULL, strand = NULL, start = NULL, end = NULL)
#'
#' @param sites a GRanges object where each row represents one integrated
#' element.
#' @param seqnames a character vector of chromosomes for integrated elements,
#' i.e. chr1, chr3, chrX.
#' @param strand orientation or strand on which the integrated element
#' is on (+, -, *)
#' @param start the lower numerical position of the integrated element,
#' considered the start for "+" integrated elements.
#' @param end the greater numerical position of the integrated element,
#' considered the end for "+" integrated elements.
#'
#' @examples
#' chr <- c("chr1", "chr3", "chrX")
#' strands <- c("+", "-", "+")
#' starts <- c(900231, 13254892, 603292)
#' ends <- c(900431, 13255292, 603592)
#' ranges <- IRanges(start = starts, end = ends)
#' gr <- GRanges(seqnames = chr, ranges = ranges, strand = strands)
#'
#' generate_posid(sites = gr)
#'
#' generate_posid(seqnames = chr, strand = strands, start = starts, end = ends)
#'
#' @author Christopher Nobles, Ph.D.
#'
#' @export

generate_posid <- function(sites=NULL, seqnames=NULL, strand=NULL, start=NULL,
                           end=NULL){
  if(!is.null(sites) & length(sites) != 0){
    if(class(sites) == "GRanges"){
      chr <- as.character(GenomicRanges::seqnames(sites))
      strand <- as.vector(GenomicRanges::strand(sites))
      pos <- ifelse(
        strand == "+", GenomicRanges::start(sites), GenomicRanges::end(sites))
      posID <- paste0(chr, strand, pos)
    }else{
      message("Sites provided not a GRanges object,
              please use alternative inputs.")
      stop()
    }
  }else if(!is.null(sites) & length(sites) == 0){
    posID <- character()
  }else if(is.null(sites)){
    if(!is.null(seqnames) & !is.null(strand) & !is.null(start) & !is.null(end)){
      chr <- as.character(seqnames)
      strand <- as.vector(strand)
      start <- as.integer(start)
      end <- as.integer(end)
      sites.df <- data.frame(chr, strand, start, end)
      sites.df$pos <- ifelse(strand == "+", sites.df$start, sites.df$end)
      posID <- paste0(sites.df$chr, sites.df$strand, sites.df$pos)
    }else{
      message("Please supply seqnames, strand, start, and end info.")
      stop()
    }}
  return(posID)
}
