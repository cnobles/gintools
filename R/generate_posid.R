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
#' generate_posid(seqnames = NULL, strand = NULL, start = NULL, end = NULL, delim = "")
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
#' @param delim a single character to delimit the data within the string.
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
                           end=NULL, delim = ""){
  if(!is.null(sites)){
    if(class(sites) == "GRanges"){
      chr <- as.character(GenomicRanges::seqnames(sites))
      strand <- as.vector(GenomicRanges::strand(sites))
      pos <- ifelse(
        strand != "-", GenomicRanges::start(sites), GenomicRanges::end(sites))
    }else{
      stop("Sites provided not a GRanges object, use alternative inputs.")
    }
  }else if(is.null(sites)){
    if(!is.null(seqnames) & !is.null(strand) & !is.null(start) & !is.null(end)){
      if(length(unique(sapply(
        list(seqnames, strand, start, end), length))) != 1){
          stop("Input seqnames, strand, start, and end are not equal length.")
      }
      chr <- as.character(seqnames)
      strand <- as.vector(strand)
      start <- as.integer(start)
      end <- as.integer(end)
      pos <- ifelse(strand != "-", start, end)
    }else{
      chr <- strand <- pos <- delim <- c()
    }
  }

  if(length(chr) == 0 & length(strand) == 0 & length(pos) == 0){
    return(character())
  }else{
    return(paste0(chr, delim, strand, delim, pos))
  }
}

