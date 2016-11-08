#' Normalize integration site positions between samples
#'
#' \code{normalize_intsite_positions} identifies sites within a gap distance and
#' assigns each of the sites in the group the most frequent position as well as
#' a positionID (pos.clus).
#'
#' @description Integration sites across multiple samples may be identified at
#' slightly different nucleotide positions. This function identifies sites
#' within the gap distance from eachother and normalized their position to the
#' most frequent position. Similar to standardizing the integration sites, but
#' not as logic heavy.
#'
#' @usage
#' normalize_intsite_positions(sites)
#'
#' normalize_intsite_positions(sites, gap = 0L)
#'
#' @param sites GRanges object where each row represents one integration range
#' or site.
#'
#' @param gap integer for the distance between sites to consider.
#'
#' @examples
#'
#' gr <- GRanges(
#'   seqnames = rep("chr1", 10),
#'   ranges = IRanges(
#'     start = rnorm(10, 20, 1),
#'     width = rep(15, 10)),
#'   strand = rep("+", 10))
#'
#' normalize_intsite_positions(gr)
#'
#' normalize_intsite_positions(gr, gap = 1L)
#'
#' @author Christopher Nobles, Ph.D.
#'

normalize_intsite_positions <- function(sites, gap = 0L){
  sites$called.pos <- ifelse(strand(sites) == "+", start(sites), end(sites))
  sites$called.bp <- ifelse(strand(sites) == "+", end(sites), start(sites))
  sites <- flank(sites, -1, start = TRUE)

  graph.gap <- graph_clusters(sites, gap = gap)
  clusters.gap <- clusters(graph.gap)
  clusters_to_normalize <- grep(TRUE, clusters.gap$csize > 1)

  sites$pos.clus <- clusters.gap$membership

  norm.clus.dfl <- lapply(clusters_to_normalize, function(pos.clus){
    obs <- grep(TRUE, sites$pos.clus == pos.clus)
    dfr <- data.frame("strand" = strand(sites[obs]),
                      "called.pos" = sites[obs]$called.pos,
                      "pos.clus" = sites[obs]$pos.clus,
                      stringsAsFactors = FALSE)
    dfr
  })

  #Generate a key of the clusters and normalized positions
  norm.key <- as.data.frame(bind_rows(lapply(norm.clus.dfl, function(pos.dfr){
    pos.freq <- as.data.frame(table(pos.dfr$called.pos))
    top.freq <- as.numeric(as.character(
      pos.freq[pos.freq$Freq == max(pos.freq$Freq), "Var1"]
    ))
    norm.pos <- ifelse(unique(pos.dfr$strand) == "+", min(top.freq), max(top.freq))
    key.row <- data.frame("norm.pos" = norm.pos,
                          "pos.clus" = unique(pos.dfr$pos.clus),
                          stringsAsFactors = FALSE)
    key.row
  })))
  unnorm.sites <- grep(TRUE, !sites$pos.clus %in% clusters_to_normalize)
  unnorm.key <- data.frame("norm.pos" = sites[unnorm.sites]$called.pos,
                           "pos.clus" = sites[unnorm.sites]$pos.clus,
                           stringsAsFactors = FALSE)
  clus.pos.key <- as.data.frame(bind_rows(norm.key, unnorm.key))
  row.names(clus.pos.key) <- clus.pos.key$pos.clus
  sites$norm.pos <- clus.pos.key[as.character(sites$pos.clus), "norm.pos"]

  ranges <- IRanges(start = ifelse(strand(sites) == "+",
                                   sites$norm.pos,
                                   sites$called.bp),
                    end = ifelse(strand(sites) == "+",
                                 sites$called.bp,
                                 sites$norm.pos))
  normalized.sites <- GRanges(seqnames = seqnames(sites),
                              ranges = ranges,
                              strand = strand(sites),
                              seqinfo = seqinfo(sites))
  mcols(normalized.sites) <- mcols(sites)
  normalized.sites
}
