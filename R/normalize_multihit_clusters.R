#' Normalize multihit cluster IDs from multiple samples.
#'
#' \code{normalize_multihit_clusters} will normalize multihit clusterIDs across
#' multiple samples so that multihit sites can be identified across time points,
#' cell types, etc.
#'
#' @description As the INSPIIRED pipeline calls multihits, or integration sites
#' that can not be placed in a single location on the reference genome, it
#' assigns multihitID's to the various locations that the integration site may
#' exist. As each replicate is individually analyzed, multihitIDs for each
#' replicate are different, even through they may refer to the same integration
#' site. For this reason, normalize_multihit_clusters uses the previously
#' assigned multihitID and genomic positions to reassign multihitIDs across
#' multiple samples. Input for the function needs to be a GRanges object with a
#' metadata column labeled as "multihitid".
#'
#' @usage
#' normalize_multihit_clusters(multihits.gr)
#'
#' normalize_multihit_clusters(multihits.gr, gap = 5L, grouping = NULL)
#'
#' @param multihits.gr GRanges object with a column named 'multihitid'.
#'
#' @param gap integer designating the range to which consider sites identical.
#'
#' @param grouping Character, name of the column used to assign groups that will
#' not be compared to one another. Such as 'patient'.
#'
#' @examples
#' dfr <- data.frame(
#'   "chr" = c("chr1", "chr2", "chr2", "chr3"),
#'   "position" = c(5379927, 92775920, 2719573, 7195924),
#'   "breakpoint" = c(5380070, 92775995, 2719450, 7195890),
#'   "strand" = c("+", "+", "-", "-"),
#'   "sampleName" = rep("GTSP1234-1", 4),
#'   stringsAsFactors = FALSE)
#'
#' gr1 <- granges(db_to_granges(dfr))
#' gr2 <- gr1
#' gr1$multihitid <- c(1, 1, 2, 3)
#' gr2$multihitid <- c(4, 4, 5, 6)
#' gr1$patient <- rep(1, 4)
#' gr2$patient <- rep(2, 4)
#' gr <- normalize_intsite_positions(c(gr1, gr2))
#'
#' normalize_multihit_clusters(gr)
#'
#' # Group by patient will keep the two samples from being normalized to
#' # eachother
#'
#' normalize_multihit_clusters(gr, grouping = 'patient')
#'
#' @author Christopher Nobles, Ph.D.
#' @export

normalize_multihit_clusters <- function(multihits.gr, gap = 5L, grouping = NULL){
  #Multihits must be standardized and have clusterID info
  if(is.null(grouping)){
    multihits.gp <- GRangesList(multihits.gr)
  }else if(grouping %in% names(mcols(multihits.gr))){
    groups <- mcols(multihits.gr)[
      grep(grouping, names(mcols(multihits.gr)))]
    multihits.gr$groups <- groups[,1]
    multihits.gp <- split(multihits.gr, multihits.gr$groups)
  }else{
    stop("Grouping partitioning failed. Make sure grouping is either NULL or
         refering to the correct column in GRanges object.")
  }

  base_clus <- 0 # Start count for multihitid's

  norm.multi.grl <- lapply(multihits.gp, function(gr){
    key <- unique(data.frame(
      "multihitid" = as.character(gr$multihitid),
      "clusid" = gr$pos.clus,
      stringsAsFactors = FALSE
    ))

    split.key <- split(key, key$clusid)
    edgelist <- matrix(c(
      as.character(Rle(
        values = sapply(split.key, function(x) x$multihitid[1]),
        lengths = sapply(split.key, nrow)
      )),
      unlist(sapply(split.key, function(x) x$multihitid))),
      ncol = 2
    )
    graph <- graph.edgelist(edgelist, directed = FALSE)
    norm.multi.id <- as.numeric(clusters(graph)$membership) + base_clus
    names(norm.multi.id) <- names(clusters(graph)$membership)
    key$norm.multihitid <- norm.multi.id[key$multihitid]
    key$clusid <- NULL
    key <- unique(key)
    base_clus <<- max(norm.multi.id)

    norm.key <- data.frame(
      row.names = key$multihitid,
      "norm.multihitid" = key$norm.multihitid
    )

    gr$norm.multihitid <- norm.key[as.character(gr$multihitid), "norm.multihitid"]
    gr
  })

  norm.multi.gr <- unlist(GRangesList(norm.multi.grl))
  if(!is.null(grouping)) mcols(norm.multi.gr)$groups <- NULL
  norm.multi.gr
}

