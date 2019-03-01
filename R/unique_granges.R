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
#' unique_granges(sites, sum.cols = FALSE, rm.cols = NULL, rm.dup.cols = NULL)
#'
#' @param sites A GRanges object with or without metadata columns.
#' @param sum.cols a logical or character vector of column name(s) in the 
#' metadata of the input GRanges object to sum across unique observations. 
#' These column(s) will not be considered when identifying unique observations. 
#' Default is FALSE, if set to TRUE, then will attempt to sum the column 
#' 'counts'.
#' @param rm.cols a character vector of column name(s) to remove from the 
#' metadata of the input GRanges object. Removal of columns will occur prior to 
#' identifying unique observations, and therefore these columns will not only 
#' be dropped from the output, but will not be considered in identifying unique 
#' observations.
#' @param rm.dup.cols logical or a regular expression to identify duplicate 
#' columns and remove the duplicates. This will only remove columns that are 
#' true duplicates, meaning they have identical content. Removal of duplicated 
#' columns will occur prior to identifying unique observations, and therefore 
#' identified duplicates will not only be dropped from the output, but will not 
#' be considered in identifying unique observations.
#'
#' @examples
#' gr <- gintools:::generate_test_granges(
#'   n_sites = 1,
#'   n_reads_p_site = 12,
#'   site_range = 1:20,
#'   read_width_range = 20:30
#' )
#' gr <- refine_breakpoints(gr)
#' gr <- standardize_sites(gr)
#' gr$sample <- rep(c("A","B"), 6)
#' gr$sample.1 <- rep(c("A","B"), 6)
#' gr$counts <- rep(1:4, c(3,3,3,3))
#' 
#' # Calling unique() on gr returns a miss interpreted data set
#' unique(gr)
#'
#' # Using unique_granges() without options returns all distinct rows
#' unique_granges(gr)
#'
#' # Using the options for sum.cols, rows 'counts' are added when combined 
#' # together.
#' unique_granges(gr, sum.cols = TRUE)
#' 
#' # Or multiple columns can be added simultaneously.
#' gr$tags <- rep(1:2, c(6,6))
#' unique_granges(gr, sum.cols = c("counts", "tags"))
#' 
#' # Remove specific columns by passing a character vector to 'rm.cols'
#' unique_granges(gr, sum.cols = c("counts", "tags"), rm.cols = "sample.1")
#' 
#' # Remove any column that may be a duplicate of another with a specific 
#' # name pattern. Column content will be checked for identity before removal.
#' unique_granges(gr, sum.cols = c("counts", "tags"), rm.dup.cols = "[\\w]+")
#'
#' @author Christopher L. Nobles, Ph.D.
#' @importFrom magrittr %>%
#' @export

unique_granges <- function(sites, sum.cols = FALSE, 
                           rm.cols = NULL, rm.dup.cols = NULL){
  
  # Checks and balance
  if( !class(sites) == "GRanges" ){
    stop("\n  Sites object is not a GRanges class.\n")
  }
  
  if( any(sum.cols != FALSE) ){
    
    if( all(sum.cols == TRUE) ){
      counts_col <- "counts"
    }else{
      counts_col <- sum.cols
    }

    if( any(!counts_col %in% names(GenomicRanges::mcols(sites))) ){
      missing_cols <- counts_col[
        which(!counts_col %in% names(GenomicRanges::mcols(sites)))
      ]
      stop(
        "\n  Could not find column name(s) in sites object: ",
        paste(missing_cols, collapse = ", "), 
        "\n"
      )
    }
    
  }
  
  if( !is.null(rm.cols) ){
    
    if( any(!rm.cols %in% names(GenomicRanges::mcols(sites))) ){
      
      missing_cols <- rm.cols[
        which(!rm.cols %in% names(GenomicRanges::mcols(sites)))
        ]
      
      message(
        "Could not find column name(s) to remove in sites object: ",
        paste(missing_cols, collapse = ", ")
      )
      
      rm.cols <- rm.cols[rm.cols %in% missing_cols]
      
    }
    
    if( length(rm.cols) > 0 ){
      GenomicRanges::mcols(sites)[, rm.cols] <- NULL
    }
    
  }
  
  if( !is.null(rm.dup.cols) ){
    
    col_names <- names(GenomicRanges::mcols(sites))
    
    if( is.character(rm.dup.cols) ){
      col_names <- stringr::str_extract(col_names, rm.dup.cols)
    }
    
    dup_cols <- names(table(col_names))[table(col_names) > 1]
    dup_idxs <- lapply(dup_cols, function(x) which(col_names == x))
    
    rm_col_idx <- unlist(lapply(dup_idxs, function(x){
      
      col_digests <- unlist(lapply(x, function(y){
        digest::digest(GenomicRanges::mcols(sites)[,y, drop = TRUE])
      }))
      
      x[which(duplicated(col_digests))]
      
    }))
    
    GenomicRanges::mcols(sites)[,rm_col_idx] <- NULL
    
    if( is.character(rm.dup.cols) ){
      
      names(GenomicRanges::mcols(sites)) <- stringr::str_extract(
        names(GenomicRanges::mcols(sites)), rm.dup.cols
      )
      
    }

  }

  # Convert sites to a data.frame and remove duplicates
  if( !length(names(sites)) == length(unique(names(sites))) ){
    
    message("Dropping rownames for data.frame conversion.")
    df <- GenomicRanges::as.data.frame(sites, row.names = NULL)
    
  }else{
    
    df <- GenomicRanges::as.data.frame(sites)
    
  }
  
  cols <- names(df)

  if( any(sum.cols != FALSE) ){
    counts_pos <- match(counts_col, cols)
  }

  # Sum counts if needed
  if( !any(sum.cols != FALSE) ){
    
    df <- dplyr::distinct(df)
    
  }else{
    
    groups <- lapply(cols[-counts_pos], as.symbol)
    
    df <- dplyr::group_by(df, .dots = groups) %>%
      dplyr::summarise_all(.funs = "sum") %>%
      dplyr::ungroup()
    
    df <- df[,cols]
    
  }

  # Rebuild GRanges object
  gr <- GenomicRanges::GRanges(
    seqnames = df$seqnames,
    ranges = IRanges::IRanges(start = df$start, end = df$end),
    strand = df$strand,
    seqinfo = GenomicRanges::seqinfo(sites)
  )

  GenomicRanges::mcols(gr) <- dplyr::select(df, 6:length(df))
  gr
}
