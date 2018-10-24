#' Collapse row contents of a data.frame or matrix into single vector.
#'
#' \code{vcollapse} returns a single vector from input data.frame or matrix 
#' where row contents have been combined or collapsed, as with 
#' `paste(..., collaspe = "")`.
#'
#' @description Similar to python zip, `vzip` takes input vectors and merges
#' them together by their input order and index. A simple example is two numeric
#' vectors, A = c(1,1,1) and B = c(2,2,2). The output of vzip(A,B) would simply
#' be a single vector of c(1,2,1,2,1,2). Any number of vectors can be input, but
#' each input vector must be of the same length. Output vector class depends on
#' input vector consensus.
#'
#' @usage
#' vcollapse(d)
#' vcollapse(d, sep = "-", fill = "NA")
#'
#' @param d data.frame or matrix or object coercible to a matrix. Row contents
#' will be combined into a single output vector.
#' @param sep character used to separate collapsed contents.
#' @param fill character used to fill empty values within the coerced object.
#'
#' @examples
#' df <- data.frame(
#'   "A" = letters[1:5],
#'   "B" = 3:7,
#'   "C" = LETTERS[2:6])
#' vcollapse(df)
#' vcollapse(df, sep = "-")
#' vcollapse(df, sep = "-", fill = "z")
#'
#' @author Christopher Nobles, Ph.D.
#' @export
#'

vcollapse <- function(d, sep = "", fill = NULL){
  if(is.vector(d)){
    stop("Function vcollapse() is not used on vectors, use paste(collapse = ...).")
  }
  if(any(sapply(seq_len(ncol(d)), function(i) class(d[,i])) == "factor")){
    fct_idx <- which(
      sapply(seq_len(ncol(d)), function(i) class(d[,i])) == "factor")
    mod_env <- new.env()
    mod_env$d <- d
    null <- lapply(fct_idx, function(i){
      mod_env$d[,i] <- as.character(d[,i])
    })
    d <- mod_env$d
  }
  if(class(d) != "matrix") d <- gsub(" ", "", as.matrix(d))
  mat <- d
  if(!is.null(fill)) mat <- ifelse(is.na(mat), fill, mat)
  mat <- do.call(
    cbind, 
    lapply(seq_len(ncol(mat)), function(i){
      if(i < ncol(mat)){
        cbind(mat[,i], rep(sep, nrow(mat)))
      }else{
        mat[,i]
      } 
    }))
  mat <- cbind(mat, rep(">|<", nrow(mat)))
  if(any(is.na(mat))){
    stop("NA values detected in object, please use fill param for vcollapse.")
  }
  div_str <- stringr::str_c(t(mat), collapse = "")
  unlist(strsplit(div_str, split = ">\\|<"))
}
