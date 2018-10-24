#' A Binary Ambiguous Nucleotide scoring Matrix (BAN Mat)
#' 
#' \code{banmat} returns a binary ambiguous nucleodide matrix.
#' 
#' @description Constructed based on NUC4.4 and designed for comparing ambiguous
#' sequences against "A", "T", "G", "C", and "N" containing sequences. 
#' Currently matches between ambiuous nucleotides are considered mismatch.
#' 
#' @usage 
#' banmat()
#' 
#' @example 
#' banmat()
#'  
#' @author Christopher Nobles, Ph.D.
#' @export

banmat <- function(){
  matrix(c(
    1,0,0,0,0,1,1,0,0,1,0,1,1,1,1,0,
    0,1,0,0,0,1,0,1,1,0,1,0,1,1,1,0,
    0,0,1,0,1,0,1,0,1,0,1,1,0,1,1,0,
    0,0,0,1,1,0,0,1,0,1,1,1,1,0,1,0,
    0,0,1,1,1,0,0,0,0,0,0,0,0,0,1,0,
    1,1,0,0,0,1,0,0,0,0,0,0,0,0,1,0,
    1,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,
    0,1,0,1,0,0,0,1,0,0,0,0,0,0,1,0,
    0,1,1,0,0,0,0,0,1,0,0,0,0,0,1,0,
    1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,
    0,1,1,1,0,0,0,0,0,0,1,0,0,0,1,0,
    1,0,1,1,0,0,0,0,0,0,0,1,0,0,1,0,
    1,1,0,1,0,0,0,0,0,0,0,0,1,0,1,0,
    1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,0,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
    ncol = 16,
    nrow = 16,
    byrow = TRUE,
    dimnames = list(
      c("A", "T", "G", "C", "S", "W", "R", "Y", 
        "K", "M", "B", "V", "H", "D", "N", "?"),
      c("A", "T", "G", "C", "S", "W", "R", "Y", 
        "K", "M", "B", "V", "H", "D", "N", "?")))
}
