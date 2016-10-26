#' Calculation for descriptors of populations
#'
#' \code{pop_calcs} Calculations for describing features of populations.
#'
#' @description Given a population of integration sites, TCR/IGH sequences,
#' etc., a numerical value can be used to describe the diveristy (shannon index
#' or entropy), or the clonality (gini or clonality). These numerical values
#' allow for a comparison to be made between populations. Entropy and Clonality
#' are factors used by Adaptive Biotechnologies to describe their TCR / IGH
#' sequence data.
#'
#' @usage
#' pop_calcs(x, calc)
#'
#' @param x numeric a vector of abundances or proportions / frequencies.
#'
#' @param calc character one of the various population calcuations. Choices are:
#' "shannon", "gini", "entropy", or "clonality".
#'
#' @examples
#' x <- sample(1:100, 50, replace = TRUE)
#'
#' calc_shannon(x)
#' calc_gini(x)
#' calc_entropy(x)
#' calc_clonality(x)
#'
#' @author Christopher Nobles, Ph.D.
#' @export

pop_calcs <- function(x, calc){
  if(calc == "shannon"){
    calc_shannon(x)
  }else if(calc == "gini"){
    calc_gini(x)
  }else if(calc == "entropy"){
    calc_entropy(x)
  }else if(calc == "clonality"){
    calc_clonality(x)
  }
}

#' @describeIn pop_calcs Calculate Shannon Diversity.
calc_shannon <- function(x, base = exp(1)){
  x <- x[!is.na(x)]
  x <- x/sum(x)
  shannon <- -x*log(x, base)
  shannon <- sum(shannon)
  shannon
}
#' @describeIn pop_calcs Calculate gini index.
calc_gini <- function(x){
  stopifnot(require(reldist))
  x <- x[!is.na(x)]
  x <- x/sum(x)
  gini <- gini(x)
  gini
}

#' @describeIn pop_calcs Calculate Entropy as defined by Adaptive Biotech.
calc_entropy <- function(x){
  x <- x[!is.na(x)]
  x <- x/sum(x)
  entropy <- -sum(x*log(x, base = 2))
  entropy
}

#' @describeIn pop_calcs Calcuate Clonality as defined by Adaptive Biotech.
calc_clonality <- function(x){
  x <- x[!is.na(x)]
  x <- x/sum(x)
  clonality <- 1+sum(x*log(x, base = length(x)))
  clonality
}
