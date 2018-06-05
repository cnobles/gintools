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
#' "shannon", "gini", "entropy", "clonality", or "uc50".
#'
#' @examples
#' x <- sample(1:100, 50, replace = TRUE)
#'
#' pop_calcs(x, calc = "shannon")
#' pop_calcs(x, calc = "gini")
#' pop_calcs(x, calc = "entropy")
#' pop_calcs(x, calc = "clonality")
#' pop_calcs(x, calc = "uc50")
#'
#' @author Christopher Nobles, Ph.D.
#'
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
  }else if(calc == "uc50"){
    calc_uc50(x)
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
calc_gini <- function(x, wt = NULL){
  if(is.null(wt)) wt <- rep(1, length(x))
  
  not_na <- which(!is.na(x))
  wt <- wt[not_na]
  x <- x[not_na]

  wt <- wt[order(x)]/sum(wt)  
  x <- x[order(x)]/sum(x)

  cum_wt <- cumsum(wt)
  cum_prod <- cumsum(x * wt)
  
  rel_prod <- cum_prod / cum_prod[length(cum_prod)]
  sum(rel_prod[-1] * cum_wt[-length(cum_wt)]) -
    sum(rel_prod[-length(rel_prod)] * cum_wt[-1])
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

#' @describeIn pop_calcs Calculate UC50 or the number of must abundant clones
#' which make up half the population.
calc_uc50 <- function(x){
  stopifnot(is.vector(x) & is.numeric(x))
  x <- x[order(x)]
  accum <- sapply(1:length(x), function(i){sum(x[1:i])})
  length(accum[accum >= sum(x)/2])
}
