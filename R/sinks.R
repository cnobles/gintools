#' Identify sinks in a directed graph
#'
#' \code{sinks} returns a numerical vector of sink nodes.
#'
#' @description From a directed graph, this function returns all sink nodes, or
#' nodes which only act as tail nodes and not as head nodes.
#'
#' @usage
#' sinks(graph)
#'
#' @param graph a directed graph (igraph).
#'
#' @examples
#' g <- make_graph(edges = c(1,2, 2,3, 1,3), directed = TRUE)
#' sinks(g)
#'
#' @author Christopher Nobles, Ph.D.
#'

sinks <- function(graph){
  if(!igraph::is_directed(graph)){
    message("Graph provided is not a directed graph.")
    snks <- c()
  }else{
    snks <- which(Matrix::rowSums(
      igraph::get.adjacency(graph, sparse = TRUE)) == 0)
  }
  snks
}
