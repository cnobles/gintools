#' Identify sources in a directed graph
#'
#' \code{sources} returns a numerical vector of source nodes.
#'
#' @description From a directed graph, this function returns all source nodes, or
#' nodes which only act as head nodes and not as tail nodes.
#'
#' @usage
#' sources(graph)
#'
#' @param graph a directed graph (igraph).
#'
#' @examples
#' g <- make_graph(edges = c(1,2, 2,3, 1,3), directed = TRUE)
#' sources(g)
#'
#' @author Christopher Nobles, Ph.D.
#'

sources <- function(graph){
  if(!igraph::is_directed(graph)){
    message("Graph provided is not a directed graph.")
    srcs <- c()
  }else{
    srcs <- which(Matrix::colSums(
      igraph::get.adjacency(graph, sparse = TRUE)) == 0)
  }
  srcs
}
