#' Internal function to cluster igraph data
#'
#' @param data \emph{igraph} object
#' @param method greedy, louvain or walktrap
#'
#' @return \emph{data.frame} with two columns: \emph{Source} and \emph{Target}
#'
#' @import igraph
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom data.table fread
cluster_knnn <- function(data,method)
{

  if(method =="greedy")
  {
    cluster <- cluster_fast_greedy(data)
    return(data.frame(Source = cluster$names, Cluster = cluster$membership))
  }else if(method =="louvain")
  {
    cluster <- cluster_louvain(data)
    return(data.frame(Source = cluster$names, Cluster = cluster$membership))
  }else if(method =="walktrap")
  {
    cluster <- cluster_walktrap(data)
    return(data.frame(Source = cluster$names, Cluster = cluster$membership))
  }else{
    stop("Unrecognized method for igraph object")
  }
}
