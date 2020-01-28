#' Similarity Network
#'
#' The function creates a similarity network connecting genomes
#' with 'distances' lower than the especified threshold.
#'
#' @param data An object of class \emph{accnet} or \emph{mash}
#' @param threshold Maximun distance to include in the network.
#'
#' @details The distance must to be addecuated for the type of data.
#' MASH values rangin between 0-0.05 intra-specie mean while accnet distances
#' could be higher.
#'
#' @return An \emph{igraph} object.
#' @export
#'
#' @import dplyr
#' @import igraph
#' @import parallelDist
#'
#' @examples
similarity_network <- function(data, threshold)
{
  if(is(data,"mash"))
  {

    gr <- data$table %>%
      filter(Dist <= threshold) %>%
      graph_from_data_frame(directed = FALSE) %>%
      igraph::simplify(remove.multiple = TRUE,remove.loops = TRUE)
    return(gr)

  }else if (class(data) == "accnet")
  {
    gr <- data$matrix %>% parallelDist(method ="binary") %>% as.matrix()
    gr = gr < threshold
    gr = graph_from_adjacency_matrix(gr, mode ="upper")
  }

}
