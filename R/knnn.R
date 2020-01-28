#' K-Nearest Neighbour Network
#'
#' This function creates a network with the \emph{k} best neigbours
#' of each genome (or pangenome). The network can be buid with the
#' best neigbours with or without repetitions. The option \emph{repeats}
#' establish if the best \emph{k} neigbours includes bi-directional links
#' or not. Let be G_i and G_j two genomes, if G_i is one of the best
#' k-neighbours of G_j and G_j is one of the best k-neighbours of G_i if the
#' \emph{repeats} option is TRUE then G_i is removed from the k-neighbour and
#' substituted be the next best neighbour in other case G_i keeps in the list.
#'
#'
#' @param data An Accnet or Mash object
#' @param n_neigh The number of best K-Neighbours
#' @param repeats \emph{Boolean} Include repetitions?
#'
#' @return Returns an \emph{igraph} object.
#' @export
#'
#' @examples
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @import data.table
#' @import igraph
#' @import parallelDist

knnn <- function(data, n_neigh, repeats = TRUE)
{
  if(is(data,"accnet"))
  {

    matrix <- data$matrix %>%
      column_to_rownames("Source")%>%
      parallelDist(., method = "binary") %>%
      as.matrix() %>%
      as.data.frame()
    List <- matrix %>%
      rownames_to_column("Source") %>%
      gather(Target,Dist, -Source) %>%
      filter(Dist >0)

  }else if(is(data,"mash"))
  {
    List <- data$table
    matrix <- data$matrix

  }else{
    stop("data must be a mash or accnet object")

  }

  if(repeats)
  {
    gr <- igraph::graph_from_data_frame(List %>%
                                  group_by(Source) %>%
                                  top_n(-n_neigh,Dist) %>%
                                  mutate(weight = 2-Dist) %>%
                                  select(-Dist))
    gr <- igraph::as.undirected(gr)
    gr <- igraph::simplify(gr,remove.multiple = TRUE,remove.loops = TRUE)

  }else
  {
    matrix <- matrix %>%
      as.matrix()
    matrix[lower.tri(matrix)] <- 0
    matrix <- as.data.frame(matrix) %>% rownames_to_column("Source")

    List.U<- matrix %>%
      gather(Target,Dist, -Source) %>%
      filter(Dist >0)

    gr <- igraph::graph_from_data_frame(List.U %>%
                                  group_by(Source) %>%
                                  top_n(-n_neigh,Dist)%>%
                                  mutate(weight = 2-Dist) %>%
                                  select(-Dist), directed = FALSE)
    gr <- igraph::simplify(gr,remove.multiple = TRUE,remove.loops = TRUE)
  }
  return(gr)
}
