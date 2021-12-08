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
#' @param threshold Minimum value to create an edge.
#'
#' @return Returns an \emph{igraph} object.
#' @export
#'
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @importFrom data.table fread
#' @import igraph
#' @import parallelDist

knnn <- function(data, n_neigh, repeats = TRUE, threshold = 1)
{
  if(is(data,"accnet"))
  {

    if(is.null(data$dist))
    {
      matrix <- data$matrix %>%
        column_to_rownames("Source")%>% as.matrix() %>%
        parallelDist(., method = "binary") %>%
        as.matrix() %>%
        as.data.frame()
    }else{
      matrix = data$dist %>% as.data.frame()
    }
    List <- matrix %>%
      rownames_to_column("Source") %>%
      gather(Target,Dist, -Source)

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
                                          filter(Dist < threshold) %>%
                                  group_by(Source) %>%
                                  slice_min(Dist, n= n_neigh+1,with_ties =F) %>%
                                  mutate(weight = 2-Dist) %>%
                                  select(-Dist))

    gr <- igraph::as.undirected(gr)
    gr <- igraph::simplify(gr,remove.multiple = TRUE,remove.loops = F)

  }else
  {
    matrix <- matrix %>%
      as.matrix()
    matrix[lower.tri(matrix, diag = F)] <- NA
    matrix <- as.data.frame(matrix) %>% rownames_to_column("Source")

    List.U<- matrix %>%
      pivot_longer(names_to = "Target",values_to = "Dist", -Source) %>%
      filter(!is.na(Dist)) ## remove the lower.tri

    gr <- igraph::graph_from_data_frame(List.U %>%
                                  filter(Dist < threshold) %>%
                                  group_by(Source) %>%
                                  slice_min(Dist, n= n_neigh+1,with_ties =F) %>%
                                  mutate(weight = 2-Dist) %>%
                                  select(-Dist), directed = FALSE)
    gr <- igraph::simplify(gr,remove.multiple = TRUE,remove.loops = F)
  }
  return(gr)
}
