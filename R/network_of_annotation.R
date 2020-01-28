#' Network of Annotation
#'
#' This function builds a \emph{igraph} object with the results from \emph{annotate} function.
#' Allows to filter the results by identity and/or Evalue.
#'
#' @param data \emph{annotate} result.
#' @param min_identity Minimun identity.
#' @param max_evalue Max Evalue.
#'
#' @return \emph{igraph} object
#' @export
#'
#' @examples
#' @seealso \code{\link{annotate}}
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @import igraph
#' @import data.table
#'
network_of_annotation <- function(data, min_identity = 0.95, max_evalue = 1e-25)
{

  data %>%
    group_by(Genome,Protein) %>%
    summarise_all(first) %>%
    filter(pident >= min_identity) %>%
    filter(evalue <= max_evalue) %>%
    group_by(Genome,Gene) %>%
    summarise(pident = max(pident)) %>%
    rename(from = Genome, to = Gene, weight = pident) %>%
    graph_from_data_frame() %>% return()


}
