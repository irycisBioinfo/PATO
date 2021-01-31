#' Find coincidents gene/proteins (clustering)
#'
#' This function find groups of coincidents genes/proteins in an \emph{accnet}
#' object. Unlike \emph{accnet_enrichment_analysis()} \emph{concidents} does not need
#' to define genome/pangenome clusters. The methos perform a multi-dimensional scaling
#' using \emph{umap} and then find clusters with \emph{hdbscan}.
#'
#' @param data An \emph{accnet} object
#' @param min_freq Minimun frequency of gene/protein in the population.
#' @param min_size Minimun size of the reported clusters.
#'
#' @return A data.frame with the protein/gene ID ("Target") and the membership (cluster number)
#' @export
#'
#' 
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @importFrom data.table fread
#' @import uwot
#' @import dbscan
#'
coincidents <- function(data, min_freq = 5, min_size = 4)
{

  if(!is(data,"accnet"))
  {
    stop("data must be an Accnet object")
  }


  tmp <- data$list %>% group_by(Target) %>%
    mutate(Freq = n()) %>% filter(Freq > min_freq) %>%
    mutate(Value =1) %>%
    pivot_wider(names_from = Source, values_from = Value, values_fill = list(Value =0))

  um <- umap(tmp %>% column_to_rownames("Target"),metric = "manhattan", n_components = 2)
  hdb <- hdbscan(um,minPts = min_size)

  return(data.frame(Target = tmp$Target,Cluster = hdb$cluster))

}
