#' Find singles (genes/proteins) in a AccNET object
#'
#' \emph{singles()} return the genes/proteins that are only present in each
#' genome or cluster of genomes.
#'
#'
#' @param data An \emph{accnet} object.
#' @param cluster A data.frame with two columns. The first one must match with
#' \emph{Source} column of Accnet object (usually the genomes/samples). Second one
#' must contain a categorical variable (cluster, group, ST, origin, source ...)
#'
#' @return A data.frame with Genome/cluster, the genes/proteins belong to
#' Genome/Cluster and its annotation.
#' @export
#'
#' 
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @importFrom data.table fread
#'
singles <- function(data, cluster)
{
  if(!is(data,"accnet"))
  {
    stop("data must be an Accnet object")
  }

  if(missing(cluster))
  {
    table <-  data$list %>%
      group_by(Target) %>% mutate(Freq = n()) %>%
      filter(Freq ==1)
    return(table %>%
             select(Source,Target) %>%
             inner_join(data$annot %>%
                          rename(Target = ID)))
  }else{
    colnames(cluster) <- c("Source","Cluster")
    cluster <- cluster %>%
      distinct() %>%
      group_by(Cluster) %>%
      mutate(ClusterSize = n()) %>%
      ungroup()
    table <- data$list %>%
      inner_join(cluster, by ="Source") %>%
      group_by(Target) %>%
      mutate(TotalFreq =n()) %>% ungroup() %>%
      group_by(Cluster,Target) %>%
      mutate(ProtFreq = n()) %>%
      filter(ProtFreq == ClusterSize & ProtFreq == TotalFreq)
    return(table %>%
             select(Cluster,Target) %>%
             inner_join(data$annot %>%
                          rename(Target = ID)) %>%
             distinct())
  }


}
