#' Creates Pangenomes from mmseqs object
#'
#' Creates clusters of genomes (pan-genome) using a list of pre-defined clusters.
#' The list of cluster can be from external source (e.g. mlst or BAPS) or from an
#' internal function (e.g \emph{clustering} function). Moreover, pan-genomes can be
#' filtered by genome frequency (i.e. the number of genomes that belong to
#' the pangenome). pangenomes builds a new \emph{accnet} object that relates the pangenome
#' with the genes/proteins. The user can filter the proteins that are included by
#' the presence frequency.
#'
#'
#' @param data A mmseq object
#' @param cluster A data.frame with two columns: strain label (ID) and cluster.
#' Colnames are not required. Cluster column can be any categorical variable.
#' @param min_freq Minimun gene/protein frequency to include in the pangenomes.
#' @param max_freq Maximun gene/protein frequency to include in the pangenomes.
#' @param min_pangenome_size Minimun number of genomes to create a new pangenome.
#'
#' @return An \emph{accnet} object.
#' @export
#'
#' @seealso \code{\link{clustering}}
#'
#' @examples
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @import data.table
#'

pangenomes_from_mmseqs <-function(data,cluster,min_freq = 0,max_freq = 1,min_pangenome_size = 1)
{
  if(!is(data,"mmseq"))
  {
    stop("mmseqs must be a 'mmseq' object")
  }

  colnames(cluster) = c("Source", "Cluster")

  accnet <- data$table %>% as_tibble() %>% select(Source = Genome_genome, Target = Prot_prot) %>%
    distinct() %>%
    inner_join(cluster) %>%
    mutate(dummy = "Pangenome") %>%
    unite(dummy,
          Cluster,
          col = "Pange",
          sep = "_",
          remove = FALSE)
  n_genomes <- accnet %>% distinct(Source) %>% nrow(.)

  correspondence_table <- accnet %>% select(Source,Pange) %>% distinct()

  accnet <- accnet %>% group_by(Pange) %>%
    mutate(PangeSize = n_distinct(Source)) %>%
    ungroup() %>%
    group_by(Pange, Target, PangeSize) %>%
    summarise(ProtFreq = n()) %>%
    distinct(.keep_all = FALSE) %>%
    ungroup() %>% group_by(Pange) %>% mutate(PangeDegree = n()) %>%
    ungroup() %>% group_by(Target) %>% mutate(ProtDegree = n()) %>%
    ungroup() %>%
    filter(ProtDegree >= min_freq) %>%
    filter(PangeSize >= min_pangenome_size) %>%
    filter(ProtDegree <= max_freq * n_genomes)

  accnet.matrix <- accnet %>%
    mutate(value = 1) %>%
    select(Pange, Target, ProtFreq) %>%
    spread(Target, ProtFreq, fill = 0) %>% rename(Source = Pange)


  Annotation <-
    data$annot %>% semi_join(accnet %>% select(Target) %>% distinct() %>% rename(Prot_prot = Target),
                             by = "Prot_prot")
  results <- list(
    list = accnet %>% rename(Source = Pange) %>% select(Source, Target, PangeSize, ProtFreq, PangeDegree, ProtDegree),
    matrix = accnet.matrix,
    annot = Annotation %>% rename(ID = Prot_prot,),
    membership = correspondence_table
  )
  class(results) <- "accnet"
  return(results)
}
