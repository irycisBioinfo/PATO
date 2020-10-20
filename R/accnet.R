#' Build accessory genome network (AccNET)
#'
#' Create an \emph{accnet} object.
#'
#' This function extract the accessory genome from a \emph{mmseq} object.
#' The \emph{mmseq} contains all the relationship among genomes (or other
#' genomic unit) and genes/proteins. \emph{accnet} extract the accsessory
#' genome removing the sof-coregenome. Hard-coregenome is the set of
#' proteins present in the 100\\% of the genomes. Soft-coregenome is the set
#' of proteins present in the \emph{threshold} of the genomes.
#'
#'
#' @param mmseqs An object od class \emph{mmseq}
#' @param threshold Remove proteins with a frequency higher than threshold (soft-coregenome).
#' @param singles keep proteins presence in just one sample?
#'
#' @return An accnet object.\cr
#'   \emph{accnet} object contains\cr
#'   \itemize{
#'     \item{\emph{list}}: Data.frame with the relationship of genome-protein
#'     \item{\emph{matrix}}: Adjacency matrix (presence/absence)
#'     \item{\emph{annot}}: Annotation table with the protein functional annotation
#'     }
#'
#' @seealso \emph{\link{mmseqs}}
#'
#' @references AcCNET (Accessory Genome Constellation Network): comparative genomics software for accessory genome analysis using bipartite networks
#' VF Lanza, F Baquero, F de la Cruz, TM Coque - Bioinformatics, 2017
#' @export
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @import data.table
#'
accnet <- function(mmseqs,threshold = 0.8, singles = TRUE)
{


  if(!is(mmseqs,"mmseq"))
  {
    stop("mmseqs must be a 'mmseq' object")
  }

  accnet <- mmseqs$table %>% as_tibble()
  accnet <- accnet %>% select(Prot_prot,Genome_genome) %>% distinct()
  n_genomes <- accnet %>% distinct(Genome_genome) %>% count() %>% as_tibble()

  n_genomes <- n_genomes$n

  accnet <- accnet %>% group_by(Prot_prot) %>% mutate(freq = n()) %>% filter(freq < threshold*n_genomes)

  if(!singles)
  {
    accnet <- accnet %>% filter(freq > 1)
  }

  accnet.matrix <- accnet %>%
    mutate(value =1) %>%
    select(Prot_prot,Genome_genome,value) %>%
    spread(Prot_prot,value, fill = 0) %>% rename(Source = Genome_genome)


  Annotation <- mmseqs$annot  %>% semi_join(accnet %>% select(Prot_prot) %>% distinct(), by = "Prot_prot")
  results <- list(list = accnet %>% rename(Target = Prot_prot,
                                          Source = Genome_genome,
                                          degree = freq) %>% select(Source,Target,degree) ,
                 matrix = accnet.matrix,
                 annot = Annotation %>%
                         rename(ID = Prot_prot) %>%
                         group_by(ID) %>%
                         summarise(Annot = first(Annot)) %>%
                         ungroup()
                 )
  class(results) <- "accnet"
  return(results)
}
#
