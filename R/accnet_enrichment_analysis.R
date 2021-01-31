#' Enrichment analysis of accessory genome
#'
#' Perform a statistical test (Hypergeometric test: \code{phyper}) of the genes/proteins.
#' It takes the cluster definition from the \code{cluster} table and compute
#' the statistical test assuming the frequency of the gene/protein in the cluster vs
#' the frequency of the gene/protein in the population.
#'
#' @param data An \code{accnet} object
#' @param cluster Dataframe with two colums (Source,Cluster)
#' @param max_pvalue Maximun p-value of retrived data.
#' @param padj_method Method adjust the p-value for multiple comparisons:
#' "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none" (see p.adjust)
#'
#' @return Return a data.frame with the columns:
#' \itemize{
#'   \item Target: Protein
#'   \item Source: Genome
#'   \item Cluster: Cluster of genome
#'   \item perClusterFreq: Percentage (\%) of the protein in this cluster
#'   \item ClusterFreq: Frequency [0-1] of the protein in this cluster
#'   \item ClusterGenomeSize: Number of genomes in the cluster
#'   \item perTotalFreq: Percentage (\%) of the protein in the population
#'   \item TotalFreq: Frequency [0-1] of the protein in the population
#'   \item OdsRatio: Ods ration of the protein
#'   \item pvalue: p-value of the hypergeometric test (see \code{phyper})
#'   \item padj: Adjusted p-value using \code{p.adjust}
#'   \item AccnetGenomeSize: Total number of genomes
#'   \item AccnetProteinSize: Total number of proteins
#'   \item Annot: Funtional annotation of the protein.
#'   }
#'
#' @seealso \code{\link[stats]{phyper}}
#' @seealso \code{\link[stats]{p.adjust}}
#' @seealso \code{\link{accnet}}
#' @export
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
accnet_enrichment_analysis <- function(data,cluster, padj_method = "BY")
{

  if(!is(data,"accnet"))
  {
    stop("Data must be an 'accnet' object")
  }



  colnames(cluster) <- c("Source","Cluster")
  accnet <- full_join(cluster,data$list)


  accnet <-accnet %>% filter(!is.na(Target)) %>% filter(!is.na(Source)) %>% select(-degree)





  tmp = accnet %>% select(Source) %>% distinct() %>% count()
  accnet$AccnetGenomeSize = tmp$n

  tmp = accnet %>% ungroup() %>% select(Target) %>% distinct() %>% count()
  accnet$AccnetProteinSize = tmp$n


  accnet <-accnet %>%
    select(Target, Source) %>% distinct() %>% group_by(Target) %>% mutate(TotalFreq = n()) %>%
    ungroup() %>% group_by(Source) %>% mutate(Degree = n()) %>%
    full_join(accnet, by = c("Target", "Source"))

  accnet <-accnet %>% select(Cluster, Source) %>% distinct() %>%
    group_by(Cluster) %>% summarise(ClusterGenomeSize = n()) %>%
    full_join(accnet, by = "Cluster")

  accnet <-accnet %>% select(Cluster, Target) %>% distinct() %>%
    group_by(Cluster) %>% summarise(ClusterProteinSize = n()) %>%
    full_join(accnet, by = "Cluster")

  accnet <-accnet %>% group_by(Cluster, Target) %>% mutate(ClusterFreq = n())

  accnet <-accnet %>%
    mutate(
      perClusterFreq = ClusterFreq / ClusterGenomeSize,
      perTotalFreq = TotalFreq / AccnetGenomeSize
    )

  accnet <-accnet %>%
    select(Source,
      Target,
      Cluster,
      TotalFreq,
      ClusterFreq,
      AccnetProteinSize,
      AccnetGenomeSize,
      ClusterGenomeSize,
      perClusterFreq,
      perTotalFreq
    ) %>%
    distinct() %>%
    mutate(
      OddsRatio =  perClusterFreq / perTotalFreq ,
      pvalue = phyper(
        ClusterFreq,
        TotalFreq,
        AccnetGenomeSize - TotalFreq,
        ClusterGenomeSize,
        lower.tail = FALSE
      )
    )

  accnet <-accnet %>% select(Cluster, Target, pvalue) %>%
    distinct() %>% group_by(Cluster) %>% mutate(padj = p.adjust(pvalue, method = padj_method)) %>%
    full_join(accnet, by = c("Cluster", "Target", "pvalue"))

  accnet <-accnet %>%  select(
    Target,
    Source,
    Cluster,
    perClusterFreq,
    ClusterFreq,
    ClusterGenomeSize,
    perTotalFreq,
    TotalFreq,
    OddsRatio,
    pvalue,
    padj,
    AccnetGenomeSize,
    AccnetProteinSize
  ) %>% full_join(data$annot, by = c("Target"="ID"))

  result <- accnet %>% ungroup() %>% as_tibble()
  class(result) <- append("accnet_enr",class(result))
  return(result)
}
