#' Creates Pangenomes from file list
#'
#' This function perform a clustering using a designated maximun distance among samples
#' and creates a set of pangenomes using it.
#' Pan-genomes can be filtered by genome frequency (i.e. the number of genomes that belong to
#' the pangenome). \emph{pangenomes_from_list()} builds a new \emph{accnet} object that relates
#' the pangenome with the genes/proteins. The user can filter the proteins that are included
#' by the presence frequency in the set.
#'
#' @param files A data.frame of one column whit the path of the files.
#' @param distance Maximun distance (mash distance) among genomes in each pangenome. (exlude \emph{cluster})
#' @param cluster (optional, exclude distance). Data.frame of two columns (file, cluster) describing the clutering of the files
#' @param coverage Minimun coverage (length) to cluster.
#' @param identity Minimun Identity.
#' @param evalue Maximun Evalue.
#' @param n_cores Number of cores to use.
#' @param cluster_mode 	Cluster mode:\itemize{
#' \item{0: Setcover}
#' \item{1: connected component}
#' \item{2: Greedy clustering by sequence length}
#' \item{3: Greedy clustering by sequence length (low mem)}
#' }
#' @param cov_mode Coverage mode:\itemize{
#' \item 0: Coverage of query and target
#' \item 1: Coverage of target
#' \item 2: coverage of query
#' \item 3: target seq. length needs be at least x% of query length
#' \item 4: query seq. length needs
#' }
#' @param file_type Type of fasta file (prot or nucl)
#' @param min_pange_size Minimun number of genomes per pan-genome
#' @param min_prot_freq Mininum frequency of a protein (minimun number of pan-genomes in where is present)
#'
#' @return An \emph{accnet} with an extra membership table.
#' @export
#'
#' @examples
pangenomes_from_files <- function(files, min_pange_size = 10, min_prot_freq = 2, file_type = 'prot', distance, cluster, coverage = 0.8, identity = 0.8, evalue = 1e-6, n_cores = 5, cov_mode = 0, cluster_mode = 0)
{

  if(!missing(distance))
  {
    cl <- cluster_files_from_distance(files, file_type, distance, n_cores)
  }else if (!missing(cluster))
  {
    cl <-  cluster
  } else{
    stop("Error: distance or cluster must be provided")
  }

  cluster <- as_tibble(cl)
  colnames(cluster) <- c("Source","Cluster")

  cluster <- cluster %>% group_by(Cluster) %>% mutate(Size = n()) %>% filter(Size >= min_pange_size)

  cluster$Cluster <- as.factor(cluster$Cluster)

  system("rm pangenome_*_*")
  results <- data.frame()
  members <- data.frame()
  for (lv in levels(cluster$Cluster))
  {
    print(lv)
    tmp <- pangenomes_from_files_mmseqs(files = cluster %>% filter(Cluster ==lv),
                                 lv,
                                 coverage,
                                 identity,
                                 evalue,
                                 n_cores,
                                 cov_mode,
                                 cluster_mode)

    results <- bind_rows(results,tmp$table)

    members <- bind_rows(members,tmp$members)

  }
  mm <- pangenomes_mmseqs(results %>%
                            select(file) %>% mutate(file = gsub("_cluster.tsv","_rep_seq.fasta",file)) %>%
                            distinct(), coverage, identity, evalue, n_cores, cov_mode, cluster_mode)
  results <- results %>%
    mutate(Genome_genome = gsub("_cluster.tsv","_rep_seq.fasta",file)) %>%
    rename(Genome_prot = cluster) %>%
    inner_join(mm$table) %>% as_tibble() %>% group_by(Prot_prot) %>% mutate(degree = n()) %>% filter(degree >= min_prot_freq)

  results <- list(list = results %>% rename(Target = Prot_prot, Source = Genome_genome, Weight = freq) %>%
                    select(Source,Target,Weight),
                   matrix = results %>% rename(Target = Prot_prot,
                                               Source = Genome_genome,
                                               Weight = freq) %>%
                                         group_by(Source,Target) %>%
                                         summarise(Weight = max(Weight)) %>%
                                         spread(Target,Weight, fill = 0),
                   annot = mm$annot %>% as_tibble() %>%
                                         semi_join(results %>% select(Prot_prot) %>% distinct(), by = "Prot_prot") %>%
                                         rename(ID = Prot_prot) %>%
                                         group_by(ID) %>%
                                         summarise(Annot = first(Annot)) %>%
                                         ungroup(),
                  members = members %>% rename(path =file) %>% mutate(file = basename(path))
                )

  class(results) = c("accnet","pangenome")
  return(results)




}
