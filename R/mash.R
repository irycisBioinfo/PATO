#' MASH distance estimation.
#'
#' MASH (Fast genome and metagenome distance estimation using MinHash)
#' is a fast sequence distance estimator that uses the MinHash
#' algorithm and is designed to work with genomes and metagenomes in the
#' form of assemblies or reads (https://mash.readthedocs.io/). This function
#' is a wrapper to execute \emph{mash} in the background and import to R as a
#' \emph{mash} object.
#'
#' @param file_list Data frame with the full path to the genome files (gene or protein multi-fasta).
#' @param n_cores Number of cores to use.
#' @param sketch Number of sketches to use for distance estimation.
#' @param kmer Kmer size.
#' @param type Type of sequence 'nucl' (nucleotides) or 'prot' (aminoacids)
#'
#' @return A \emph{mash} object
#'
#' @note A \emph{mash} is a list of two element.
#'
#' The first one contains a rectangular and simetric \emph{matrix} with the distances among genomes.
#' As a matrix has genomes as rownames and colnames
#'
#' The second one is a \emph{data.table/data.frame} with all the distancies as list.
#' The table has the columns \emph{c("Source","Target","Dist")}
#' @export
#'
#' @references Mash: fast genome and metagenome distance estimation using MinHash. Ondov BD, Treangen TJ, Melsted P, Mallonee AB, Bergman NH, Koren S, Phillippy AM. Genome Biol. 2016 Jun 20;17(1):132. doi: 10.1186/s13059-016-0997-x.
#' @references Mash Screen: High-throughput sequence containment estimation for genome discovery. Ondov BD, Starrett GJ, Sappington A, Kostic A, Koren S, Buck CB, Phillippy AM. BioRxiv. 2019 Mar. doi: 10.1101/557314
#' @examples
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @import data.table



mash <- function(file_list, n_cores =4, sketch = 1000, kmer = 21, type = "prot")
{


  mashPath = system.file("mash",package = "pato")
  system("rm all.msh input_mash.txt")

  write.table(file_list[,1],"input_mash.txt", quote = F, col.names = FALSE, row.names = FALSE)

  if(type == "prot")
  {
    cmd1 <- paste(mashPath," sketch -p ",n_cores," -s ",sketch," -k ",kmer," -l input_mash.txt"," -a -o all.msh", sep = "", collapse = "")
  }else if(type =="nucl")
  {
    cmd1 <- paste(mashPath," sketch -p ",n_cores," -s ",sketch," -k ",kmer," -l input_mash.txt"," -o all.msh", sep = "", collapse = "")
  } else{
    stop("Error in type options. Only prot or nucl options are allowed")
  }

  system(cmd1)

  cmd3 <- paste(mashPath," dist -p ",n_cores," -t all.msh all.msh > Dist.tab", sep = "", collapse = "")
  system(cmd3)

  mash.matrix <- data.table::fread("Dist.tab", header = T) %>% as_tibble()
  colnames(mash.matrix) <- gsub("#","",colnames(mash.matrix))

  mash.list <- mash.matrix %>%
    mutate(Genome = basename(query)) %>%
    select(-query)%>%
    gather(Target,Dist, -Genome) %>%
    rename(Source = Genome) %>%
    mutate(Target = basename(Target))

  mash.matrix <- mash.matrix %>%
         mutate(Genome = basename(query)) %>%
         select(-query) %>%
         column_to_rownames("Genome") %>%
         as.matrix()

  colnames(mash.matrix) <- rownames(mash.matrix)
  results <- list(matrix = mash.matrix, table = mash.list)
  class(results) <- "mash"
  return(results)

}
