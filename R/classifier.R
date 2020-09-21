#' Classifier
#'
#' Classifier take a list of genome files (nucleotide o protein) and identify
#' the most similar specie to each file.
#' Classsifier uses all reference and representative genomes from NCBI Refseq
#' database and search, using mash, the best hit for each genome file in the input list.
#'
#' @param file_list Data frame with the full path to the genome files (gene or protein multi-fasta).
#' @param re_use You can re-use the sketch step if you have executed a previous \emph{mash()} function
#' @param n_cores Number of cores to use.
#' @param type Type of sequence 'nucl' (nucleotides) or 'prot' (aminoacids)
#' @param max_dist Maximun distance to report (1-Average Nucleotide Identity). Usually all species have 0.05 distance among all each memebers
#'
#' @return Classifier returns a data.frame with the best hit for each input genome.
#' @export
#'
#' @examples
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @import data.table

classifier <- function(file_list, re_use = TRUE, n_cores, type ="nucl", max_dist = 0.06)
{
  if(missing(n_cores))
  {
    n_cores <- detectCores()-1
  }

  if(type == "prot")
  {
    reference <- system.file("classifier/ReferenceProt.msh",package = "pato")
  }else if(type == "nucl")
  {
    reference <- system.file("classifier/ReferenceNucl.msh",package = "pato")
  }else{
    stop("Error in type options. Only prot or nucl options are allowed")
  }


  mashPath <- system.file("mash",package = "pato")
  if(!re_use)
  {
    system("rm all.msh input_mash.txt")
    write.table(file_list[,1],"input_mash.txt", quote = F, col.names = FALSE, row.names = FALSE)

    if(type == "prot")
    {
      cmd1 <- paste(mashPath," sketch -p ",n_cores," -l input_mash.txt"," -a -o all.msh", sep = "", collapse = "")
    }else if(type =="nucl")
    {
      cmd1 <- paste(mashPath," sketch -p ",n_cores," -l input_mash.txt"," -o all.msh", sep = "", collapse = "")
    } else{
      stop("Error in type options. Only prot or nucl options are allowed")
    }
    system(cmd1)
  }


  cmd3 <- paste(mashPath," dist -p ",n_cores," -d ",max_dist," ",reference," all.msh > classification.tab", sep = "", collapse = "")
  system(cmd3)

  class.table <- data.table::fread("classification.tab", header = F) %>% as_tibble()
  colnames(class.table) <- c("Source","Target","Dist","pvalue","sketch")


  header <- data.table::fread(system.file("classifier/headers.tsv",package = "pato"), sep = "\t", header = T, stringsAsFactors = F)



  class.table %>%
    mutate(Target = basename(Target)) %>%
    separate(Source,c("kk","acc1","acc2"), sep = "_") %>%
    unite(assembly_accession,kk,acc1,sep = "_") %>%
    mutate(Dist = 1-Dist) %>%
    left_join(header) %>%
    select(Target,Similarity = Dist,organism_name,species_taxid,assembly_accession,refseq_category) %>%
    group_by(Target) %>% slice_max(order_by = Similarity,n = 1) %>% ungroup() %>%
    return()


}
