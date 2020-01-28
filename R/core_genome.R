#' Core-Genome Alignment
#'
#' Find and creates a core-genome alignment. Unlike \emph{core_plots()}
#' this function find the hard core-genome (genes presence in 100% of genomes
#' and without repetitions). The function takes a \emph{mmseqs()} output, so the deffinition of
#' the paralogous genes of the core-genome (similarity, coverage and/or e-value)
#' depends on the \emph{mmseqs()} parameters.
#'
#' The function perform a pseudo-msa per each paralog using the function
#' \emph{result2msa} of \emph{mmseqs2}.This approach is much faster than
#' classical MSA (clutal, mafft or muscle) but is less accurate. Taking into account
#' that most of the phylogenetic inference software only takes variant columns with
#' no insertions or deletion, there are not to many difference in the final phylogenetic trees.
#'
#' \emph{core_genome()} can build a core-genome alingment of thusands of genomes in minutes.
#'
#' @param data An \emph{mmseqs} object
#'
#' @return A core_genome object (a data.frame with two columns: fasta header and sequence)
#' @export
#'
#' @examples
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @import data.table
#'
core_genome <- function(data)
{
  if(!is(data,"mmseq"))
  {
    stop("Error, 'data' must be a mmseq object")
  }

  if(sum(grep("avx2",system("cat /proc/cpuinfo",intern = TRUE),ignore.case = TRUE)))
  {
    mmseqPah = system.file("mmseqs.avx2", package = "PATO")
  }else{
    mmseqPah = system.file("mmseqs.sse41", package = "PATO")
  }

  if(file.exists("all.mmseq") & file.exists("all.cluster.index"))
  {

    nGenomes  <-  data$table %>% distinct(Genome_genome) %>% count()
    table <-  data$table %>% as_tibble() %>% group_by(Prot_prot) %>% mutate(Nprots = n()) %>%
      group_by(Prot_prot,Nprots) %>% summarise(Ngenomes = n_distinct(Genome_genome)) %>% ungroup() %>%
      filter(Nprots == nGenomes$n & Ngenomes == nGenomes$n) %>% ungroup()

    lookup <- data.table::fread("all.mmseq.lookup") %>% as_tibble()
    colnames(lookup) = c("ID","head","value")

    core_table <-  data$table %>%
      semi_join(table) %>%
      unite(head, Prot_genome, Prot_prot, sep = "#", remove = FALSE) %>%
      select(head) %>%
      distinct() %>% inner_join(lookup) %>%
      select(ID) %>% write.table("IDS.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

    system("rm -r f_core core*")
    system("mkdir f_core")

    paste(mmseqPah," createsubdb IDS.txt all.cluster all.cluster.subset -v 0", sep = "", collapse = "") %>%
      system(.,intern = F,ignore.stdout = T, ignore.stderr = T)
    paste(mmseqPah," result2msa all.mmseq all.mmseq all.cluster.subset all.msa --summarize --max-seq-id 1 --diff ",nGenomes$n+1," --qsc -50 --omit-consensus -v 0", sep = "", collapse = "") %>%
      system(.,intern = F,ignore.stdout = T, ignore.stderr = T)
    system("strings all.msa > all.msa.flat")
    system("csplit -z all.msa.flat /#cl-/ {*} -f core --suppress-matched")

    system("mv core* f_core")

    seqs <-  data.frame()
    for (i in(dir("./f_core")))
    {
      tmp <- readLines(paste("./f_core/",i,sep = "",collapse = ""))
      tmp2 <- data.frame(Head = tmp[seq(1,length(tmp)-1,2)],Seq = tmp[seq(2,length(tmp),2)])
      seqs <- bind_rows(seqs,tmp2)
    }
    print(colnames(seqs))
    seqs <- seqs %>%
      separate(Head,c("Genomes","Prot"),sep="#") %>%
      group_by(Genomes) %>%
      summarise(Seq = paste(Seq,sep = "",collapse = ""))

    print("Number of hard core-genes (100% presence):")
    print(nrow(table))
    class(seqs) <- "core_genome"
    return(seqs)



  } else{
    stop("You must execute mmseqs() before core_genome()")
  }






}
