#' Individual pangenome clustering from file-list
#'
#' This is an internal function for the workflow of pangenomes
#'
#' @param files a list of files (genomes of the pangenome)
#' @param i The pangenome number
#' @param identity min identitity for homologous cluster
#' @param coverage min coverage  for homologous cluster
#' @param evalue max E-value for homologous cluster
#' @param n_cores number of cores to use
#' @param cov_mode coverage mode
#' @param cluster_mode cluster mode
#'
#' @return A \emph{list} with two tables, the membership of the pangenome, 
#' and the gene/protein frequency.
#'
#' @examples
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @import data.table
#'
pangenomes_from_files_mmseqs <- function(files,i, coverage, identity, evalue, n_cores = 5, cov_mode,cluster_mode)
{

  if(sum(grep("avx2",system("cat /proc/cpuinfo",intern = TRUE),ignore.case = TRUE)))
  {
    mmseqPah <-  system.file("mmseqs.avx2", package = "PATO")
  }else{
    mmseqPah <-  system.file("mmseqs.sse41", package = "PATO")
  }


  members <- data.frame()

  system("rm tmp.fasta*")
  num <- 0


  for(f in files$Source)
  {
    num <-  num+1
    members <- bind_rows(members,data.frame(pangenome = i, file = f, number = num))
    if(grep("gz",f))
    {
      system(paste("zcat ",f," | sed 's/>/>",i,".",num,"|/' >> tmp.fasta",sep ="",collapse = ""))
    }else{
      system(paste("sed 's/>/>",i,".",num,"|/' ",f," >> tmp.fasta",sep ="",collapse = ""))
    }
  }
  cmd <- paste(mmseqPah,
        " easy-linclust tmp.fasta pangenome_",i,
        " tmpDir ",
        " --threads ",n_cores,
        " -e ",evalue,
        " --min-seq-id ",identity,
        " -c ",coverage,
        " --cov-mode ", cov_mode,
        " --cluster-mode ",cluster_mode,
        " -v 0",
        sep = "",collapse = "")


  system(cmd,intern = T)

  cluster_table <- fread(paste("pangenome_",i,"_cluster.tsv", sep = "", collapse = ""),header = F, sep = "\t")
  colnames(cluster_table) <- c("cluster","prot")
  cluster_table$pangenome <-  paste("pangenome_",i, sep = "", collapse = "")
  cluster_table$file <-  paste("pangenome_",i,"_cluster.tsv", sep = "", collapse = "")

  list( table = cluster_table %>%
          group_by(file,pangenome,cluster) %>%
          summarise(freq = n()/nrow(files)) %>%
          ungroup() %>% mutate(pangenome = as.character(i)),
        members = members) %>% return()

}
