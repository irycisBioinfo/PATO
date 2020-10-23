#' Individual pangenome clustering from file-list
#'
#' This is an internal function for the workflow of pangenomes
#'
#' @param file_list a list of file_list (genomes of the pangenome)
#' @param i The pangenome number
#' @param identity min identitity for homologous cluster
#' @param coverage min coverage  for homologous cluster
#' @param evalue max E-value for homologous cluster
#' @param n_cores number of cores to use
#' @param cov_mode coverage mode
#' @param cluster_mode cluster mode
#' @param folder
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
pangenomes_from_files_mmseqs <- function(file_list,i, coverage, identity, evalue, n_cores, cov_mode,cluster_mode,folder)
{

  if(missing(n_cores))
  {
    n_cores = parallel::detectCores()-1
  }

  if(grepl('linux',Sys.getenv("R_PLATFORM"))) ## Linux
  {
    proc_cpu = readLines("/proc/cpuinfo")

    if(sum(grep("avx2",proc_cpu,ignore.case = TRUE)))
    {
      mmseqPath = system.file("mmseqs.avx2", package = "pato")
    }else{
      mmseqPath = system.file("mmseqs.sse41", package = "pato")
    }
  }else if(grepl('apple',Sys.getenv("R_PLATFORM"))){ ##MacOS


    if(grepl("AVX2",system("sysctl -a | grep 'AVX2'", intern = T)))
    {
      mmseqPath = system.file("mmseqs.macos.avx2", package = "pato")
    }else{
      mmseqPath = system.file("mmseqs.macos.sse41", package = "pato")
    }
  }else{
    stop("Error, OS not supported.")
  }



  members <- data.frame()


  num <- 0




  for(f in file_list[,1])
  {

    num <-  num+1
    members <- bind_rows(members,data.frame(pangenome = i, file = f, number = num))
    if(grepl("gz",f))
    {
      print(paste("zcat ",f," | sed 's/>/>",i,".",num,"|/' >> ",folder,"/tmp.fasta",sep ="",collapse = ""))
      system(paste("zcat ",f," | sed 's/>/>",i,".",num,"|/' >> ",folder,"/tmp.fasta",sep ="",collapse = ""))
    }else{
      print(paste("sed 's/>/>",i,".",num,"|/' ",f," >> ",folder,"/tmp.fasta",sep ="",collapse = ""))
      system(paste("sed 's/>/>",i,".",num,"|/' ",f," >> ",folder,"/tmp.fasta",sep ="",collapse = ""))
    }
  }

  cmd <- paste(mmseqPath,
               " easy-linclust ",folder,"/tmp.fasta ",folder,"/pangenome_",i,
               " ",folder,"/tmpDir ",
               " --threads ",n_cores,
               " -e ",evalue,
               " --min-seq-id ",identity,
               " -c ",coverage,
               " --cov-mode ", cov_mode,
               " --cluster-mode ",cluster_mode,
               " -v 0",
               sep = "",collapse = "")

  print(cmd)
  system(cmd)

  cluster_table <- fread(paste(folder,"/pangenome_",i,"_cluster.tsv", sep = "", collapse = ""),header = F, sep = "\t")
  colnames(cluster_table) <- c("cluster","prot")
  cluster_table$pangenome <-  paste("pangenome_",i, sep = "", collapse = "")
  cluster_table$file <-  paste("pangenome_",i,"_cluster.tsv", sep = "", collapse = "")

  list( table = cluster_table %>%
          group_by(file,pangenome,cluster) %>%
          summarise(freq = n()/nrow(file_list)) %>%
          ungroup() %>% mutate(pangenome = as.character(i)),
        members = members) %>% return()

}
