#' Internal Function of the workflow pangenomes_from_file
#'
#' @param files data.frame with the path of the selected genomes
#' @param file_type Type of fasta files 'prot' or 'nucl'
#' @param distance threshold distances of the clustering proccess
#' @param n_cores Num of cores to use with MASH
#' @param folder
#'
#' @return A data frame with the membership of each genome.
#'
#'
#' @import dplyr
#' @import dtplyr
#' @import tidyr
#' @import tibble
#' @importFrom data.table fread
#' @import igraph
#'
cluster_files_from_distance <- function(files, file_type, distance, n_cores, folder)
{


  if(grepl('linux',Sys.getenv("R_PLATFORM"))) ## Linux
  {
    mashPath = system.file("mash",package = "pato")
  }else if(grepl('apple',Sys.getenv("R_PLATFORM"))){ ##MacOS
    mashPath = system.file("mash.macos",package = "pato")
  }else{
    stop("Error, OS not supported.")
  }





  write.table(files[,1],"input_mash.txt", quote = F, col.names = FALSE, row.names = FALSE)


  if(!file.exists("all.msh"))
  {
    if(file_type == "prot")
    {

      cmd1 <- paste(mashPath," sketch -p ",n_cores," -l input_mash.txt -a -o ",folder,"/all.msh", sep = "", collapse = "")


    }else if(file_type =="nucl")
    {
      cmd1 <- paste(mashPath," sketch -p ",n_cores," -l input_mash.txt -o ",folder,"/all.msh", sep = "", collapse = "")
    } else{
      stop("Error in type options. Only prot or nucl options are allowed")
    }
    print(cmd1)
    system(cmd1)
  }




  if(!file.exists("pangeDist.tab"))
  {
    cmd3 <- paste(mashPath," dist -p ",n_cores," -d ",distance," ",folder,"/all.msh ",folder,"/all.msh > ",folder,"/pangeDist.tab", sep = "", collapse = "")
    system(cmd3)
  }

  table <- fread(paste(folder,"/pangeDist.tab",sep = "",collapse = ""), header = F) %>% as_tibble()
  #table <- vroom::vroom("pangeDist.tab", col_names = F)
  colnames(table) <-  c("from","to","weigth","pvalue","sketch")

  cl <- table %>%
    select(from,to,weigth) %>%
    mutate(weigth = 1-weigth) %>%
    graph_from_data_frame() %>%
    simplify(remove.multiple = T, remove.loops = T) %>%
    as.undirected() %>% cluster_louvain()

  return(data.frame(Source = cl$names, Cluster = cl$membership))
}
