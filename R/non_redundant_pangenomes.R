#' Title
#'
#' @param files
#' @param distance
#' @param type
#' @param n_cores
#' @param sketch
#' @param kmer

#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @import data.table
#' @import igraph
#
#' @return
#' @export
#'
#' @examples
non_redundant_pangenomes <- function(files, distance, type = "prot", n_cores,sketch = 1000, kmer = 21)
{


  if(grepl('linux',Sys.getenv("R_PLATFORM"))) ## Linux
  {
    mashPath = system.file("mash",package = "pato")
  }else if(grepl('apple',Sys.getenv("R_PLATFORM"))){ ##MacOS
    mashPath = system.file("mash.macos",package = "pato")
  }else{
    stop("Error, OS not supported.")
  }



  folderName = paste(getwd(),"/",md5(paste(files[,1], sep = "",collapse = "")),"_nr_mash",sep = "",collapse = "")

  if(!dir.exists(folderName))
  {
    dir.create(folderName,)
  }

  write.table(files[,1],paste(folderName,"/input_mash.txt",sep = "",collapse = ""),
              quote = F, col.names = FALSE, row.names = FALSE)
  if(!file.exists(paste(folderName,"/all.msh",sep = "", collapse = ""))){
    if(type == "prot")
    {
      cmd1 <- paste(mashPath," sketch -p ",n_cores," -s ",sketch," -k ",kmer," -l ",folderName,"/input_mash.txt"," -a -o ",folderName,"/all.msh", sep = "", collapse = "")
    }else if(type =="nucl")
    {
      cmd1 <- paste(mashPath," sketch -p ",n_cores," -s ",sketch," -k ",kmer," -l ",folderName,"/input_mash.txt"," -o ",folderName,"/all.msh", sep = "", collapse = "")
    } else{
      stop("Error in type options. Only prot or nucl options are allowed")
    }
    system(cmd1)
  }


  cmd3 <- paste(mashPath," dist -p ",n_cores," -d ",distance," ",folderName,"/all.msh ",folderName,"/all.msh > ",folderName,"/Dist.tab", sep = "", collapse = "")
  system(cmd3)

  mash.table <- data.table::fread(paste(folderName,"/Dist.tab",sep = "",collapse = ""), header = F) %>% as_tibble()

  colnames(mash.table) <- c("Source","Target","Dist","p-value","sketch")

  gr.tmp <- mash.table %>%
    select(Source,Target) %>%
    graph_from_data_frame() %>%
    simplify(., remove.multiple = TRUE, remove.loops = TRUE) %>%
    as.undirected()
  cluster <- components(gr.tmp)

  Nc <- as.numeric(max(cluster$membership))
  cent <- centralization.degree(gr.tmp)
  results = data.frame(Source = as.character(vertex.attributes(gr.tmp)$name),
                       centrality = cent$res,
                       cluster = cluster$membership, stringsAsFactors = F)
  return(results)
}
