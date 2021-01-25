#' Export accnet or igraph object to GraphML format
#'
#' This function export a network (\emph{accnet} or \emph{igraph}) to universal
#' format graphml.
#'
#' @param data \emph{accnet} or \emph{igraph} object
#' @param file Filename for network
#'
#' @export
#'
#' @note Big networks (especially accnet networks) can not be visualized in R
#' enviroment. We recomend to use external software such as Gephi or Cytoscape.
#' Gephi can manage bigger networks and the ForceAtlas2 layout can arrange
#' accnet networks up to 1.000 genomes. If \emph{accnet} object comes from
#' \emph{accnet_with_padj()} function, the graph files will include the
#' edges-weight proportional to the p-value that associate the protein with
#' the genome.
#'
#' @import igraph
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @import igraph
#' @import data.table
export_to_gephi <- function(data, file, cluster)
{
  if(is(data,"accnet"))
  {
    file.net = paste(file,".net.tsv",collapse = "",sep="")
    file.annot = paste(file,".table.tsv",collapse = "",sep="")

    if("Weight" %in% colnames(data$list))
    {
      data$list %>% select(Source,Target,edge_Weight = Weight) %>%
        distinct() %>%
        write.table(file =file.net,sep = "\t",row.names = FALSE, col.names = TRUE, quote = FALSE)
    }else{
      data$list %>% select(Source,Target) %>% distinct() %>%
        write.table(file =file.net,sep = "\t",row.names = FALSE, col.names = TRUE, quote = FALSE)
    }



    annotation <- data$annot %>% ungroup() %>%
      mutate(Type = "Protein") %>%
      bind_rows(.,data$list %>%
                  select(ID = Source) %>% distinct() %>% mutate(Type ="Genome"))


    if(!missing(cluster))
    {
       annotation <- annotation %>% ungroup() %>%
         left_join(cluster, by  = c("ID" = "Source")) %>%
         replace_na(list(Cluster =0))
    }

    annotation %>% ungroup() %>% distinct() %>%
      write.table(file =file.annot,sep = "\t",row.names = FALSE, col.names = TRUE, quote = FALSE)

  }else if(is(data,"igraph"))
  {
    file.net = paste(file,".net.tsv",collapse = "",sep="")
    file.annot = paste(file,".table.tsv",collapse = "",sep="")
    igraph::as_data_frame(data) %>% rename(Source = from, Target = to) %>% ungroup() %>%
      write.table(file =file.net,sep = "\t",row.names = FALSE, col.names = TRUE, quote = FALSE)


  }else{
    stop("Error: 'data' must be accnet or igraph object")
  }

}
