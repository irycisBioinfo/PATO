#' Non Redundant set of samples using hierarchical search.
#'
#' Creates a non-redundant sub-set of samples. The function accept an
#' \emph{accnet} or \emph{mash} object and returns a \emph{nr_list} object.
#' The sub-set can be created using a certain \emph{distance} (mash o jaccard distance),
#' a specific \emph{number} of elements or a \emph{fraction} of the whole sample set.
#' The function performes an iterative seach so sometimes the exact number of returned elements
#' could be not the same that the specified in the input. This difference can be defined
#' with the \emph{threshold} parameter. This function perform a hierarchical search. For small
#' datasets (<2000) \emph{non_redundant()} could be faster. This function comsume less memory than
#' the original so it is better for large datasets.
#'
#' @param data An \emph{accnet} or \emph{mash} object.
#' @param number The number of non-redundant samples.
#' @param fraction The fraction of the whole set of non-redundant samples.
#' @param distance Minimun distance among samples.
#' @param tolerance Percentage of error between the input number and the final number of samples.
#' @param max_iter Maximun number of search iterations.
#' @param fast If fast is TRUE the clustering process uses "components" in other case use "fast_greddy".
#' @param partitions Number of partitions to hiaerarchical search.
#'
#' @return \emph{nr_list} object
#' @export
#'
#' @seealso \code{\link{extract_non_redundant}}
#'
#' @examples
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @import igraph
#' @import data.table
#' @import parallelDist
#'
non_redundant_hier <- function(data, number, fraction, distance, tolerance = 0.05, partitions = 10,max_iter = 10000, fast =FALSE)
{
  if(is(data,"accnet"))
  {
    m.matrix <-  data$matrix %>%
      as_tibble() %>%
      column_to_rownames("Source") %>%
      as.matrix() %>%
      parallelDist(., method = "binary") %>%
      as.matrix()
    m.list <- m.matrix %>%
      rownames_to_column("Source") %>%
      spread(Target,Dist,-Source)

  }else if (is(data,"mash"))
  {
    m.matrix <-  data$matrix
    m.list <- data$table
  }else {
    print("Error: data must be accnet or mash object")
  }

### With distance

  if (!missing(distance)){
    gr.tmp <- m.list %>% filter(Dist <= distance) %>% graph_from_data_frame(.,directed = FALSE)
    gr.tmp <- simplify(gr.tmp, remove.multiple = TRUE, remove.loops = TRUE) %>% as.undirected()

    if(fast ==TRUE) {
      cluster <- components(gr.tmp)
    }else{
      cluster <- cluster_louvain(gr.tmp)
    }

    Nc <- as.numeric(max(cluster$membership))

    if(is.infinite(Nc))
    {
      cluster <- components(gr.tmp)
      Nc <- as.numeric(max(cluster$membership))
    }


    cent <- centralization.degree(gr.tmp)
    results <- data.frame(Source = as.character(vertex.attributes(gr.tmp)$name),
                          centrality = cent$res,
                          cluster = cluster$membership)
    class(results) <- "nr_list"
    return(results)
  }

## With numbers

  if(!missing(number))
  {

    min <- 0
    max <- 1
    Th <- 0.05
  }
  else if (!missing(fraction))
  {
    number <-  m.list %>% select(Source) %>% distinct() %>% count()
    number <- number[1] * fraction
    min <- 0
    max <- 1
    Th <- 0.05
  }else{
    print("Error: number, fraction or distance must be provided")
  }


  steps <- round(seq(1,nrow(m.matrix), length.out = partitions))
#####
  table.tmp <- data.frame()
  for(i in 1:(partitions-1))   ### First Iteration
  {

    m.tmp <-  m.matrix[steps[i]:steps[i+1],steps[i]:steps[i+1]]
    if(fast==TRUE){
      cm <- graph_from_adjacency_matrix(m.tmp < Th, mode = "upper") %>% as.undirected()%>% components()
    }else{
      cm <- graph_from_adjacency_matrix(m.tmp < Th, mode = "upper") %>% as.undirected() %>% cluster_louvain()
    }
    table.tmp <-  cm$membership %>%
      as.data.frame() %>%
      rownames_to_column("Source") %>%
      rename(cluster = 2) %>%
      group_by(cluster) %>%
      summarise(Source = first(Source)) %>% bind_rows(table.tmp)
  }


  gr.tmp <- m.list %>% as_tibble() %>%
    semi_join(table.tmp %>%
                select(Source), by ="Source") %>%
    semi_join(table.tmp %>%
                rename(Target = Source), by ="Target") %>%
    filter(Dist <= Th) %>%
    graph_from_data_frame()

  if(fast == TRUE)
  {
    cm.all <- gr.tmp %>% as.undirected() %>% components()
  }else{
    cm.all <- gr.tmp %>% as.undirected() %>% cluster_louvain()
  }

  Nc <- as.numeric(max(cm.all$membership))

###### END first Iteration

  count <-0

  while (count <= max_iter)
  {
    print(Th)
    if (abs(1-(Nc/number)) < tolerance){
      gr.tmp <- m.list %>% filter(Dist <= Th) %>% graph_from_data_frame(.,directed = FALSE)
      cluster <- simplify(gr.tmp, remove.multiple = TRUE, remove.loops = TRUE) %>%
        as.undirected() %>%
        cluster_louvain()
      cent <- centralization.degree(gr.tmp)
      results <- data.frame(Source = as.character(vertex.attributes(gr.tmp)$name),
                            centrality = cent$res,
                            cluster = cluster$membership)
      class(results) <- "nr_list"
      return(results)
    }else if (Nc > number){
      min <- Th
      Th <-(max + Th)/2
    }else if (Nc < number){
      max <- Th
      Th <- (Th - min)/2
    }

    table.tmp <- data.frame()
    for(i in 1:(partitions-1))
    {
      m.tmp <-  m.matrix[steps[i]:steps[i+1],steps[i]:steps[i+1]]
      if(fast==TRUE){
        cm <-  graph_from_adjacency_matrix(m.tmp < Th, mode = "upper") %>%
          as.undirected() %>%
          components()
      }else{
        cm <- graph_from_adjacency_matrix(m.tmp < Th, mode = "upper") %>%
          as.undirected() %>%
          cluster_louvain()
      }
      table.tmp <-  cm$membership %>%
        as.data.frame() %>%
        rownames_to_column("Source") %>%
        rename(cluster = 2) %>%
        group_by(cluster) %>%
        summarise(Source = first(Source)) %>%
        bind_rows(table.tmp)
    }
    gr.tmp <- m.list %>%
      semi_join(table.tmp %>%
                  select(Source), by = "Source") %>%
      semi_join(table.tmp %>%
                    rename(Target = Source), by = "Target") %>%
                    filter(Dist <= Th) %>%
                    graph_from_data_frame()

    if(fast == TRUE)
    {
      cm.all <- gr.tmp %>% as.undirected() %>% components()
    }else{
      cm.all <- gr.tmp %>% as.undirected() %>% cluster_louvain()
    }

    Nc <- as.numeric(max(cm.all$membership))
    print(paste("Th =",Th,"  Nc=",Nc,sep = "",collapse = ""))

  }
  stop("Max iter reached")


}

