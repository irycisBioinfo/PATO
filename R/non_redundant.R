#' Non Redundant set of samples.
#'
#' Creates a non-redundant sub-set of samples. The function accept an
#' \emph{accnet} or \emph{mash} object and returns a \emph{nr_list} object.
#' The sub-set can be created using a certain \emph{distance} (mash o jaccard distance),
#' a specific \emph{number} of elements or a \emph{fraction} of the whole sample set.
#' The function performes an iterative seach so sometimes the exact number of returned elements
#' could be not the same that the specified in the input. This difference can be defined
#' with the \emph{threshold} parameter.
#'
#' @param data An \emph{accnet} or \emph{mash} object
#' @param number The number of non-redundant samples.
#' @param fraction The fraction of the whole set of non-redundant samples.
#' @param distance Minimun distance among samples
#' @param tolerance Percentage of error between the input number and the final number of samples.
#' @param max_iter Maximun number of search iterations.
#' @param fast If fast is TRUE the clustering process uses "components" in other case use "louvain"
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
non_redundant <- function(data, number, fraction, distance, tolerance = 0.05, max_iter = 10000, fast =TRUE, snps)
{

 if(is(data,"accnet"))
 {
   m.list = data$matrix %>%
     as_tibble() %>%
     column_to_rownames("Source") %>%
     as.matrix() %>%
     parallelDist(., method = "binary") %>%
     as.matrix() %>%
     as.data.frame() %>%
     rownames_to_column("Source") %>%
     gather(Target,Dist, -Source)

 }else if (is(data,"mash")){
   m.list = data$table
 }else if (is(data,"matrix")){
   if(!missing(snps))
   {
     m.list = data %>% as_tibble() %>% rownames_to_column("Source") %>% gather(Target,Dist, -Source)
     gr.tmp <- m.list %>% filter(Dist <= snps) %>% graph_from_data_frame(.,directed = FALSE)
     gr.tmp <- simplify(gr.tmp, remove.multiple = TRUE, remove.loops = TRUE) %>% as.undirected()

     if(fast ==TRUE) {
       cluster <- components(gr.tmp)
     }else{
       cluster <- cluster_louvain(gr.tmp)
     }

     Nc <- as.numeric(max(cluster$membership))
     cent <- centralization.degree(gr.tmp)
     results = data.frame(Source = as.character(vertex.attributes(gr.tmp)$name),
                          centrality = cent$res,
                          cluster = cluster$membership, stringsAsFactors = F)
     class(results) <- append(class(results),"nr_list")
     return(results)
   }else{
     m.list = data %>%
       as_tibble() %>%
       rownames_to_column("Source") %>%
       gather(Target,Dist, -Source)
   }
 }else {
   print("Error: data must be accnet or mash object or snps matrix")
 }

  if(!missing(number))
  {
    min <- min(m.list$Dist)
    max <- max(m.list$Dist)
    Th <- mean(m.list$Dist)


    gr.tmp <- m.list %>% filter(Dist <= Th) %>% graph_from_data_frame(.,directed = FALSE)
    gr.tmp <- simplify(gr.tmp, remove.multiple = TRUE, remove.loops = TRUE) %>% as.undirected()

    if(fast ==TRUE)
    {
      cluster <- components(gr.tmp)
    }else{
      cluster <- cluster_louvain(gr.tmp)
    }

    Nc <- as.numeric(max(cluster$membership))
    count =0

    while (count <= max_iter)
    {
      if (abs(1-(Nc/number)) < tolerance){
        cent <- centralization.degree(gr.tmp)
        results = data.frame(Source = as.character(vertex.attributes(gr.tmp)$name),
                             centrality = cent$res,
                             cluster = cluster$membership, stringsAsFactors = F)
        class(results) <- "nr_list"
        return(results)
      }else if (Nc > number){
        min <- Th
        Th <-(max + Th)/2
      }else if (Nc < number){
        max <- Th
        Th <- (Th - min)/2
      }


      gr.tmp <- m.list %>% filter(Dist <= Th) %>% graph_from_data_frame(.,directed = FALSE)
      gr.tmp <- simplify(gr.tmp, remove.multiple = TRUE, remove.loops = TRUE) %>% as.undirected()

      if(fast ==TRUE) {
        cluster <- components(gr.tmp)
      }else{
          cluster <- cluster_louvain(gr.tmp)
      }
      Nc <- as.numeric(max(cluster$membership))
      count <- count+1

    }
    stop("Max iter reached")

  }else if (!missing(fraction)){
    number = m.list %>% select(Source) %>% distinct() %>% count()
    number = number[1] * fraction

    min <- min(m.list$Dist)
    max <- max(m.list$Dist)
    Th <- mean(m.list$Dist)

    gr.tmp <- m.list %>% filter(Dist <= Th) %>% graph_from_data_frame(.,directed = FALSE)
    gr.tmp <- simplify(gr.tmp, remove.multiple = TRUE, remove.loops = TRUE) %>% as.undirected()

    if(fast ==TRUE) {
      cluster <- components(gr.tmp)
    }else{
        cluster <- cluster_louvain(gr.tmp)
    }

    Nc <- as.numeric(max(cluster$membership))
    count =0

    while (count <= max_iter)
    {

      if (abs(1-(Nc/number)) < tolerance){
        cent <- centralization.degree(gr.tmp)
        results = data.frame(Source = as.character(vertex.attributes(gr.tmp)$name),
                             centrality = cent$res,
                             cluster = cluster$membership,
                             stringsAsFactors = F)
        class(results) <- "nr_list"
        return(results)
      }else if (Nc > number){
        min <- Th
        Th <-(max + Th)/2
      }else if (Nc < number){
        max <- Th
        Th <- (Th - min)/2
      }

      gr.tmp <- m.list %>% filter(Dist <= Th) %>% graph_from_data_frame(.,directed = FALSE)
      gr.tmp <- simplify(gr.tmp, remove.multiple = TRUE, remove.loops = TRUE) %>% as.undirected()

      if(fast ==TRUE) {
        cluster <- components(gr.tmp)
      }else{
          cluster <- cluster_louvain(gr.tmp)
      }
      Nc <- as.numeric(max(cluster$membership))
      count <- count+1

    }
    stop("Max iter reached")

  }else if (!missing(distance)){
    gr.tmp <- m.list %>% filter(Dist <= distance) %>% graph_from_data_frame(.,directed = FALSE)
    gr.tmp <- simplify(gr.tmp, remove.multiple = TRUE, remove.loops = TRUE) %>% as.undirected()

    if(fast ==TRUE) {
      cluster <- components(gr.tmp)
    }else{
        cluster <- cluster_louvain(gr.tmp)
    }

    Nc <- as.numeric(max(cluster$membership))
    cent <- centralization.degree(gr.tmp)
    results = data.frame(Source = as.character(vertex.attributes(gr.tmp)$name),
                         centrality = cent$res,
                         cluster = cluster$membership, stringsAsFactors = F)
    class(results) <- append(class(results),"nr_list")
    return(results)

  }else{
    print("Error: number, fraction or distance must be provided")
  }
}





