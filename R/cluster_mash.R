#' Internal function to cluster mash data
#'
#' @param data \emph{mash} object
#' @param method Clustering method
#' @param n_cluster number of cluster (if applicable)
#' @param d_reduction boolean. Dimensional reduction using umap
#' @param ... additional parameters
#'
#' @return \emph{data.frame} with two columns: \emph{Source} and \emph{Target}
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import uwot
#' @import mclust
#' @importFrom data.table fread
#' @import dbscan
#'
cluster_mash <- function(data, method, n_cluster,d_reduction,...)
{
  Dist <-  as.dist(data$matrix)
  if (method == "mclust")
  {
    if (!d_reduction)
    {
      min <- 2
      max <- 9
      max_clust <- 9

      for (i in 1:10)
      {
        cluster <- Mclust(Dist, min:max)
        max_clust <-  max(cluster$classification)
        if (max_clust == max)
        {
          min = max - 1
          max = max + 10
        } else{
          break
        }

      }

      cluster <-
        cluster$classification %>% as.data.frame() %>% tibble::rownames_to_column("Source")
      colnames(cluster) <- c("Source", "Cluster")
    } else{
      min <- 2
      max <- 9
      max_clust <- 9


      umap <-  uwot::umap(as.matrix(Dist))
      rownames(umap) <- rownames(as.matrix(Dist))
      for (i in 1:10)
      {
        cluster <- Mclust(umap, min:max)
        max_clust <-  max(cluster$classification)
        if (max_clust == max)
        {
          min = max - 1
          max = max + 10
        } else{
          break
        }

      }

      cluster <-
        cluster$classification %>% as.data.frame() %>% rownames_to_column("Source")
      colnames(cluster) <- c("Source", "Cluster")
    }

  } else if (method == "upgma")
  {
    tree <-  hclust(Dist, method = "average")
    cluster <- cutree(tree, n_cluster)
    cluster <-
      cluster %>% as.data.frame() %>% tibble::rownames_to_column("Source")
    colnames(cluster) <- c("Source", "Cluster")

  } else if (method == "ward.D2")
  {
    tree <-  hclust(Dist, method = "ward.D2")
    cluster <- cutree(tree, n_cluster)
    cluster <-
      cluster %>% as.data.frame() %>% tibble::rownames_to_column("Source")
    colnames(cluster) <- c("Source", "Cluster")


  } else if (method == "hdbscan")
  {
    umap <-  uwot::umap(as.matrix(Dist))
    rownames(umap) <- rownames(as.matrix(Dist))
    cluster <- dbscan::hdbscan(umap, minPts = 0.05 * nrow(umap))
    cluster <-
      cluster$cluster %>%  as.data.frame() %>% mutate(Source = rownames(umap))
    colnames(cluster) <- c("Cluster", "Source")
    return(cluster %>% select(Source, Cluster))

  } else{
    stop("Unrecognized method")
  }
  return(cluster)
}
