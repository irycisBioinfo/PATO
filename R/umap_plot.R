#' Plot UMAP representation of a data set.
#'
#' UMAP (Uniform Manifold Approximation and Projection for Dimension Reduction) is a dimension
#' reduction technique that can be used for visualisation similarly to t-SNE, but also for
#' general non-linear dimension reduction. This function uses UMAP to represent the population structure
#' of the data set (\emph{accnet} or \emph{mash})
#'
#' @param data An \emph{accnet} or \emph{mash} object
#' @param plot print plot?
#' @param cluster Two column data.frame with the nodes (column 1) and the cluster
#' or classification(column 2)
#' @param ... Passed to uwot::umap function.
#'
#' @return A matrix of optimized coordinates.
#' @export
#'
#' @seealso \code{\link{uwot}}
#' @references McInnes, L, Healy, J, UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction, ArXiv e-prints 1802.03426, 2018
#' @references James Melville (2019). uwot: The Uniform Manifold Approximation and Projection (UMAP) Method for
#' Dimensionality Reduction. R package version 0.1.4. https://CRAN.R-project.org/package=uwot
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import ggplot2
#' @importFrom data.table fread
#' @importFrom randomcoloR distinctColorPalette
umap_plot <- function(data, plot = TRUE, cluster,...)
{
  if(is(data,"accnet"))
  {
    umap <- umap(data$matrix %>% column_to_rownames("Source"),...)
    umap <- as.data.frame(umap)
    colnames(umap) <- c("X","Y")
    umap$Source <- data$matrix$Source

  }else if(is(data,"mash"))
  {
    umap <- umap(data$matrix %>% as.matrix(),...)
    umap <- as.data.frame(umap)
    colnames(umap) <- c("X","Y")
    umap$Source <- rownames(data$matrix)
  }else{
    stop("data must be a mash or accnet object")
  }

  if(!missing(cluster))
  {
    colnames(cluster) <- c("Source","Cluster")
    cluster <- cluster %>% mutate(Cluster = as.factor(Cluster))
    umap <- umap %>% left_join(cluster)
  }
  if(plot)
  {
    if(missing(cluster))
    {
      print(umap %>% ggplot(aes(x = X, y = Y))+
              geom_point() +
              theme_bw() +
              xlab("umap_x") +
              ylab("umap_y"))
    }else{
      print(umap %>% ggplot(aes(x = X, y = Y, color = Cluster))+
              geom_point() +
              theme_bw() +
              xlab("umap_x") +
              ylab("umap_y"))

    }
  }

  return(umap)
}
