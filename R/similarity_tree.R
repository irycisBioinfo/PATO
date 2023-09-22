#' Similarity tree
#'
#' Build a similarity tree (pseudo-phylogenetic) of the samples,
#' using whole-genome (\emph{mash}) or accessory genome (\emph{accnet}).
#'
#' @param data An \emph{accnet} or \emph{mash} object.
#' @param method Tree method c("fastme", "NJ", "UPGMA"):
#' \itemize{
#' \item{fastme: Tree Estimation Based on the Minimum Evolution Algorithm}
#' \item{NJ: Neighbor-Joining Tree Estimation}
#' \item{UPGMA: Unweighted pair group method with arithmetic mean}
#' \item{ward.D2: Ward's minimum variance metho}
#' \item{median: WPGMC Weighted Pair Group Method with Arithmetic Mean}
#' \item{complete: complete linkage  }
#' }
#'
#' @seealso \code{\link[ape]{fastme}}
#' @seealso \code{\link[ape]{nj}}
#' @seealso \code{\link[stats]{hclust}}
#'
#' @return An object of class \emph{phylo}
#' @export
#'
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom data.table fread
#' @import ape
#' @import parallelDist


similarity_tree <- function(data, method ="NJ")
{

  if(is(data,"accnet"))
  {
    if(is.null(data$dist))
    {
      dist <- data$matrix %>%
        column_to_rownames("Source") %>%
        as.matrix()%>%
        parallelDist(., method = "binary")
    }else{
      dist <- data$dist %>% as.dist()
    }


  }else if (is(data,"mash"))
  {
    dist = as.dist(data$matrix)

  }else{
    stop("data must be a mash or accnet object")
  }

  if(method == "fastme")
  {
    tree = ape::fastme.bal(dist)
    return(phangorn::midpoint(tree))
  }else if (method == "NJ")
  {
    tree = ape::bionj(dist)
    return(phangorn::midpoint(tree))
  }else if (method =="UPGMA")
  {
    tree = hclust(dist, method ="average")
    return(ape::as.phylo(tree))
  }else if (method =="ward.D2")
  {
    tree = hclust(dist, method ="ward.D2")
    return(ape::as.phylo(tree))
  }
  else if (method =="median")
  {
    tree = hclust(dist, method ="median")
    return(ape::as.phylo(tree))

  }else if (method =="complete")
  {
    tree = hclust(dist, method ="complete")
    return(ape::as.phylo(tree))

  }else{
    stop("invalid method")
  }

}
