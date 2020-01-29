#' Plot K-Nearest Neighbour Network
#'
#' This function uses three.js javascript library (https://threejs.org/)
#' (R package -> threejs: https://bwlewis.github.io/rthreejs) to draw
#' the k-nnn. You can select the layout algorithm to arrange the network.
#' Moreover, the network can visualize any kind of clustering (internal or external)
#' an color the nodes by cluster.
#' Finally the funtion return a list with the network, the layout and the colors
#' (to use with other functions)
#'
#' @param net An \emph{igraph} object (preferably knnn result)
#' @param layout String with the layout option c("fr","kk","DrL","mds")
#' \itemize{
#' \item{fr :  Fruchterman-Reingold layout. Place vertices on the plane
#' using the force-directed layout algorithm by Fruchterman and Reingold. (default) }
#' \item{kk: The Kamada-Kawai layout algorithm. Place the vertices on the plane,
#' or in the 3d space, based on a phyisical model of springs. }
#' \item{DrL: The DrL graph layout generator. DrL is a force-directed graph
#' layout toolbox focused on real-world large-scale graphs, developed by Shawn
#' Martin and colleagues at Sandia National Laboratories. }
#' \item{mds: Graph layout by multidimensional scaling. Multidimensional scaling
#' of some distance matrix defined on the vertices of a graph. }}
#' @param dim Number of network dimensions (2 or 3). Default = 3.
#' @param cluster Two column data.frame with the nodes (column 1) and the cluster
#' or classification(column 2)
#' @param ... Passed to layout_with_
#'
#' @return A list with the elements "graph", "layout" and "colors"
#' @export
#'
#' @seealso \code{\link{igraph}}
#'
#' @examples
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @import threejs
#' @import data.table
#' @import manipulateWidget
#' @importFrom randomcoloR distinctColorPalette
plot_knnn_network <- function(net, layout ="fr", dim = 3, cluster, ...)
{
  if(!is(net,"igraph"))
  {
	  stop("net must be an igraph object")
  }
  if (layout == "fr")
  {
    coords <- igraph::layout_with_fr(net, dim = dim)

  } else if (layout == "kk")
  {
    coords <- igraph::layout_with_kk(net, dim = dim)
  } else if (layout == "DrL")
  {
    coords <- igraph::layout_with_drl(net, dim = dim)
  } else if (layout == "mds")
  {
    coords <- igraph::layout_with_mds(net, dim = dim)
  }
  else{
    stop("Unrecognised layout algorithm")
  }
  if (dim == 2)
  {
    coords_plot <- cbind(coords, rep(5, nrow(coords)))
  } else{
    coords_plot <- coords
  }

  if (missing(cluster))
  {
    print(
      threejs::graphjs(
        net,
        coords_plot,
        vertex.size = 0.2,
        edge.color = "lightGrey",
        vertex.label = vertex.attributes(net)$name,
        ...,
      )
    )
    return(list(graph = net, layout = coords))
  } else{

    colnames(cluster) <- c("Source", "Cluster")
    cluster = cluster %>% full_join(data.frame(
      Cluster = unique(cluster$Cluster),
      color = distinctColorPalette(length(unique(cluster$Cluster)))
    ))

    print(
      combineWidgets(
        threejs::graphjs(
            net,
            coords_plot,
            vertex.size = 0.2,
            vertex.color = cluster$color,
            edge.color = "lightGrey",
            vertex.label = vertex.attributes(net)$name, ...,
          ),
        tags$div(
          tags$b("CLUSTERS", style= "color:black;font-size:9px"),
          tags$div("1", style = "color:red;font-size:15px"),
          tags$div("2", style ="color:blue;font-size:15px"))
      )
    )

    return(list(
      graph = net,
      layout = coords,
      colors = cluster
    ))
  }

}












