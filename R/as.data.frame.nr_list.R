#' Wrapper to change format between nr_list and data.frame
#'
#' @param data A `nr_list` object
#'
#' @return A `data.frame`
#' @export
#'
#' @examples
as.data.frame.nr_list <- function(data)
{
  return(data.frame(
      Source = as.character(data$Source),
    centrality = data$centrality,
    cluster = data$cluster))
  UseMethod("as.data.frame")
  }


