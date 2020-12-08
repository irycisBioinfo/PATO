#' Title
#'
#' @param ref
#' @param query
#'
#' @return
#' @export
#'
#' @examples
overlap_mafs <- function(ref,query)
{
  result = data.frame()
  for(i in 1:nrow(query))
  {
    for(j in 1:nrow(ref))
    {
      if(query$src[i] == ref$src[j])
      {
        if((query$start[i] <= ref$start[j]) && (query$end[i] >= ref$start[j]) && (query$end[i] <= ref$end[j]))
        {
          result <- bind_rows(result,data.frame(src = query$src[i], start = ref$start[j], end = query$end[i]))

        }else if ((query$start[i] >= ref$start[j]) && (query$start[i] <= ref$end[j]) && (query$end[i] >= ref$end[j]))
        {
          result <- bind_rows(result,data.frame(src = query$src[i], start = query$start[i], end = ref$end[j]))

        }else if ((query$start[i] >= ref$start[j]) && (query$end[i] <= ref$end[j]))
        {
          result <- bind_rows(result,data.frame(src = query$src[i], start = query$start[i], end = query$end[i]))

        }else if ((query$start[i] <= ref$start[j]) && (query$end[i] >= ref$query[j]))
        {
          result <- bind_rows(result,data.frame(src = query$src[i], start = ref$start[j], end = ref$end[j]))
        }
      }
    }
  }
  return(result)
}
