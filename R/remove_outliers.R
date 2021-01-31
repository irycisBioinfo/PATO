#' Remove the outliers of a data set (mash or accnet)
#'
#' Remove the outliers (or any desire element) of an \emph{accnet} or \emph{mash} object.

#' @param data An \emph{accnet} or \emph{mash} object
#' @param outliers A one data frame with a column "Source".
#' It must contains the elements to remove
#'
#' @return An \emph{accnet} or \emph{mash} object without the removed elements.
#' @export
#'
#' @seealso \code{\link{outliers}}
#' 
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @importFrom data.table fread
#'
remove_outliers <- function(data,outliers)
{
  if(is(data,"accnet"))
  {
    list <- data$list %>% as_tibble() %>% anti_join(outliers)
    matrix <- list %>%
      mutate(value =1) %>%
      select(Source,Target,value) %>%
      spread(Target,value, fill = 0)
    annot <- data$annot %>% as_tibble() %>% anti_join(list %>%
                                        select(Target) %>%
                                        distinct() %>% rename(ID = Target))
    results <-  list(list = list, matrix = matrix, annot = annot, path = data$path)
    class(results) <- append(class(results),"accnet")
    return(results)
  }else if (is(data,"mash"))
  {
    list <- data$table %>% as_tibble() %>% anti_join(outliers)
    list <- list %>% anti_join(outliers %>% rename(Target = Source))
    matrix <- list %>%
      select(Source,Target,Dist) %>%
      spread(Target,Dist, fill = 0) %>% column_to_rownames("Source")
    results <- list(table =list, matrix = matrix,path = data$path)
    class(results) <- append(class(results),"mash")
    return(results)

  } else{
    print("Object accnet or mash must by provided")
  }

}
