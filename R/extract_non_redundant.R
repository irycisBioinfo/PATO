#' Create a accnet or mash object with the representatives members of each cluster
#'
#' This function creates a non-redundant representation of the input object
#' (\emph{accnet} or \emph{mash}).
#'
#' @param data \emph{accnet} or \emph{mash} object
#' @param nr_list \emph{data.frame} with two columns (Source, Target) or
#' \emph{nr_list} object created with \emph{non_redundant()} function. For
#' \emph{data.frame}, each first member of the cluster will be selected as
#' representative. In the caso of \emph{nr_list} object. The member with higher
#' centrality will be selected as representative.
#'
#' @return \emph{accnet} or \emph{mash} object
#' @export
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @import data.table
#'
extract_non_redundant <- function(data, nr_list)
{
  if (is(nr_list,"nr_list"))
  {
    members = as.data.frame(nr_list)
    members <- members %>% group_by(cluster) %>% top_n(1, centrality)

    if (is(data,"accnet"))
    {
      list <-  data$list %>% as_tibble() %>% semi_join(members)
      matrix <-  list %>%
        mutate(value = 1) %>%
        select(Source, Target, value) %>%
        as.data.table() %>%
        spread(Target, value, fill = 0)
      tmp <-  list %>%
        select(Target) %>%
        distinct() %>%
        rename(ID = Target)
      annot <-
        data$annot %>% as_tibble() %>%  semi_join(tmp, by = c("ID"))
      results  <-  list(list = list,
                        matrix = matrix,
                        annot = annot)
      class(results) <- "accnet"
      return(results)


    } else if (is(data,"mash"))
    {
      table <-
        data$table %>% semi_join(members) %>% semi_join(members %>% rename(Target = Source))
      matrix <- table %>% spread(Target, Dist) %>% column_to_rownames("Source")
      results <- list(table = table, matrix = matrix)
      class(results) = "mash"
      return(results)

    }
  } else{
    if (grep("Source", colnames(nr_list)) &
        grep("cluster", colnames(nr_list)))
    {
      nr_list <-
        nr_list %>% group_by(cluster) %>% summarise(Source = first(Source))
      if (is(data,"accnet"))
      {
        list <-  data$list %>% as_tibble() %>% semi_join(nr_list)
        matrix <-  list %>%
          mutate(value = 1) %>%
          select(Source, Target, value) %>%
          as.data.table() %>%
          spread(Target, value, fill = 0)
        tmp <-  list %>%
          select(Target) %>%
          distinct() %>%
          rename(ID = Target)
        annot <-
          data$annot %>% as_tibble() %>%  semi_join(tmp, by = c("ID"))
        results  <-
          list(list = list,
               matrix = matrix,
               annot = annot)
        class(results) <- "accnet"
        return(results)


      } else if (is(data,"mash"))
      {
        table <-
          data$table %>%
          semi_join(nr_list) %>%
          semi_join(nr_list %>%
                      rename(Target = Source))

        matrix <- table %>% spread(Target, Dist) %>% column_to_rownames("Source")
        results <- list(table = table, matrix = matrix)
        class(results) = "mash"
        return(results)

      }
    }
  }
}
