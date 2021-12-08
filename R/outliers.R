#' Find outliers in the dataset.
#'
#' Find outliers in the dataset using a threshold distance (mash o jaccard distance) or by
#' standard deviation times. This function also plot a manhattan plot to visualize the
#' structure of the distance distribution of the input object.
#'
#' @param data An \emph{accnet} or \emph{mash} object
#' @param threshold threshold for outliers (outliers >= threshold)
#' @param n_sd Number of standard deviation as threshold.
#' @param plot \emph{Boolean} plot manhattan?
#'
#' @return a data.frame with the set of outliers.
#' @export
#'
#' @seealso \code{\link{remove_outliers}}
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @import ggplot2
#' @import magrittr
#' @import data.table
#' @import parallelDist
#'
#'
outliers <- function(data, threshold, n_sd, plot = TRUE)
{
  if (missing(threshold) & missing(n_sd))
  {
    stop("threshold or n_sd attributes must be provided")
  }

  if (is(data,"mash"))
  {
    m.list <- data$matrix %>% as.data.frame() %>%
      rownames_to_column("Source") %>%
      as.data.table() %>%
      gather(Target, Dist,-Source) %>%
      group_by(Source) %>%
      summarise(mean = mean(Dist))
  } else if (is(data,"accnet"))
  {
    if(is.null(data$dist))
    {
      m.list <- data$matrix %>%
        column_to_rownames("Source") %>%
        as.matrix() %>%
        parallelDist(., method = "binary") %>%
        as.matrix() %>%
        as.data.frame() %>%
        rownames_to_column("Source") %>%
        as.data.table() %>%
        gather(Target, Dist,-Source) %>%
        group_by(Source) %>%
        summarise(mean = mean(Dist))
    }else{
      m.list <- data$dist %>%
        as.data.frame() %>%
        rownames_to_column("Source") %>%
        as.data.table() %>%
        gather(Target, Dist,-Source) %>%
        group_by(Source) %>%
        summarise(mean = mean(Dist))
    }
  } else{
    stop("Error: object mash or accnet must be provided")

  }


  if (!missing(threshold))
  {
    if (plot)
    {
      print(
        m.list %>% mutate(color = ifelse(mean >= threshold, "red","black" )) %>%
          ggplot(aes(x = Source, y = mean, color = color)) + scale_color_identity()+
          geom_point(show.legend = FALSE) +
          geom_hline(yintercept = threshold, color = "coral") +
          geom_text(
            aes(
              y = threshold,
              x = 1,
              label = "Threshold",
              color = "coral"
            ),
            vjust = "inward",
            hjust = "inward",
            show.legend = FALSE,
            inherit.aes = FALSE
          ) +
          theme(axis.ticks.x = element_blank(),
                axis.text.x = element_blank())
      )
    }
    return(m.list %>% filter(mean >= threshold) %>% select(Source))

  } else if (!missing(n_sd))
  {
    meanAbs <- mean(m.list$mean)
    sdAbs <- sd(m.list$mean)
    if (plot)
    {
      print(
        m.list %>% mutate(color = ifelse(mean >= meanAbs, "red","black" )) %>%
          ggplot(aes(x = Source, y = mean, color = color)) +
          geom_point(show.legend = FALSE) +
          geom_hline(yintercept = meanAbs,
                     color = "red") +
          geom_text(
            aes(
              y = meanAbs,
              x = 1,
              label = "Mean",
              color = "coral"
            ),
            vjust = "inward",
            hjust = "inward",
            show.legend = FALSE,
            inherit.aes = FALSE
          ) +
          geom_hline(yintercept = meanAbs + (n_sd * sdAbs),
                     color = "coral") +
          geom_text(
            aes(
              y = meanAbs + (n_sd * sdAbs),
              x = 1,
              label = paste(n_sd, " x SD", sep = "", collapse = ""),
              color = "coral"
            ),
            vjust = "inward",
            hjust = "inward",
            show.legend = FALSE,
            inherit.aes = FALSE
          ) +
          theme(axis.ticks.x = element_blank(),
                axis.text.x = element_blank())
      )

    }
    return(m.list %>% filter(mean >= (meanAbs + (n_sd * sdAbs))) %>% select(Source))
  }
}
