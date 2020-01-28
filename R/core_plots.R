#' Plot the size of core-genome, accessory-genome and pan-genome
#'
#' Takes the a \emph{mmseq} object and plot the size of core, accessory and
#' pan-genome data in the specified points and with the specified number of
#' replicates.
#'
#' @param data A \emph{mmseq} object.
#' @param steps number or vector of points to calculate sizes.
#' @param reps number of replicates for each point.
#' @param threshold Minimun frequency (0-1) of the core-genome genes/proteins.
#' @param type Output in format \emph{pato} or \emph{roary} style.
#'
#'
#' @return A list with a \emph{ggplot2} object and a \emph{data.frame} with the values of the final step.
#'
#' @note It's supossed the a gene/protein of the core-genome must be present in all
#' the genomes of the dataset. Nevertheless, by random selection and/or errors
#' in sequencing proccess (sequencing, assembly, ORF finding etc..) some
#' genes/proteins could be missing. Pato define the genome as, pangeno (the sum of all genomes),
#' core-genome.
#'
#' @export
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @import ggplot2
#' @import data.table
#'
core_plots <- function(data, steps = 10, reps = 10, threshold = 0.98, type = "pato")
{
  if (!is(data,"mmseq"))
  {
    stop("datamust be a 'mmseq' object")
  }



  results = data.frame(Data="Core",N=0,Iter=0,Value=0)
  results = results[-1,]
  results$Data = as.character(results$Data)

  data_plots <-  data$table %>% select(Prot_prot, Genome_genome) %>% distinct() %>% as_tibble()
  genomes <- data$table %>% distinct(Genome_genome) %>% as_tibble()
  nGenomes  <-  data_plots %>% distinct(Genome_genome) %>% count()
  if(length(steps) <=1)
  {
    intervals <- seq(10,nGenomes$n, by = (nGenomes$n - 10) / steps)
  }else {
    intervals <- steps
  }


  if(type == "pato")
  {

    for (i in intervals)
    {
        for (j in 1:reps) {
          tmp <-  data_plots %>%
            inner_join(genomes %>% sample_n(i), by = "Genome_genome") %>%
            group_by(Prot_prot) %>%
            summarise(freq = n())
          core <-  tmp %>%
            filter(freq >= i * threshold) %>% count()

          pange <-  tmp %>% filter(freq > 1) %>% count()

          acc <-  pange$n - core$n


          results <-  bind_rows(results,
                              data.frame(
                                Data = "Core",
                                N = i,
                                Iter = j,
                                Value = core$n
                              ))
          results  <-  bind_rows(results,
                              data.frame(
                                Data = "Pangenome",
                                N = i,
                                Iter = j,
                                Value = pange$n
                              ))
          results <-  bind_rows(results,
                              data.frame(
                                Data = "Accessory",
                                N = i,
                                Iter = j,
                                Value = acc
                              ))
        }
      }

  }else if (type =="roary")
  {
    for (i in intervals)
    {
      for (j in 1:reps) {
        tmp <-  data_plots %>%
          inner_join(genomes %>% sample_n(i), by = "Genome_genome") %>%
          group_by(Prot_prot) %>%
          summarise(freq = n())


        core <-  tmp %>% filter(freq >= i * 0.99) %>% count()
        softcore <-  tmp %>% filter(freq >= i * 0.95) %>% filter(freq < i * 0.99) %>% count()
        shellgenes <- tmp %>% filter(freq >= i * 0.15) %>% filter(freq < i * 0.95) %>% count()
        cloudgenes <- tmp %>% filter(freq < i * 0.15) %>% count()




        results <-  bind_rows(results,
                              data.frame(
                                Data = "Core",
                                N = i,
                                Iter = j,
                                Value = core$n
                              ))
        results  <-  bind_rows(results,
                               data.frame(
                                 Data = "SoftCore",
                                 N = i,
                                 Iter = j,
                                 Value = softcore$n
                               ))
        results <-  bind_rows(results,
                              data.frame(
                                Data = "ShellGenes",
                                N = i,
                                Iter = j,
                                Value = shellgenes$n
                              ))
        results <-  bind_rows(results,
                              data.frame(
                                Data = "CloudGenes",
                                N = i,
                                Iter = j,
                                Value = cloudgenes$n
                              ))

      }
    }
  }

  data = results %>% group_by(N, Data) %>% summarise(n_mean = mean(Value), sd = sd(Value))
  p <- data %>% ggplot(aes(x = N,y = n_mean,group = Data,color = Data)) +
    geom_line() +
    geom_point() + geom_errorbar(aes(ymin = n_mean-sd, ymax = n_mean+sd)) +
    theme_light() +
    xlab("Number of Genomes") +
    ylab("Number of Proteins")

  print(p)
  return(list(plot = p, data = data))
}
