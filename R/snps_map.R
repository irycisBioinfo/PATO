#' SNPs Map
#'
#' snps_map is a function that makes a similar plot to manhattan plots.
#' This function plot the cumulative number of snps per each contig-position pair.
#'
#' @param core_snp_genome A \emph{core_snp_genome} object.
#'
#' @return A \emph{ggplot2} object
#' @export
#'
#' @import microseq
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#'
#' @examples
snps_map <- function(core_snp_genome)
{

  mhtt <- microseq::readFasta(core_snp_genome$reference) %>%
    mutate(cont_len = str_length(Sequence)) %>%
    arrange(desc(cont_len)) %>%
    mutate(tot=cumsum(cont_len)-cont_len) %>%
    rename(CHROM = Header) %>% select(CHROM,tot,cont_len) %>%
    inner_join(core_snp_genome$vcf,.,by=c("CHROM"="CHROM")) %>%
    mutate(POScum=POS+tot) %>%
    group_by(CHROM,POScum,cont_len) %>%
    summarise(Nsnps = n())

  axisdf <- mhtt %>% group_by(CHROM) %>% summarize(center=( max(POScum) + min(POScum) ) / 2 )


  plot <- mhtt %>% ggplot(aes(x = POScum, y = Nsnps, color = CHROM)) +
    geom_point(alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("grey", "skyblue"), 50))+
    scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center )+
    labs(x ="", y= "Number of SNPs") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1), legend.position = "none")

  print(plot)

  return(plot)


}
