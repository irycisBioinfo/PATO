## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  install.packages("devtools")

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  devtools::install_github("https://github.com/irycisBioinfo/PATO.git", build_vignettes = TRUE)

## ----setup, echo=FALSE--------------------------------------------------------
knitr::opts_chunk$set(
	message = FALSE,
	warning = FALSE,
	include = FALSE
)

## ----eval=FALSE---------------------------------------------------------------
#  my_mash <- mash(my_list, type ="prot")

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  my_mash <- mash(my_list, n_cores = 20, sketch = 1000, kmer = 21, type = "prot")

## ----eval = FALSE, include=TRUE-----------------------------------------------
#  my_mmseq <- mmseqs(my_list)

## ----eval=FALSE---------------------------------------------------------------
#  my_mmseq <- mmseqs(my_list, coverage = 0.8, identity = 0.8, evalue = 1e-6, n_cores = 20, cov_mode = 0, cluster_mode = 0)
#  

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  my_accnet <- accnet(my_mmseq, threshold = 0.8, singles = TRUE)

## ----eval=FALSE, include= TRUE------------------------------------------------
#  my_accnet <- mmseqs(my_files) %>% accnet(threshold = 0.8)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  library(PATO)
#  # We strongly recommend to load tidyverse meta-package
#  library(tidyverse)
#  setwd("/myFolder/")

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  files <- read.table("my_list.txt",header = FALSE)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  efaecium_mash_all <- mash(files, n_cores = 20, type = 'prot')
#  efaecium_mm<- mmseqs(files, coverage = 0.8, identity = 0.8, evalue = 1e-6, n_cores = 20)

## ----eval=FALSE, include=FALSE------------------------------------------------
#  efaecium_accnet_all <- accnet(efaecium_mm,threshold = 0.8, singles = FALSE)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  outl <- outliers(efaecium_mash_all,threshold = 0.06)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  efaecium_mash_all <-remove_outliers(efaecium_mash_all, outl)
#  efaecium_accnet_all <-remove_outliers(efaecium_accnet_all, outl)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  files <-  anti_join(files, outl, by=c("V1"="Source"))

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  # For select 800 non redundant samples
#  nr_list <- non_redundant(efaecium_mash_all,number = 800)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  # For select a subset of samples with 99.99% of indentity
#  nr_list <- non_redundant(efaecium_mash_all, distance = 0.0001)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  efaecium_accnet <- extract_non_redundant(efaecium_accnet_all, nr_list)
#  efaecium_mash <- extract_non_redundant(efaecium_mash_all, nr_list)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  nr_list <- non_redundant_hier(efaecium_mash_all,800, partitions = 10)
#  efaecium_accnet <- extract_non_redundant(efaecium_accnet_all, nr_list)
#  efaecium_mash <- extract_non_redundant(efaecium_mash_all, nr_list)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  cp <- core_plots(efaecium_mm,reps = 10, threshold = 0.95, steps = 10)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  core <- mmseqs(files) %>% core_genome()
#  export_core_to_fasta(core,"core.aln")
#  

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  core_tree = read.tree("core.aln.treefile")
#  core_tree %>% midpoint.root() %>% ggtree()

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  efaecium_accnet_all <- accnet(efaecium_mm,threshold = 0.8, singles = FALSE)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  ef_nr<- non_redundant(efaecium_mash,number = 800 )
#  efaecium_accnet_nr <- extract_non_redundant(efaecium_accnet, ef_nr)
#  ef_800_cl <- clustering(efaecium_accnet_nr, method = "mclust", d_reduction = TRUE)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  export_to_gephi(efaecium_accnet_nr, "accnet800", cluster = ef_800_cl)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  accnet_enr_result <- accnet_enrichment_analysis(efaecium_accnet_nr, cluster = ef_800_cl)
#  
#  accnet_enr_result
#  
#  # A tibble: 854,887 x 14
#     Target Source Cluster perClusterFreq ClusterFreq ClusterGenomeSi… perTotalFreq TotalFreq OdsRatio  pvalue
#     <chr>  <chr>    <dbl>          <dbl>       <int>            <int>        <dbl>     <int>    <dbl>   <dbl>
#   1 WP_00… GCF_0…       1          0.248          53              214        0.145       172     1.71 2.44e-6
#   2 WP_00… GCF_0…       1          0.248          53              214        0.145       172     1.71 2.44e-6
#   3 WP_00… GCF_0…       1          0.248          53              214        0.145       172     1.71 2.44e-6
#   4 WP_00… GCF_0…       1          0.248          53              214        0.145       172     1.71 2.44e-6
#   5 WP_00… GCF_0…       1          0.248          53              214        0.145       172     1.71 2.44e-6
#   6 WP_00… GCF_0…       1          0.248          53              214        0.145       172     1.71 2.44e-6
#   7 WP_00… GCF_0…       1          0.248          53              214        0.145       172     1.71 2.44e-6
#   8 WP_00… GCF_0…       1          0.248          53              214        0.145       172     1.71 2.44e-6
#   9 WP_00… GCF_0…       1          0.248          53              214        0.145       172     1.71 2.44e-6
#  10 WP_00… GCF_0…       1          0.248          53              214        0.145       172     1.71 2.44e-6
#  # … with 854,877 more rows, and 4 more variables: padj <dbl>, AccnetGenomeSize <int>, AccnetProteinSize <int>,
#  #   Annot <chr>

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  accnet_with_padj(accnet_enr_result) %>% export_to_gephi("accnet800.padj", cluster = ef_800_cl)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  mash_tree_fastme <- similarity_tree(efaecium_mash)
#  mash_tree_NJ <- similarity_tree(efaecium_mash, method = "NJ")
#  mash_tree_upgma <- similarity_tree(efaecium_mash,method = "UPGMA")
#  accnet_tree_upgma <- similarity_tree(efaecium_accnet,method = "UPGMA")

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  ggtree(mash_tree_NJ) + geom_tippoint()

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  library(dendextend)
#  tanglegram(ladderize(mash_tree_upgma), ladderize(accnet_tree_upgma), fast = TRUE)
#  

## ----eval=FALSE, include=FALSE------------------------------------------------
#  mmseqs(files) %>% accnet() %>% export_accnet_aln(.,file ="my_accnet_aln.fasta")

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  acc_tree = read.tree("acc.aln.treefile")
#  acc_tree %>% midpoint.root() %>% ggtree()

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  
#  # K-NNN with 10 neighbours and with repetitions
#  knnn_mash_10_w_r <- knnn(efaecium_mash,n=10, repeats = TRUE)
#  # K-NNN with 25 neighbours and with repetitions
#  knnn_mash_25_w_r <- knnn(efaecium_mash,n=25, repeats = TRUE)
#  # K-NNN with 50 neighbours and with repetitions
#  knnn_mash_50_w_r <- knnn(efaecium_mash,n=50, repeats = TRUE)
#  

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  export_to_gephi(knnn_mash_50_w_r,file = "knnn_50_w_r.tsv")

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  plot_knnn_network(knnn_mash_50_w_r)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  ef_cl_mclust_umap <- clustering(efaecium_mash, method = "mclust",d_reduction = TRUE)

## ----eval =FALSE, include=TRUE------------------------------------------------
#  ef_cl_knnn <-clustering(knnn_mash_50_w_r, method = "louvain")

## ----eval =FALSE, include=TRUE------------------------------------------------
#  umap_mash <- umap_plot(efaecium_mash)

## ----eval =FALSE, include=TRUE------------------------------------------------
#  umap_mash <- umap_plot(efaecium_mash, cluster = ef_cl_mclust_umap)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  cl_louvain = clustering(knnn_mash_25_wo_r, method = "louvain")
#  plot_knnn_network(knnn_mash_25_wo_r, cluster = cl_louvain, edge.alpha = 0.1)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  library(tidyverse)
#  
#  annotation <- annotate(efaecium_mm, type = "prot",database = c("AbR","VF_A"))
#  
#  annotation %>%
#    filter(pident > 0.95 ) %>% #remove all hits with identity lower than 95%
#    filter(evalue < 1e-6) %>%  #remove all hits with E-Value greater than 1e-6
#    group_by(Genome,Protein) %>% top_n(1,target) #select only the best hit for each protein-genome.
#  

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  heatmap_of_annotation(annotation %>% filter(DataBase =="AbR"), #We select only "AbR" results
#                        min_identity = 0.99)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  network_of_annotation(annotation %>% filter(DataBase =="AbR"), min_identity = 0.99) %>% export_to_gephi("annotation_Network")

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  res <- pangenomes_from_files(files,distance = 0.03,min_pange_size = 10,min_prot_freq = 2)
#  
#  export_to_gephi(res,"/storage/tryPATO/firmicutes/pangenomes_gephi")

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  assembly = data.table::fread("ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt",sep = "\t",skip = 1, quote = "")
#  
#  colnames(assembly) = gsub("# ","",colnames(assembly))
#  
#  annot = res$members %>%
#    mutate(file = basename(path)) %>%
#    separate(file,c("GCF","acc","acc2"),sep = "_", remove = FALSE) %>%
#    unite(assembly_accession,GCF,acc,sep = "_") %>%
#    left_join(assembly) %>%
#    separate(organism_name,c("genus","specie"), sep = " ") %>%
#    group_by(pangenome,genus,specie) %>%
#    summarise(N = n()) %>%
#    distinct() %>%
#    group_by(pangenome) %>%
#    top_n(1,N) %>%
#    mutate(ID = paste("pangenome_",pangenome,"_rep_seq.fasta", sep = "",collapse = "")) %>%
#    write_delim("/storage/tryPATO/firmicutes/pangenomes_gephi_extra_annot.tsv",delim = "\t", col_names = TRUE)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  sh <-  sharedness(res)
#  sh  <-  sh %>% as.data.frame() %>%
#    rownames_to_column("ID") %>%
#    inner_join(annot) %>%
#    unite(Name,genus,specie,pangenome, sep = "_") %>%
#    select(-ID,-N)%>%
#    column_to_rownames("Name")
#  
#  colnames(sh) = rownames(sh)
#  
#  pheatmap::pheatmap(sh,clustering_method = "ward.D2", clustering_distance_cols = "correlation", clustering_distance_rows = "correlation")

