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
#  library(pato)
#  # We strongly recommend to load tidyverse meta-package
#  library(tidyverse)
#  setwd("/myFolder/")

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  files <- read.table("my_list.txt",header = FALSE)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  species <- classifier(files,n_cores = 20,type = 'prot')
#  
#  species %>% group_by(organism_name) %>% summarise(Number = n())
#  
#  # A tibble: 8 x 2
#    organism_name                             Number
#    <chr>                                      <int>
#  1 Citrobacter gillenii                           1
#  2 Escherichia coli O157:H7 str. Sakai          323
#  3 Escherichia coli str. K-12 substr. MG1655   2520
#  4 Escherichia marmotae                           4
#  5 Shigella boydii                               18
#  6 Shigella dysenteriae                           7
#  7 Shigella flexneri 2a str. 301                 45
#  8 Shigella sonnei                               40

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  files = species %>% filter(!grepl("Citrobacter",organism_name)) %>% filter(!grepl("marmotae",organism_name))

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  ecoli_mash_all <- mash(files, n_cores = 20, type = 'prot')
#  ecoli_mm<- mmseqs(files, coverage = 0.8, identity = 0.8, evalue = 1e-6, n_cores = 20)

## ----eval=FALSE, include=FALSE------------------------------------------------
#  ecoli_accnet_all <- accnet(ecoli_mm,threshold = 0.8, singles = FALSE)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  outl <- outliers(ecoli_mash_all,threshold = 0.06)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  ecoli_mash_all <-remove_outliers(ecoli_mash_all, outl)
#  ecoli_accnet_all <-remove_outliers(ecoli_accnet_all, outl)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  files <-  anti_join(files, outl, by=c("V1"="Source"))

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  # For select 800 non redundant samples
#  nr_list <- non_redundant(ecoli_mash_all,number = 800)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  # For select a subset of samples with 99.99% of indentity
#  nr_list <- non_redundant(ecoli_mash_all, distance = 0.0001)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  ecoli_accnet <- extract_non_redundant(ecoli_accnet_all, nr_list)
#  ecoli_mash <- extract_non_redundant(ecoli_mash_all, nr_list)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  nr_list <- non_redundant_hier(ecoli_mash_all,800, partitions = 10)
#  ecoli_accnet <- extract_non_redundant(ecoli_accnet_all, nr_list)
#  ecoli_mash <- extract_non_redundant(ecoli_mash_all, nr_list)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  cp <- core_plots(ecoli_mm,reps = 10, threshold = 0.95, steps = 10)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  nr_files = nr_list %>%
#    as.data.frame() %>%
#    group_by(cluster) %>%
#    top_n(1,centrality) %>%
#    summarise_all(first) %>%
#    select(Source) %>%
#    distinct()

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  core <- mmseqs(nr_files) %>% core_genome(type = "prot")
#  export_core_to_fasta(core,"core.aln")
#  

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  library(phangorn)
#  library(ggtree)
#  core_tree = read.tree("core.tree")
#  
#  annot_tree = species %>%
#    filter(!grepl("Citrobacter",organism_name)) %>%
#    filter(!grepl("marmotae",organism_name)) %>%
#    select(Target,organism_name) %>% distinct()
#  core_tree %>% midpoint %>% ggtree(layout = "circular") %<+% annot_tree + geom_tippoint(aes(color = organism_name))

## ----eval = FALSE, include=TRUE-----------------------------------------------
#  var_matrix = core_snps_matrix(core, norm = TRUE)
#  pheatmap::pheatmap(var_matrix,show_rownames = F, show_colnames = F)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  ecoli_accnet_all <- accnet(ecoli_mm,threshold = 0.8, singles = FALSE)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  ec_nr<- non_redundant(ecoli_mash,number = 800 )
#  ecoli_accnet_nr <- extract_non_redundant(ecoli_accnet, ef_nr)
#  ec_800_cl <- clustering(ecoli_accnet_nr, method = "mclust", d_reduction = TRUE)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  export_to_gephi(ecoli_accnet_nr, "accnet800", cluster = ec_800_cl)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  accnet_enr_result <- accnet_enrichment_analysis(ecoli_accnet_nr, cluster = ef_800_cl)
#  
#  # A tibble: 1,149,041 x 14
#     Target             Source                   Cluster perClusterFreq ClusterFreq ClusterGenomeSize perTotalFreq TotalFreq OdsRatio   pvalue    padj AccnetGenomeSize AccnetProteinSi… Annot
#     <chr>              <chr>                      <dbl>          <dbl>       <int>             <int>        <dbl>     <int>    <dbl>    <dbl>   <dbl>            <int>            <int> <chr>
#   1 1016|NZ_CP053592.… GCF_000009565.1_ASM956v…       1         0.0436          13               298       0.0173        20     2.52  3.62e-5 0.00130             1156            74671 ""
#   2 1016|NZ_CP053592.… GCF_000022665.1_ASM2266…       1         0.0436          13               298       0.0173        20     2.52  3.62e-5 0.00130             1156            74671 ""
#   3 1016|NZ_CP053592.… GCF_000023665.1_ASM2366…       1         0.0436          13               298       0.0173        20     2.52  3.62e-5 0.00130             1156            74671 ""
#   4 1016|NZ_CP053592.… GCF_000830035.1_ASM8300…       1         0.0436          13               298       0.0173        20     2.52  3.62e-5 0.00130             1156            74671 ""
#   5 1016|NZ_CP053592.… GCF_000833145.1_ASM8331…       1         0.0436          13               298       0.0173        20     2.52  3.62e-5 0.00130             1156            74671 ""
#   6 1016|NZ_CP053592.… GCF_001039415.1_ASM1039…       1         0.0436          13               298       0.0173        20     2.52  3.62e-5 0.00130             1156            74671 ""
#   7 1016|NZ_CP053592.… GCF_001596115.1_ASM1596…       1         0.0436          13               298       0.0173        20     2.52  3.62e-5 0.00130             1156            74671 ""
#   8 1016|NZ_CP053592.… GCF_009663855.1_ASM9663…       1         0.0436          13               298       0.0173        20     2.52  3.62e-5 0.00130             1156            74671 ""
#   9 1016|NZ_CP053592.… GCF_009832985.1_ASM9832…       1         0.0436          13               298       0.0173        20     2.52  3.62e-5 0.00130             1156            74671 ""
#  10 1016|NZ_CP053592.… GCF_013166955.1_ASM1316…       1         0.0436          13               298       0.0173        20     2.52  3.62e-5 0.00130             1156            74671 ""
#  # … with 1,149,031 more rows

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  accnet_with_padj(accnet_enr_result) %>% export_to_gephi("accnet800.padj", cluster = ec_800_cl)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  mash_tree_fastme <- similarity_tree(ecoli_mash)
#  mash_tree_NJ <- similarity_tree(ecoli_mash, method = "NJ")
#  mash_tree_upgma <- similarity_tree(ecoli_mash,method = "UPGMA")
#  accnet_tree_upgma <- similarity_tree(ecoli_accnet,method = "UPGMA")

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  mash_tree_fastme %>% midpoint %>% ggtree(layout = "circular") %<+% annot_tree + geom_tippoint(aes(color = organism_name))

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  library(dendextend)
#  tanglegram(ladderize(mash_tree_upgma), ladderize(accnet_tree_upgma), fast = TRUE, main = "Mash vs AccNET UPGMA")
#  

## ----eval=FALSE, include=FALSE------------------------------------------------
#  mmseqs(files) %>% accnet() %>% export_accnet_aln(.,file ="my_accnet_aln.fasta")

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  acc_tree = read.tree("acc.aln.treefile")
#  acc_tree %>% midpoint.root() %>% ggtree()

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  
#  # K-NNN with 10 neighbours and with repetitions
#  knnn_mash_10_w_r <- knnn(ecoli_mash,n=10, repeats = TRUE)
#  # K-NNN with 25 neighbours and with repetitions
#  knnn_mash_25_w_r <- knnn(ecoli_mash,n=25, repeats = TRUE)
#  # K-NNN with 50 neighbours and with repetitions
#  knnn_mash_50_w_r <- knnn(ecoli_mash,n=50, repeats = TRUE)
#  

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  export_to_gephi(knnn_mash_50_w_r,file = "knnn_50_w_r.tsv")

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  plot_knnn_network(knnn_mash_50_w_r)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  ec_cl_mclust_umap <- clustering(ecoli_mash, method = "mclust",d_reduction = TRUE)

## ----eval =FALSE, include=TRUE------------------------------------------------
#  ec_cl_knnn <-clustering(knnn_mash_50_w_r, method = "louvain")

## ----eval =FALSE, include=TRUE------------------------------------------------
#  umap_mash <- umap_plot(ecoli_mash)

## ----eval =FALSE, include=TRUE------------------------------------------------
#  umap_mash <- umap_plot(ecoli_mash, cluster = ef_cl_mclust_umap)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  cl_louvain = clustering(knnn_mash_25_wo_r, method = "louvain")
#  plot_knnn_network(knnn_mash_25_wo_r, cluster = cl_louvain, edge.alpha = 0.1)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  library(tidyverse)
#  
#  annotation <- annotate(files, type = "prot",database = c("AbR","VF_A"))
#  
#  annotation %>%
#    filter(pident > 0.95 ) %>% #remove all hits with identity lower than 95%
#    filter(evalue < 1e-6) %>%  #remove all hits with E-Value greater than 1e-6
#    group_by(Genome,Protein) %>% arrange(desc(bits)) %>% slice_head() #select only the best hit for each protein-genome.
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
#    #mutate(file = basename(path)) %>%
#    separate(path,c("GCF","acc","acc2"),sep = "_", remove = FALSE) %>%
#    unite(assembly_accession,GCF,acc,sep = "_") %>%
#    left_join(assembly) %>%
#    separate(organism_name,c("genus","specie"), sep = " ") %>%
#    group_by(pangenome,genus,specie) %>%
#    summarise(N = n()) %>%
#    distinct() %>% ungroup() %>%
#    group_by(pangenome) %>%
#    slice_head() %>%
#    mutate(ID = paste("pangenome_",pangenome,"_rep_seq.fasta", sep = "",collapse = ""))
#  
#  annot <- annot %>%
#    mutate(genus = gsub("\\[","",genus)) %>%
#    mutate(genus = gsub("\\]","",genus)) %>%
#    mutate(specie = gsub("\\[","",specie)) %>%
#    mutate(specie = gsub("\\]","",specie)) %>%
#    unite("Label",genus,specie, remove = F)
#  
#  write_delim(annot,"../pangenomes_gephi_extra_annot.tsv",delim = "\t", col_names = TRUE)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  sh <-  sharedness(res)
#  sh  <-  sh %>% as.data.frame() %>%
#    rownames_to_column("ID") %>%
#    inner_join(annot) %>%
#    unite(Name,genus,specie,pangenome, sep = "_") %>%
#    unite("Row_n",Label,N) %>%
#    select(-ID)%>%
#    column_to_rownames("Row_n")
#  
#  colnames(sh) = rownames(sh)
#  
#  

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  pheatmap::pheatmap(log2(sh+1),
#                     clustering_method = "ward.D2",
#                     clustering_distance_cols = "correlation",
#                     clustering_distance_rows = "correlation")

