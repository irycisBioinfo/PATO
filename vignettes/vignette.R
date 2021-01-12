## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  install.packages("devtools")

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  setRepositories()
#  ## Select options 1 (CRAN) and 2 (Bioconductor Software)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  devtools::install_github("https://github.com/irycisBioinfo/PATO.git", build_vignettes = TRUE)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  install.packages("devtools")

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  setRepositories()
#  ## Select options 1 (CRAN) and 2 (Bioconductor Software)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  devtools::install_github("https://github.com/irycisBioinfo/PATO.git", build_vignettes = TRUE)

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
#  files <- dir(pattern = ".faa", full.names = T)

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
#  
#  setwd("~/examplePATO/)
#  

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  
#  gff_files <- dir("~/examplesPATO/", pattern = "\\.gff", full.names =T)  ##Creates a file-list
#  gffs <- load_gff_list(gff_files)                                        ##Load the genomes
#  

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  strain_names <- read.table("names.txt")%>%                #Load the files
#    rename(Genome = V1, Name = V2) %>%                      #Rename the columns
#    mutate(Name = gsub("-","",Name)) %>%                    #delete the '-' character
#    mutate(Sample = str_sub(Name,1,4)) %>%                  #Extract the first 4 character as Sample name
#    mutate(Source = str_sub(Name,5)) %>%                    #Extract the final characters as Source
#    mutate(Genome = str_replace(Genome,"_genomic.gff",""))  #delete the final part of the filename

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  mash <- mash(gffs, type ="wgs", n_cores = 20)
#  mash_tree <- similarity_tree(mash,method = "fastme") %>% phangorn::midpoint()
#  

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  # remove '_genomic.fna' from the tip.labels to fix with the strain_names table
#  mash_tree$tip.label <- gsub("_genomic.fna","",mash_tree$tip.label)
#  ggtree(mash_tree) %<+% strain_names + geom_tippoint(aes(color=Source)) + geom_tiplab(aes(label = Sample))
#  

## ----eval=F, include=T--------------------------------------------------------
#  mash$table %>%
#    mutate(Source = gsub("_genomic.fna","",Source)) %>%   #Remove extra characters
#    mutate(Target = gsub("_genomic.fna","",Target)) %>%
#    inner_join(strain_names %>%                           #Add real names to Source
#                 rename(Source2 = Name) %>%
#                 select(Genome,Source2),
#               by= c("Source" = "Genome")) %>%
#    inner_join(strain_names %>%                           #Add real names to Target
#                 rename(Target2 = Name) %>%
#                 select(Genome,Target2),
#               by= c("Target" = "Genome")) %>%
#    filter(Source != Target) %>%                          #Remove diagonal
#    mutate(ANI = (1-Dist)*100) %>%                        #Compute the ANI
#    filter(ANI > 99.9) %>%                                #Filter hits higher than 99.9% identity
#    as.data.frame()                                       #Transform to data.frame to show all digits of ANI
#  
#    Source2 Target2      ANI
#  1   HH08H   HH29S 99.94674
#  2  HH15CH  HH29CH 99.91693
#  3  HH29CH  HH15CH 99.91693
#  4   HH29S   HH08H 99.94674
#  

## ----eval = FALSE, include=TRUE-----------------------------------------------
#  
#  mm <- mmseqs(gffs, type = "nucl")
#  core <- core_genome(mm,type = "nucl", n_cores = 20)
#  export_core_to_fasta(core,file = "pato_roary_like.fasta")
#  
#  #Externaly compute the phylogenetic tree
#  system("fasttreeMP -nt -gtr pato_roary_like.fasta > pato_roary_like.tree")
#  
#  #Import and rooting the tree
#  pato_roary_tree <- ape::read.tree("pato_roary_like.tree") %>% phangorn::midpoint()
#  #Fix the tip labels
#  pato_roary_tree$tip.label <- gsub("_genomic.ffn","",pato_roary_tree$tip.label)
#  
#  #Draw the tree with the annotation.
#  ggtree(pato_roary_tree) %<+% strain_names + geom_tippoint(aes(color=Source)) + geom_tiplab(aes(label = Sample))

## ----eval=F, include=T--------------------------------------------------------
#  
#  pato_roary_m = core_snps_matrix(core, norm = T, rm.gaps = T) #compute the SNP matrix removing columns with gaps
#  
#  #Fix the colnames and rownames to put the real names.
#  
#  colnames(pato_roary_m) <- rownames(pato_roary_m) <- gsub("_genomic.ffn","",rownames(pato_roary_m))
#  tmp = pato_roary_m %>%
#    as.data.frame() %>%
#    rownames_to_column("Genome") %>%
#    inner_join(strain_names) %>%
#    column_to_rownames("Name") %>%
#    select(-Genome,-Source,-Sample) %>%
#    as.matrix()
#  colnames(tmp) = rownames(tmp)
#  
#  #Plot the matrix as a heatmap
#  pheatmap::pheatmap(tmp,
#           display_numbers = matrix(ifelse(tmp < 200,tmp,""), nrow(tmp)), #To Display only values lower than 200 SNPs
#           number_format = '%i',
#           annotation_row = strain_names %>% select(-Genome,-Sample) %>% column_to_rownames("Name"),
#           annotation_col = strain_names %>% select(-Genome,-Sample) %>% column_to_rownames("Name"),
#           show_rownames = T,
#           show_colnames = T,
#  )
#  

## ----eval=F, include=T--------------------------------------------------------
#  tmp %>%
#    as.data.frame() %>%
#    rownames_to_column("Source") %>%
#    pivot_longer(-Source,names_to = "Target", values_to = "SNPs") %>%
#    filter(Source != Target) %>%
#    filter(SNPs < 200)
#  
#  # A tibble: 10 x 3
#     Source Target  SNPs
#     <chr>  <chr>  <dbl>
#   1 HH29S  HH08H      1
#   2 HH29CH HH15CH    21
#   3 HH24H  HH24CH     0
#   4 HH24CH HH24H      0
#   5 HH24C  HH20C     72
#   6 HH20C  HH24C     72
#   7 HH19S  HH19CH     1
#   8 HH19CH HH19S      1
#   9 HH15CH HH29CH    21
#  10 HH08H  HH29S      1

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  
#  core_s <- core_snp_genome(gffs,type = "wgs")
#  export_core_to_fasta(core_s,file = "pato_snippy_like.fasta")
#  #Externaly compute the phylogenetic tree
#  system("fasttreeMP -nt -gtr pato_snippy_like.fasta > pato_snippy_like.tree")
#  
#  #Import and rooting the tree
#  pato_snippy_tree <- ape::read.tree("pato_snippy_like.tree") %>% phangorn::midpoint()
#  #Fix the tip labels
#  pato_snippy_tree$tip.label <- gsub("_genomic.fna","",pato_snippy_tree$tip.label)
#  
#  #Draw the tree with the annotation.
#  ggtree(pato_snippy_tree) %<+% strain_names + geom_tippoint(aes(color=Source)) + geom_tiplab(aes(label = Sample))
#  

## ----eval = F, include=T------------------------------------------------------
#  
#  colnames(pato_snippy_m) <- rownames(pato_snippy_m) <- gsub("_genomic.fna","",rownames(pato_snippy_m))
#  tmp = pato_snippy_m %>%
#    as.data.frame() %>%
#    rownames_to_column("Genome") %>%
#    inner_join(strain_names) %>%
#    column_to_rownames("Name") %>%
#    select(-Genome,-Source,-Sample) %>%
#    as.matrix()
#  
#  
#  
#  colnames(tmp) = rownames(tmp)
#  
#  pheatmap(tmp,
#           display_numbers = matrix(ifelse(tmp < 200,tmp,""), nrow(tmp)),
#           number_format = '%i',
#           annotation_row = strain_names %>% select(-Genome,-Sample) %>% column_to_rownames("Name"),
#           annotation_col = strain_names %>% select(-Genome,-Sample) %>% column_to_rownames("Name"),
#           show_rownames = T,
#           show_colnames = T,
#           )
#  

## ----eval = F, include=T------------------------------------------------------
#  tmp %>%
#    as.data.frame() %>%
#    rownames_to_column("Source") %>%
#    pivot_longer(-Source,names_to = "Target", values_to = "SNPs") %>%
#    filter(Source != Target) %>%
#    filter(SNPs < 200)
#  
#  # A tibble: 12 x 3
#     Source Target  SNPs
#     <chr>  <chr>  <dbl>
#   1 HH29S  HH08H      8
#   2 HH29CH HH15CH    28
#   3 HH24H  HH24CH     8
#   4 HH24CH HH24H      8
#   5 HH24C  HH20C    144
#   6 HH20C  HH24C    144
#   7 HH19S  HH19CH    11
#   8 HH19CH HH19S     11
#   9 HH16CH HH03H     91
#  10 HH15CH HH29CH    28
#  11 HH08H  HH29S      8
#  12 HH03H  HH16CH    91

## ----eval=F, include=T--------------------------------------------------------
#  
#  pw_normalized <- snps_pairwaise(gffs, type = "wgs",norm = T, n_cores = 20)
#  
#  tmp = pw_normalized %>%
#    as.data.frame() %>%
#    rownames_to_column("Genome") %>%
#    inner_join(strain_names) %>%
#    column_to_rownames("Name") %>%
#    select(-Genome,-Source,-Sample) %>%
#    as.matrix()
#  colnames(tmp) = rownames(tmp)
#  
#  pheatmap(tmp,
#           display_numbers = matrix(ifelse(tmp < 200,tmp,""), nrow(tmp)),
#           number_format = '%i',
#           annotation_row = strain_names %>% select(-Genome,-Sample) %>% column_to_rownames("Name"),
#           annotation_col = strain_names %>% select(-Genome,-Sample) %>% column_to_rownames("Name"),
#           show_rownames = T,
#           show_colnames = T,
#  )
#  
#  
#  

## ----eval=F, include=T--------------------------------------------------------
#  tmp %>%
#    as.data.frame() %>%
#    rownames_to_column("Source") %>%
#    pivot_longer(-Source,names_to = "Target", values_to = "SNPs") %>%
#    filter(Source != Target) %>%
#    filter(SNPs < 200)
#  
#  # A tibble: 10 x 3
#     Source Target  SNPs
#     <chr>  <chr>  <dbl>
#   1 HH29S  HH08H      8
#   2 HH29CH HH15CH    27
#   3 HH24H  HH24CH     2
#   4 HH24CH HH24H      2
#   5 HH24C  HH20C     97
#   6 HH20C  HH24C     97
#   7 HH19S  HH19CH    13
#   8 HH19CH HH19S     13
#   9 HH15CH HH29CH    27
#  10 HH08H  HH29S      8
#  

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

## ----eval=F, include=T--------------------------------------------------------
#  
#  library(microbenchmark)
#  
#  
#  bench <- microbenchmark(
#  
#    "roary_time" = {system("roary -p 24 -o roary_ex -e --mafft *.gff")},
#  
#    "pato_like_roary_time"= {
#      gffs <- load_gff_list(gff_files)
#      mm <- mmseqs(gffs, type = "nucl")
#      acc <- accnet(mm)
#      core <- core_genome(mm,type = "nucl", n_cores = 24)
#      core_p <- core_plots(mm,steps = 10,reps = 5)
#      export_accnet_aln(file = "accesory_aln.fa",acc)
#      system("fasttreeMP -nt accesory_aln.fa > accessory_aln.tree")
#    },
#  
#    "pato_like_snippy_time" = {
#      unlink(gffs$path,recursive = T)
#      gffs <- load_gff_list(gff_files)
#      core_s <- core_snp_genome(gffs,type = "wgs", ref = "~/examplesPATO/fna/GCF_009820385.1_ASM982038v1_genomic.fna")
#      export_core_to_fasta(core_s,file = "pato_snippy_like.fasta")
#    },
#  
#    "snippy_time" = {
#      setwd("~/examplesPATO/fna")
#  
#      # That is my $PATH. Do not try to execute. To execute Snippy from R you need to define the $PATH
#      Sys.setenv(PATH ="/usr/local/bin:/usr/bin:/bin:/home/val/bioinfo/snippy/bin/:/home/val/bioinfo/ ....")
#  
#      system("ls *.fna > tmp2")
#      system("ls $PWD/*.fna > tmp")
#      system("sed -i 's/_genomic.fna//' tmp")
#      system("paste tmp tmp2 > list_snippy.txt")
#      system("snippy-multi list_snippy.txt --ref GCF_009820385.1_ASM982038v1_genomic.fna --cpus 24 > snp.sh")
#      system("sh ./snp.sh")
#    },
#    times = 1
#  )
#  
#  

## ----eval=F, include=T--------------------------------------------------------
#  library(ape)
#  library(TreeDist)
#  library(phangorn)
#  
#  pato_snippy_tree = read.tree("pato_snippy_like.tree") %>% midpoint()
#  pato_snippy_tree$tip.label = gsub("_genomic.fna","",pato_snippy_tree$tip.label)
#  
#  pato_roary_tree = read.tree("pato_roary_like.tree") %>% midpoint()
#  pato_roary_tree$tip.label = gsub("_genomic.ffn","",pato_roary_tree$tip.label)
#  
#  
#  snippy_tree = read.tree("snippy.tree") %>% midpoint()
#  snippy_tree = TreeTools::DropTip(snippy_tree,"Reference")
#  
#  roary_tree = read.tree("roary.tree") %>% midpoint()
#  roary_tree$tip.label = gsub("_genomic","",roary_tree$tip.label)
#  
#  rand_tree = rtree(60, tip.label = snippy_tree$tip.label) %>% midpoint()
#  
#  
#  mash <- mash(gffs, type ="wgs", n_cores = 20)
#  mash_tree <- similarity_tree(mash,method = "fastme") %>% midpoint()
#  
#  mash_tree$tip.label = gsub("_genomic.fna","",mash_tree$tip.label)
#  
#  trees <- c(random_tree = rand_tree, snippy_tree = snippy_tree,roary_tree = roary_tree,pato_core = pato_roary_tree, pato_core_snp = pato_snippy_tree, mash = mash_tree)
#  
#  CompareAll(trees,RobinsonFoulds)%>% as.matrix() %>% pheatmap::pheatmap(.,display_numbers = T)
#  

## ----eval=F, include=T--------------------------------------------------------
#  
#  > sessionInfo()
#  R version 3.6.3 (2020-02-29)
#  Platform: x86_64-pc-linux-gnu (64-bit)
#  Running under: Debian GNU/Linux 10 (buster)
#  
#  Matrix products: default
#  BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.8.0
#  LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.8.0
#  
#  locale:
#   [1] LC_CTYPE=es_ES.UTF-8       LC_NUMERIC=C               LC_TIME=es_ES.UTF-8        LC_COLLATE=es_ES.UTF-8     LC_MONETARY=es_ES.UTF-8
#   [6] LC_MESSAGES=es_ES.UTF-8    LC_PAPER=es_ES.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C
#  [11] LC_MEASUREMENT=es_ES.UTF-8 LC_IDENTIFICATION=C
#  
#  attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base
#  
#  other attached packages:
#   [1] pheatmap_1.0.12      stringdist_0.9.6.3   microseq_2.1.2       rlang_0.4.8          data.table_1.13.2    TreeDist_1.2.1
#   [7] microbenchmark_1.4-7 phangorn_2.5.5       ape_5.4-1            ggtree_2.0.4         pato_1.0.1           forcats_0.5.0
#  [13] stringr_1.4.0        dplyr_1.0.2          purrr_0.3.4          readr_1.3.1          tidyr_1.1.2          tibble_3.0.4
#  [19] ggplot2_3.3.2        tidyverse_1.3.0
#  
#  loaded via a namespace (and not attached):
#    [1] Rtsne_0.15              colorspace_1.4-1        ellipsis_0.3.1          mclust_5.4.6            XVector_0.26.0          base64enc_0.1-3
#    [7] fs_1.5.0                rstudioapi_0.11         farver_2.0.3            bit64_4.0.5             fansi_0.4.1             lubridate_1.7.9
#   [13] xml2_1.3.2              codetools_0.2-16        R.methodsS3_1.8.1       doParallel_1.0.16       jsonlite_1.7.1          broom_0.7.0
#   [19] cluster_2.1.0           dbplyr_1.4.4            R.oo_1.24.0             uwot_0.1.8              shiny_1.5.0             BiocManager_1.30.10
#   [25] compiler_3.6.3          httr_1.4.2              rvcheck_0.1.8           backports_1.1.10        lazyeval_0.2.2          assertthat_0.2.1
#   [31] Matrix_1.2-18           fastmap_1.0.1           cli_2.1.0               later_1.1.0.1           htmltools_0.5.0         tools_3.6.3
#   [37] igraph_1.2.6            gtable_0.3.0            glue_1.4.2              V8_3.2.0                fastmatch_1.1-0         Rcpp_1.0.5
#   [43] cellranger_1.1.0        vctrs_0.3.4             Biostrings_2.54.0       nlme_3.1-144            crosstalk_1.1.0.1       iterators_1.0.13
#   [49] gbRd_0.4-11             rbibutils_2.0           openxlsx_4.2.2          rvest_0.3.6             mime_0.9                miniUI_0.1.1.1
#   [55] lifecycle_0.2.0         zlibbioc_1.32.0         scales_1.1.1            hms_0.5.3               promises_1.1.1          parallel_3.6.3
#   [61] RColorBrewer_1.1-2      curl_4.3                memoise_1.1.0           dtplyr_1.0.1            stringi_1.5.3           randomcoloR_1.1.0.1
#   [67] S4Vectors_0.24.4        foreach_1.5.1           tidytree_0.3.3          BiocGenerics_0.32.0     zip_2.1.1               manipulateWidget_0.10.1
#   [73] Rdpack_2.1              pkgconfig_2.0.3         parallelDist_0.2.4      lattice_0.20-40         labeling_0.4.2          treeio_1.10.0
#   [79] htmlwidgets_1.5.2       bit_4.0.4               tidyselect_1.1.0        magrittr_1.5            R6_2.4.1                IRanges_2.20.2
#   [85] generics_0.0.2          DBI_1.1.0               pillar_1.4.6            haven_2.3.1             withr_2.3.0             modelr_0.1.8
#   [91] crayon_1.3.4            utf8_1.1.4              grid_3.6.3              readxl_1.3.1            blob_1.2.1              threejs_0.3.3
#   [97] reprex_0.3.0            digest_0.6.26           webshot_0.5.2           xtable_1.8-4            R.cache_0.14.0          httpuv_1.5.4
#  [103] dbscan_1.1-5            R.utils_2.10.1          TreeTools_1.4.1         openssl_1.4.3           RcppParallel_5.0.2      stats4_3.6.3
#  [109] munsell_0.5.0           askpass_1.1             quadprog_1.5-8

