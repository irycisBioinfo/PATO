# Pangenome Analysis Toolkit (PATO)

<img src="https://github.com/irycisBioinfo/PATO/blob/master/vignettes/logo.png" width="50%" align = "right">


PATO is a R package designed to analyze pangenomes (set of genomes)
intra or inter species. It allows to analyze the core-genome, accessory 
genome and whole genome, the population structure, and the horizontal 
gene transfer dynamics. PATO uses, as core software, 
[MASH](https://mash.readthedocs.io/en/latest/) , 
[MMSeq2](https://github.com/soedinglab/MMseqs2),
[Minimap2](https://github.com/lh3/minimap2) and R.

These software can handle thousands of genomes using conventional computers      
without the necessity to use on a HPC facilities. PATO can analyze data 
in mash distance format (whole genome) or accnet format (accessory genome)
Most of the functions can handle both objets. Some functions are specific
for some kind of data. The primitive objet *mmseq* is a orthologous clustering that
is used to build accnet object, annotate the genomes and characterize 
the core-, accessory- and pan-genome size (and dynamic).

<img src="https://github.com/irycisBioinfo/PATO/blob/master/vignettes/diagram.png" width="100%">

## Installation

PATO requires R 3.6 or newer. 
To install PATO package you need to install [devtools](https://github.com/r-lib/devtools) package

```
install.packages("devtools")
```
You need to activate Bioconductor reposiroty.
```
setRepositories()
##Then select 1 and 2.
```

Once you have installed *devtools* package and activate Bioconductor repository, type:

```
devtools::install_github("https://github.com/irycisBioinfo/PATO", build_vignettes = T)
```

Some times some dependencies require system packages as `libcurl`. `libssl` or `libxml2` (this example is for Ubuntu/Debian based systems):

```
sudo apt install libcurl4-openssl-dev libssl-dev libxml2-dev libmagick++-dev libv8-dev 
```

It could take a while because PATO has the following dependencies:

```
 ggplot2 (>= 3.1.1),
 tibble (>= 2.1.3),
 htmlwidgets (>= 1.3),
 htmltools (>= 0.4.0),
 threejs (>= 0.3.1),
 igraph (>= 1.2.4.1),
 uwot (>= 0.1.3),
 tidyr (>= 0.8.3),
 dplyr (>= 1.0.0),
 dtplyr (>= 0.0.3),
 data.table (>= 1.12.2),
 phangorn (>= 2.5.5),
 randomcoloR (>= 1.1.0),
 magrittr (>= 1.5),
 Rcpp (>= 1.0.2),
 mclust (>= 5.4.5),
 dbscan (>= 1.1-4),
 ape (>= 5.3),
 parallelDist (>= 0.2.4),
 foreach (>= 1.5.0),
 doParallel (>= 1.0.15),
 parallel (>= 3.5.2),
 manipulateWidget (>= 0.10.0),
 stringr (>= 1.4.0),
 stringdist (>= 0.9.6),
 openssl (>= 1.4.3),
 Biostrings (>= 2.52.0),
 microseq (>= 2.1.2)
```

# Documentation
You can find the documentation in the [wiki](https://github.com/irycisBioinfo/PATO/wiki) area


