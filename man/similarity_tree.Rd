% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/similarity_tree.R
\name{similarity_tree}
\alias{similarity_tree}
\title{Similarity tree}
\usage{
similarity_tree(data, method = "NJ")
}
\arguments{
\item{data}{An \emph{accnet} or \emph{mash} object.}

\item{method}{Tree method c("fastme", "NJ", "UPGMA"):
\itemize{
\item{fastme: Tree Estimation Based on the Minimum Evolution Algorithm}
\item{NJ: Neighbor-Joining Tree Estimation}
\item{UPGMA: Unweighted pair group method with arithmetic mean}
\item{ward.D2: Ward's minimum variance metho}
\item{median: WPGMC Weighted Pair Group Method with Arithmetic Mean}
\item{complete: complete linkage  }
}}
}
\value{
An object of class \emph{phylo}
}
\description{
Build a similarity tree (pseudo-phylogenetic) of the samples,
using whole-genome (\emph{mash}) or accessory genome (\emph{accnet}).
}
\seealso{
\code{\link[ape]{fastme}}

\code{\link[ape]{nj}}

\code{\link[stats]{hclust}}
}
