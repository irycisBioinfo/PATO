% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_knnn.R
\name{cluster_knnn}
\alias{cluster_knnn}
\title{Internal function to cluster igraph data}
\usage{
cluster_knnn(data, method)
}
\arguments{
\item{data}{\emph{igraph} object}

\item{method}{greedy, louvain or walktrap}
}
\value{
\emph{data.frame} with two columns: \emph{Source} and \emph{Target}
}
\description{
Internal function to cluster igraph data
}
