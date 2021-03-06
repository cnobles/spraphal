% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statistics.R
\name{comp_enrichment}
\alias{comp_enrichment}
\title{Calculate the enrichment within an ordered list of verticies.}
\usage{
comp_enrichment(v, G, d = 50, f = 100, e = NULL, exact = TRUE)
}
\arguments{
\item{v}{vector of vertex identifiers. Vertices must be present in graph G.}

\item{G}{igraph object representing the graph and including vertices from v.}

\item{d}{delta or change distance for which to look for enrichment. Default
20.}

\item{f}{forward threshold for reference enrichment. Default 100.}

\item{e}{edge window to consider at the beginning. If not specified, the
default value will be equal to the delta `d`.}

\item{exact}{logical. If TRUE, the approximation of the probability is
calculated rather than approximated by `ppois()`.}
}
\description{
Calculate the enrichment within an ordered list of verticies.
}
\author{
Christopher Nobles, Ph.D.
}
