% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statistics.R
\name{comp_vertex_prob}
\alias{comp_vertex_prob}
\title{Approximate the probability of each vertex within a list belonging to the
list.}
\usage{
comp_vertex_prob(v, G, exact = TRUE)
}
\arguments{
\item{v}{vector of vertex identifiers. Vertices must be present in graph G.}

\item{G}{igraph object representing the graph and including vertices from v.}

\item{exact}{logical. If TRUE, the approximation of the probability is
calculated rather than approximated by `ppois()`.}
}
\description{
Approximate the probability of each vertex within a list belonging to the
list.
}
\author{
Christopher Nobles, Ph.D.
}
