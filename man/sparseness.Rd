% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R
\name{sparseness}
\alias{sparseness}
\title{Test graph for sparseness}
\usage{
sparseness(G, threshold = 0.1)
}
\arguments{
\item{G}{igraph object representing the graph to test.}

\item{threshold}{numeric, between 0 and 1 indicating the proportion of which
to consider much greater than. For example, if 10 is much greater than 1,
then the threshold is considered to be 0.1. Default threshold is 0.1.}
}
\description{
Test graph for sparseness
}
\author{
Christopher Nobles, Ph.D.
}
