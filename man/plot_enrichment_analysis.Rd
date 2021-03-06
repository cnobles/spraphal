% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_enrichment_analysis}
\alias{plot_enrichment_analysis}
\title{Plot enrichment analysis of ordered verticies.}
\usage{
plot_enrichment_analysis(v, G, d = 50, f = 100, e = NULL, bLim = 10,
  critVal = 0.05, exact = TRUE, ...)
}
\arguments{
\item{v}{vector of vertex identifiers. All vertex identifiers need to be
present in the graph G.}

\item{G}{igraph object of which v vertices are a part of and is considered
the background for enrichment analysis.}

\item{d}{delta distance between analysis points. Default 50.}

\item{f}{forward distance for extended analysis. Default 100.}

\item{e}{edge distance for which to consider, from the beginning. Default d.}

\item{bLim}{bottom limit for significance value, positive log of p-value. ie.
10 -> minimum p-value of 10^-10 plotted.}

\item{critVal}{critical value for which to consider significance. Default 
0.05.
#' @param exact logical. If TRUE, the approximation of the probability is
calculated rather than approximated by `ppois()`.}

\item{...}{arguments to be passed into ggplot.}
}
\description{
Plot enrichment analysis of ordered verticies.
}
\author{
Christopher Nobles, Ph.D.
}
