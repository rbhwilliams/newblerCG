% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rng.R
\name{make.rng.from.2d.binning.space}
\alias{make.rng.from.2d.binning.space}
\title{Make an RNG from two bininng variables}
\usage{
make.rng.from.2d.binning.space(bM)
}
\arguments{
\item{bM}{is a matrix with rows indexing contigs and 2 columns that index binning variables (e.g. two coverage vectors or one coverage and one GC etc). For single-sample assemblies it is generated from \code{prep.for.rbg}}
}
\value{
An \code{igraph} containing an RNG
}
\description{
This function generates a relative neighbourhood graph from two binning variables.
}


