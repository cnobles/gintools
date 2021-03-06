% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/determine_abundance.R
\name{determine_abundance}
\alias{determine_abundance}
\title{Determine the abundance of integration sites within a GRanges object.}
\usage{
determine_abundance(sites)

determine_abundance(sites, grouping = NULL, replicates = NULL,
method = "fragLen")
}
\arguments{
\item{sites}{GRanges object with flanking genomic DNA ranges.}

\item{grouping}{Character name of the metadata column with grouping
information.}

\item{replicates}{Character name of the metadata column which dictates
replicate information.}

\item{method}{Character either "fragLen" for determining abundance by unique
fragment length or "estAbund" to use the sonicLength package.}
}
\description{
Given a GRanges object containing the genomic ranges of
flanking genomic DNA from integration sites, the function returns abundance
information based on unique widths for each site (method = "fragLen",
default) or applies a more powerful tool developed in sonicLength package
(method = "estAbund"). If multiple samples are given, they can be grouped
separately by giving the metadata column to the grouping option.
}
\details{
\code{determine_abundance} returns the abundance of all distinct integration
sites within the given GRanges object based on the unique fragment length
(width) or utilizes the sonicLength package.
}
\examples{
gr <- gintools:::generate_test_granges()
std.gr <- standardize_sites(gr)

determine_abundance(std.gr)

}
\author{
Christopher Nobles, Ph.D.
}
