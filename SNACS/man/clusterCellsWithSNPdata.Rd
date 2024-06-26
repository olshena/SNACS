% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeSNACScall.R
\name{clusterCellsWithSNPdata}
\alias{clusterCellsWithSNPdata}
\title{Cluster cells based on SNP data}
\usage{
clusterCellsWithSNPdata(
  snacsObj,
  outputFormat = c("none", "", "pdf", "png"),
  verbose = FALSE
)
}
\arguments{
\item{snacsObj}{SNACSList object}

\item{outputFormat}{Character. Output file type. "" outputs to the standard output. Default is "none" which does not plot anything}

\item{verbose}{Logical. Prints information when running the method. Default is FALSE}
}
\value{
A SNACSList object
}
\description{
Cluster cells, where the number of clusters is the number of constituent samples, using the mutation data of select SNPs. Generate heatmap of the mutation data.
}
