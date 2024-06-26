% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processData.R
\name{imputeMissingMutations}
\alias{imputeMissingMutations}
\title{Impute missing mutations in SNACSList object.}
\usage{
imputeMissingMutations(snacsObj, verbose = FALSE)
}
\arguments{
\item{snacsObj}{SNACSList object}

\item{verbose}{Logical. Prints information when running the method. Default is FALSE}
}
\value{
A SNACSList object
}
\description{
Missing values in the mutation data are imputed using K-nearest neighbors (KNN) method.
}
