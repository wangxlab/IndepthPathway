% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PathwayAssociation.R
\name{pathwayAssociation}
\alias{pathwayAssociation}
\title{Title Calculation the functional associations of deambiguated pathway enrichment result}
\usage{
pathwayAssociation(
  DeambiguatePathway,
  compare.list,
  feature.list,
  preCalmatrix,
  selectn = 30,
  minsize = 10,
  rm.overfit = FALSE
)
}
\arguments{
\item{DeambiguatePathway}{up or down.deambiguation data}

\item{compare.list}{Gene Signatures provided by Molecular Signature Database (MSigDB)}

\item{feature.list}{pre-compiled molecular concepts}

\item{preCalmatrix}{Precalculated matrix to penalize the redundancies in molecular concepts}

\item{selectn}{The maximum number of pathwasy for similarity matrix.}

\item{minsize}{Default 10. Filter out concepts that have genes less than minsize}

\item{rm.overfit}{Default FALSE.}
}
\value{
A similarity matrix between top N pathways
}
\description{
Title Calculation the functional associations of deambiguated pathway enrichment result
}
\examples{
pathwayAssociation(DeambiguatePathway=up.deambiguate, compare.list, feature.list, preCalmatrix, selectn=30, minsize=10)
}
