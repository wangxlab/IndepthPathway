% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uniConSig.R
\name{CalWeight}
\alias{CalWeight}
\title{Title Calculate weight based on Ochiai Index}
\usage{
CalWeight(list1, list2, method = "Ochiai")
}
\arguments{
\item{list1}{Target gene}

\item{list2}{Genes symbols in Gene Signatures (MSigDB)}

\item{method}{"Ochiai" or "Jaccard"}
}
\value{
Calculated weight of target genes
}
\description{
Title Calculate weight based on Ochiai Index
}
\examples{
CalWeight(list1,list2,method="Ochiai")
}
