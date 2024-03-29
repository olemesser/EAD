% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/complexityMeasures.R
\name{measure_PMH}
\alias{measure_PMH}
\title{Calculate the Product Mix Heterogeneity (PMH)}
\usage{
measure_PMH(P)
}
\arguments{
\item{P}{A product matrix (e.g, RES_CONS_PAT) with products in rows and features, components, process or resources in columns.}
}
\value{
Returns the corresponding PMH value
}
\description{
The product mix heterogeneity (PMH)as defined by \insertCite{Anderson.1995;textual}{EAD} for a given product matrix is calculated with this function.
Since the PMH is calculated for each individual resource the mean and standard deviation are returned.
}
\examples{

data('csd_EAD')
P<-CSD_EAD$RES_CONS_PAT

measure_PMH(P)
}
\references{
\insertAllCited{}
}
