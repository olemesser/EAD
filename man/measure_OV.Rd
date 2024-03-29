% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/complexityMeasures.R
\name{measure_OV}
\alias{measure_OV}
\title{Calculates the option variability}
\usage{
measure_OV(P, norm = T)
}
\arguments{
\item{P}{A product matrix containing the features with products in rows and features in columns.}

\item{norm}{Boolean, default=T normalizes the HVM in a range of 0-1 where HVM=1 is the maximum possible value for a given matrix.}
}
\value{
Returns the OV value
}
\description{
This function calculates the option variability measure defined by \insertCite{MacDuffie1996;textual}{EAD}
}
\examples{

P<-matrix(c(1,0,0,1,
              0,1,0,0,
              1,0,0,0,
              0,1,1,0,
              0,1,0,1),
       nrow = 4,
       ncol = 5)
measure_OV(P)
}
\references{
\insertAllCited{}
}
