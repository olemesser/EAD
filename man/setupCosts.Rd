% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cc_Driver.R
\name{setupCosts}
\alias{setupCosts}
\title{Calculates the setup costs}
\usage{
setupCosts(P_PD, DMM_PD_PrD, DSM_PrD, DEMAND, C_hold, C_order, C_setup)
}
\arguments{
\item{C_hold}{Component's holding costs.}

\item{C_order}{Component's order costs.}

\item{C_setup}{The average setup costs.}

\item{EAD}{An EAD object created by \link[EAD]{crt_EAD}.}
}
\value{
Returns the summed setup costs as well as the lot sizes for the individual components
}
\description{
This function calculates the setup costs for a specific EAD design. The exact procedure is described in \insertCite{Meerschmidt.2024;textual}{EAD}.
See details for more information.
}
\details{
For further details see \insertCite{Meerschmidt.2024;textual}{EAD} as well as the sub models of \insertCite{Thonemann.2000;textual}{EAD} and \insertCite{Zhang.2020;textual}{EAD}.
A vignette is available by running \code{utils::vignette('setupCosts',package='EAD')}.
}
\examples{

require(EAD)
set.seed(1234)
EAD <- smallEAD
EAD$DEMAND <- c(10,50,100,30)
C_hold <- rep(0.05,4)
C_order <- rep(0.15,4)
C_setup <- 444.4

clc_setupCosts(EAD,
               C_hold,
               C_order,
               C_setup)
}
\references{
\insertAllCited{}
}
