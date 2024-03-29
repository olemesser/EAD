% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genRC.R
\name{genRC}
\alias{genRC}
\title{Generates resource costs}
\usage{
genRC(n_RES, DISP1 = 10, DISP2, TC)
}
\arguments{
\item{n_RES}{Number of resources.}

\item{DISP1}{The number of top resources which should allocate at DISP2 percent of total costs.
Therefore 1<=DISP1<NUMB_RES.}

\item{DISP2}{How much percentage of total costs are assigned to the top DISP1 resources? DISP2 0<DISP2<1 \cr}

\item{TC}{A integer of total costs.}
}
\description{
The resource costs are generated by Anand et al. (2017) procedure.
}
\examples{
# Anand's standard setting
DISP1<-10
DISP2<-0.30
TC<-10^6

RC<-genRC(n_RES=50,DISP1,DISP2,TC)
}
