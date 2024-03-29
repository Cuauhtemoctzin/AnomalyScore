% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GrangerCausalityIndeces.R
\name{mBTSRGPDC}
\alias{mBTSRGPDC}
\title{Restricted Generalized Partial Directed Coherence}
\usage{
mBTSRGPDC(xM, pmax, freqs)
}
\arguments{
\item{xM}{the matrix of K time series (variables in columns)}

\item{pmax}{maximum order(lag) of the VAR model to be considered}

\item{freqs}{frequencies at which the spectral density is estimated}

\item{responseindex}{the index of the response variable in {1,...,K}}
}
\value{
the matrix of all the Restricted Generalized Partial Directed
Coherence index across the series of a multivariate set.
#'
}
\description{
partial directed coherence matrix values based on the mBTS algorithm
for estimation of the VAR parameters
}
