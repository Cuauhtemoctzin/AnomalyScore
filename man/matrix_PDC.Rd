% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PDCmatrix.R
\name{matrix_PDC}
\alias{matrix_PDC}
\title{Partial directed coherence matrix}
\usage{
matrix_PDC(unit, ar)
}
\arguments{
\item{unit}{A Matrix containing the multivariate time series. Each column
represents a univariate time series.}

\item{ar}{Integer vector containing all the lags considered for the
vector autoregressive model}
}
\value{
An real array of dimensions, ncol(unit), ncol(unit), n, where n is the number of
frequencies at which the PDC is estimated.
}
\description{
Partial directed coherence matrix
}
\examples{
X=matrix( rnorm(2000), ncol=10  )
ar=c(1, 2)
matrix_PDC(X, ar)
}
