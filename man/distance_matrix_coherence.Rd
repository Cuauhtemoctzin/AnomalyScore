% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Distance_matrix_knn_algorithm.R
\name{distance_matrix_coherence}
\alias{distance_matrix_coherence}
\title{Distance matrix from a coherence measure}
\usage{
distance_matrix_coherence(unit, span1, span2, period)
}
\arguments{
\item{unit}{A matrix representing a multivariate time series where each
column is a univariate time series.}

\item{span1}{Odd integer giving the widths of modified Daniell
smoothers to be used to smooth the periodogram. Refers to the bandwidth of
the smoothing process.}

\item{span2}{Odd integer giving the widths of modified Daniell
smoothers to be used to smooth the periodogram. Control another level of
smoothing to the spectral density estimation without altering the peaks}

\item{period}{Integer referencing the index of the frequency to use for the
distance. It gives the Hertz or periods per unit of time; i.e., if the
sampling is per minute, and each hour cycle is the period of interest}
}
\value{
a matrix with pairwise distances
}
\description{
Pairwise distance matrix of a multivariate time series based on computing the
squared coherence and transformed it to represent a distance at a specific
frequency
}
\examples{
X=matrix( rnorm(2000), ncol=10  )
span1=2
span2=2
period=3
distance_matrix_coherence(unit=X, span1, span2, period )
}
\seealso{
\code{\link[astsa:mvspec]{astsa::mvspec()}}; https://github.com/nickpoison/astsa/.
}
