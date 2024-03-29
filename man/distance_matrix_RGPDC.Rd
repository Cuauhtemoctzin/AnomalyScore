% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Distance_matrix_knn_algorithm.R
\name{distance_matrix_RGPDC}
\alias{distance_matrix_RGPDC}
\title{Pairwise distance matrix based on the restricted generalized partial
directed coherence}
\usage{
distance_matrix_RGPDC(unit, pmax, period)
}
\arguments{
\item{unit}{A matrix representing a multivariate time series where each
column is a univariate time series.}

\item{pmax}{maximum order(lag) of the VAR model to be considered}

\item{period}{Integer referencing the index of the frequency to use for the
distance. It gives the Hertz or periods per unit of time; i.e., if the
sampling is per minute, and each hour cycle is the period of interest}
}
\value{
a matrix with pairwise distances
}
\description{
Pairwise distance matrix of a multivariate time series based on the
the restricted generalized partial directed coherence distance
between two series
}
\examples{
X=matrix( rnorm(2000), ncol=10  )
pmax=4
period=3
distance_matrix_RGPDC(unit=X, pmax, period)
}
\seealso{
Siggiridou, Elsa, Vasilios K. Kimiskidis, and Dimitris Kugiumtzis.
“Dimension Reduction of Frequency-Based Direct Granger Causality Measures
on Short Time Series.” Journal of Neuroscience Methods 289 (September 2017)
: 64–74. https://doi.org/10.1016/j.jneumeth.2017.06.021.
}
