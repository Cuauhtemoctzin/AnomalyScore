% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GrangerCausalityIndeces.R
\name{mBTS}
\alias{mBTS}
\title{modified Back-in-time Selection for vector AR parameters estimation}
\usage{
mBTS(xM, responseindex, pmax)
}
\arguments{
\item{xM}{the matrix of K time series (variables in columns)}

\item{responseindex}{the index of the response variable in \eqn{\{1,\ldots,K\}}}

\item{pmax}{maximum order(lag) of the VAR model to be considered}
}
\value{
the matrix of all explanatory lagged variables in the
DR model. The sequence of the lagged variables in 'lagM'
}
\description{
modified Back-in-time Selection for vector AR parameters estimation
}
\seealso{
I. Vlachos and D. Kugiumtzis, "Backward-in-time selection of the
order of dynamic regression prediction model," J. Forecast.,
vol. 32, pp. 685-701, 2013.
}
