% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GrangerCausalityIndeces.R
\name{multilagmatrix}
\alias{multilagmatrix}
\title{multilagmatrix}
\usage{
multilagmatrix(xM, responseindex, ordersV, indexV)
}
\arguments{
\item{xM}{the matrix of K time series (variables in columns)}

\item{responseindex}{the index of the response variable in {1,...,K}}

\item{ordersV}{vector of size 1xK of the maximum order for each of the K
variables.}

\item{indexV}{the vector of size 1 x K\emph{pmax of zeros and ones
e.g. if the component in position 2}pmax+3 is one, the
third variable, lag 3, X3(t-3), is selected.}
}
\value{
the matrix of all explanatory lagged variables in the
DR model. The sequence of the lagged variables in 'lagM'
}
\description{
multilagmatrix builds the set of explanatory variables for the dynamic
regression model.
}
\seealso{
Kugiumtzis, D. “Direct-Coupling Information Measure from Nonuniform
Embedding.” Physical Review E 87, no. 6 (June 25, 2013): 062918.
https://doi.org/10.1103/PhysRevE.87.062918.
}
