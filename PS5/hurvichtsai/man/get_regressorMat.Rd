% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{get_regressorMat}
\alias{get_regressorMat}
\title{Construct stacked regressor matrix of lags of order p}
\usage{
get_regressorMat(Z, n0, n, p)
}
\arguments{
\item{Z}{Original time series}

\item{n0}{Initial observation}

\item{n}{Number of observations}

\item{p}{Order of VAR(p)}
}
\description{
Construct stacked regressor matrix of lags of order p
}

