% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{sample_data}
\alias{sample_data}
\title{Generate a random sample of an VAR(p) m-variate time series.}
\usage{
sample_data(Phi, Sigma, p, m, n, NumLags, ...)
}
\arguments{
\item{Phi}{Transition Matrix}

\item{Sigma}{Covariance of error terms}

\item{p}{Order of VAR(p)}

\item{m}{Number of dimensions}

\item{n}{Number of observations}

\item{NumLags}{Number of pre-data obs set to zero.}
}
\value{
mx(n+NumLags) matrix
}
\description{
Generate a random sample of an VAR(p) m-variate time series.
}

