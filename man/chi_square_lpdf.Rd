% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distributions.R
\name{chi_square_lpdf}
\alias{chi_square_lpdf}
\title{Log probability density function for the chi square distribution}
\usage{
chi_square_lpdf(x, df)
}
\description{
Log probability density function for the chi square distribution
}
\details{
Equivalent to \code{sum(dchisq(x,df,log=TRUE))}. For more information see \code{\link{dchisq}}.
}
