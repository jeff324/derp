% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distributions.R
\name{logistic_lpdf}
\alias{logistic_lpdf}
\title{Log probability density function for the logistic distribution}
\usage{
logistic_lpdf(x, location, scale)
}
\description{
Log probability density function for the logistic distribution
}
\details{
Equivalent to \code{sum(dlogis(x,location,scale,log=TRUE))}. For more information see \code{\link{dlogis}}.
}
