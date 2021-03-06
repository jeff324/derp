% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distributions.R
\name{tnormal_lpdf}
\alias{tnormal_lpdf}
\title{Log probability density function for the truncated normal distribution}
\usage{
tnormal_lpdf(x, mean, sd, lower, upper)
}
\description{
Log probability density function for the truncated normal distribution
}
\details{
Equivalent to \code{sum(msm::dtnorm(x,mean,sd,lower,upper,log=TRUE))}. For more information see \code{\link[msm]{dtnorm}}.
}
