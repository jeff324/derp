% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distributions.R
\name{student_t_lpdf}
\alias{student_t_lpdf}
\title{Log probability density function for the student t distribution}
\usage{
student_t_lpdf(x, df)
}
\description{
Log probability density function for the student t distribution
}
\details{
Equivalent to \code{sum(dt(x,df,log=TRUE))}. For more information see \code{\link{dt}}.
}
