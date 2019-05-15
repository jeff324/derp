#' Log probability density function for the normal distribution
#' @details
#' Equivalent to \code{sum(dnorm(x,mean,sd,log=TRUE))}. For more information see \code{\link{dnorm}}.
#' @export
normal_lpdf = function(x,mean,sd)
{
     if (!is.null(attr(x,'eval')))
     {
         lp = sum(dnorm(x,mean,sd,log=TRUE))
         new_lp = get('lp__',envir = parent.env(parent.frame())) + lp
         assign('lp__',new_lp,envir = parent.env(parent.frame()))
     }
}

#' Log probability density function for the gamma distribution
#' @details
#' Equivalent to \code{sum(dgamma(x,shape,scale,log=TRUE)))}. For more information see \code{\link{dgamma}}.
#' @export
gamma_lpdf = function(x,shape,scale)
{
     if (!is.null(attr(x,'eval')))
     {
          lp = sum(dgamma(x,shape=shape, scale=scale, log=TRUE))
          new_lp = get('lp__',envir = parent.env(parent.frame())) + lp
          assign('lp__',new_lp,envir = parent.env(parent.frame()))
     }
}

#' Log probability density function for the truncated normal distribution
#' @details
#' Equivalent to \code{sum(msm::dtnorm(x,mean,sd,lower,upper,log=TRUE))}. For more information see \code{\link[msm]{dtnorm}}.
#' @export
tnormal_lpdf = function(x,mean,sd,lower,upper)
{
     if (!is.null(attr(x,'eval')))
     {
          lp = sum(msm::dtnorm(x=x,mean=mean,sd=sd,lower=lower,upper=upper,log=TRUE))
          new_lp = get('lp__',envir = parent.env(parent.frame())) + lp
          assign('lp__',new_lp,envir = parent.env(parent.frame()))
     }
}

#' Log probability density function for the beta distribution
#' @details
#' Equivalent to \code{sum(dbeta(x,shape1,shape2,log=TRUE))}. For more information see \code{\link{dbeta}}.
#' @export
beta_lpdf = function(x,shape1,shape2)
{
     if (!is.null(attr(x,'eval')))
     {
         lp = sum(dbeta(x=x,shape1=shape1,shape2=shape2,log=TRUE))
         new_lp = get('lp__',envir = parent.env(parent.frame())) + lp
         assign('lp__',new_lp,envir = parent.env(parent.frame()))
     }
}

#' Log probability density function for the beta distribution
#' @details
#' Equivalent to \code{sum(dbinomial(x,size,prob,log=TRUE))}. For more information see \code{\link{dbinom}}.
#' @export
binomial_lpdf = function(x,size,prob)
{
     if (!is.null(attr(x,'eval')))
     {
          lp = sum(dbinom(x=x,size=size,prob=prob,log=TRUE))
          new_lp = get('lp__',envir = parent.env(parent.frame())) + lp
          assign('lp__',new_lp,envir = parent.env(parent.frame()))
     }
}

#' Log probability density function for the cauchy distribution
#' @details
#' Equivalent to \code{sum(dcauchy(x,location,scale,log=TRUE))}. For more information see \code{\link{dcauchy}}.
#' @export
cauchy_lpdf = function(x,location,scale)
{
     if (!is.null(attr(x,'eval')))
     {
          lp = sum(dcauchy(x=x,location=location,scale=scale,log=TRUE))
          new_lp = get('lp__',envir = parent.env(parent.frame())) + lp
          assign('lp__',new_lp,envir = parent.env(parent.frame()))
     }
}

#' Log probability density function for the exponential distribution
#' @details
#' Equivalent to \code{sum(dexp(x,rate,log=TRUE))}. For more information see \code{\link{dexp}}.
#' @export
exponential_lpdf = function(x,rate)
{
     if (!is.null(attr(x,'eval')))
     {
         lp = sum(dexp(x=x,rate=rate,log=TRUE))
         new_lp = get('lp__',envir = parent.env(parent.frame())) + lp
         assign('lp__',new_lp,envir = parent.env(parent.frame()))
     }
}

#' Log probability density function for the cauchy distribution
#' @details
#' Equivalent to \code{sum(dgeom(x,prob,log=TRUE))}. For more information see \code{\link{dgeom}}.
#' @export
geometric_lpdf = function(x,prob)
{
     if (!is.null(attr(x,'eval')))
     {
          lp = sum(dgeom(x=x,prob=prob,log=TRUE))
          new_lp = get('lp__',envir = parent.env(parent.frame())) + lp
          assign('lp__',new_lp,envir = parent.env(parent.frame()))
     }
}

#' Log probability density function for the uniform distribution
#' @details
#' Equivalent to \code{sum(dunif(x,min,max,log=TRUE))}. For more information see \code{\link{dunif}}.
#' @export
uniform_lpdf = function(x,min,max)
{
     if (!is.null(attr(x,'eval')))
     {
          lp = sum(dunif(x=x,min=min,max=max))
          new_lp = get('lp__',envir = parent.env(parent.frame())) + lp
          assign('lp__',new_lp,envir = parent.env(parent.frame()))
     }
}

#' Log probability density function for the logistic distribution
#' @details
#' Equivalent to \code{sum(dlogis(x,location,scale,log=TRUE))}. For more information see \code{\link{dlogis}}.
#' @export
logistic_lpdf = function(x,location,scale)
{
     if (!is.null(attr(x,'eval')))
     {
          lp = sum(dlogis(x=x,location=location,scale=scale,log=TRUE))
          new_lp = get('lp__',envir = parent.env(parent.frame())) + lp
          assign('lp__',new_lp,envir = parent.env(parent.frame()))
     }
}

#' Log probability density function for the poisson distribution
#' @details
#' Equivalent to \code{sum(dpois(x,lambda,log=TRUE))}. For more information see \code{\link{dpois}}.
#' @export
poisson_lpdf = function(x,lambda)
{
     if (!is.null(attr(x,'eval')))
     {
          lp = sum(dpois(x=x,lambda = lambda,log=TRUE))
          new_lp = get('lp__',envir = parent.env(parent.frame())) + lp
          assign('lp__',new_lp,envir = parent.env(parent.frame()))
     }
}

#' Log probability density function for the log normal distribution
#' @details
#' Equivalent to \code{sum(dlnorm(x,meanlog,sdlog,log=TRUE))}. For more information see \code{\link{dlnorm}}.
#' @export
log_normal_lpdf = function(x,meanlog,sdlog)
{
     if (!is.null(attr(x,'eval')))
     {
          lp = sum(dlnorm(x=x,meanlog=meanlog,sdlog=sdlog,log=TRUE))
          new_lp = get('lp__',envir = parent.env(parent.frame())) + lp
          assign('lp__',new_lp,envir = parent.env(parent.frame()))
     }
}

#' Log probability density function for the student t distribution
#' @details
#' Equivalent to \code{sum(dt(x,df,log=TRUE))}. For more information see \code{\link{dt}}.
#' @export
student_t_lpdf = function(x,df)
{
     if (!is.null(attr(x,'eval')))
     {
          lp = sum(dt(x=x,df=df,log=TRUE))
          new_lp = get('lp__',envir = parent.env(parent.frame())) + lp
          assign('lp__',new_lp,envir = parent.env(parent.frame()))
     }
}

#' Log probability density function for the negative binomial distribution
#' @details
#' Equivalent to \code{sum(dnbinom(x,size,prob,log=TRUE))}. For more information see \code{\link{dnbinom}}.
#' @export
negative_binomial_lpdf = function(x,size,prob)
{
     if (!is.null(attr(x,'eval')))
     {
          lp = sum(dnbinom(x=x,size=size,prob=prob,log=TRUE))
          new_lp = get('lp__',envir = parent.env(parent.frame())) + lp
          assign('lp__',new_lp,envir = parent.env(parent.frame()))
     }
}

#' Log probability density function for the chi square distribution
#' @details
#' Equivalent to \code{sum(dchisq(x,df,log=TRUE))}. For more information see \code{\link{dchisq}}.
#' @export
chi_square_lpdf = function(x,df)
{
     if (!is.null(attr(x,'eval')))
     {
          lp = sum(dchisq(x,df,log=TRUE))
          new_lp = get('lp__',envir = parent.env(parent.frame())) + lp
          assign('lp__',new_lp,envir = parent.env(parent.frame()))
     }
}

#' Creates custom distributions
#' @examples
#'\dontrun{
#' # define normal in terms of precision
#' f = function(x,mean,prec) dnorm(x, a, 1/sqrt(prec))
#' new_f = lpdf_custom(f)
#' }
#' @export
custom_lpdf = function(f)
{
     f_new = function(...)
     {
          args = list(...)
          if (!is.null(attr(args[[1]],'eval')))
          {
               lp = sum(log(do.call(f,args)))
               new_lp = get('lp__',envir = parent.env(parent.frame())) + lp
               assign('lp__',new_lp,envir = parent.env(parent.frame()))
          }
     }
     return(f_new)
}


