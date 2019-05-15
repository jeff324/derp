utility.as_mcmc_list_ind = function(samples)
{
     theta = samples$level_1
     n_chains = dim(theta)[1]
     n_pars = dim(theta)[2]
     iter = dim(theta)[3]
     pars = colnames(theta)
     mat = lapply(1:n_chains,function(x) matrix(NA,iter,n_pars))
     for (k in 1:n_chains)
     {
          idx = 1
          column_names = NULL

          for (i in 1:n_pars)
          {
               mat[[k]][,idx] = c(theta[k,i,])
               column_names = c(column_names,pars[i])
               idx = idx + 1
          }

          colnames(mat[[k]]) = column_names
     }

     mcmc_list = lapply(mat,function(x) coda::as.mcmc(x))
     mcmc_list = coda::as.mcmc.list(mcmc_list)

     return(mcmc_list)
}

utility.as_mcmc_list = function(samples)
{
     theta = samples$level_1
     n_chains = dim(theta)[1]
     n_pars = dim(theta)[2]
     iter = dim(theta)[3]
     n_subj = dim(theta)[4]
     pars = colnames(theta)
     mat = lapply(1:n_chains,function(x) matrix(NA,iter,n_pars*n_subj))
     for (k in 1:n_chains)
     {
          idx = 1
          column_names = NULL
          for (s in 1:n_subj)
          {
               for (i in 1:n_pars)
               {
                    mat[[k]][,idx] = c(theta[k,i,,s])
                    column_names = c(column_names,paste0(pars[i],'.',s))
                    idx = idx + 1
               }
          }
          colnames(mat[[k]]) = column_names
     }

     phi = samples$level_2
     n_chains = dim(phi)[1]
     n_pars = dim(phi)[2]
     iter = dim(phi)[3]
     mat_phi = lapply(1:n_chains, function(x) matrix(NA,iter,n_pars))
     pars = colnames(phi)
     for (k in 1:n_chains)
     {
          column_names = NULL
          for (i in 1:n_pars)
          {
               mat_phi[[k]][,i] = c(phi[k,i,])
               column_names = c(column_names,pars[i])
          }
          colnames(mat_phi[[k]]) = column_names
     }

     mat_all = lapply(1:n_chains,function(x) cbind(mat[[x]],mat_phi[[x]]))
     mcmc_list = lapply(mat_all,function(x) coda::as.mcmc(x))
     mcmc_list = coda::as.mcmc.list(mcmc_list)

     return(mcmc_list)
}

#' Discard burn-in and/or thin samples
#'
#' @param mcmc_list A \code{mcmc.list} object.
#' @param burnin A \code{numeric}. All samples collected below this iteration will be discarded.
#' @param thin A \code{numeric}. The thinning interval between consecutive observations.
#'
#' @return
#' A \code{mcmc.list} object
#' @export
snip = function(mcmc_list,burnin,thin=1)
{
     end = end(mcmc_list)
     thin_seq = seq(burnin,end,by=thin)
     mcmc_list = coda::as.mcmc.list(lapply(mcmc_list,function(x) coda::mcmc(as.matrix(x[thin_seq,]))))
     return(mcmc_list)
}

#' @export
to_subset <- function(X, ...) {
     l <- X[...]
     if (is.list(l) & length(l) == 1)
     {
          l = l[[1]]
     }
     attr.names <- names(attributes(X))
     attr.names <- attr.names[attr.names != 'names']
     attributes(l)[attr.names] <- attributes(X)[attr.names]
     return(l)
}

#' @export
to_list = function(...)
{
     v = list(...)
     return(v)
}
