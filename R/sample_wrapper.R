#' Posterior sampling using differential evolution
#'
#' \code{powder} Runs power posterior sampling using differential evolution markov chain monte carlo
#' @param model A \code{function}. Defines a Bayesian model. See example below.
#' @param pars A \code{list}. Specifies the names of parameters along with their initialization bounds and blocks. See example below.
#' @param data A \code{list}. Specifies the data to be modeled. See example.
#' @param sampler A \code{character}. Specifies the name of the sampler to be used. \code{de} is the only built-in sampler.
#' @param num_samples A \code{numeric}. The total number of samples to collect.
#' @param num_chains A \code{numeric}. The number of chains to run.
#' @param migration_freq A \code{numeric}. Number of iterations to wait between each migration step. If \code{NULL}, migration is not done.
#' @param migration.start A \code{numeric}. Iteration to start migrating. This should be after chains are burned in. If \code{NULL}, migration is not done.
#' @param migration.end A \code{numeric}. Iteration to stop migrating. Migration should stop well before sampling is finished. If \code{NULL}, migration is not done.
#' @param randomize_phi A \code{logical}. Should the correlational structure between level-1 and level-2 parameters be ignored.
#' @param init_theta A list where each element contains a named vector of parameter initial start values for each subject.
#' @param init_phi A named vector parameter initial start values
#' @param parallel_backend A character vector either 'MPI', 'doParallel', or 'none' indicating backend for parallelization. Default is none.
#' @param n_cores A \code{numeric}. Number of cores when parallel_backend is specified.
#' @param benchmark A \code{logical}. Produces compute times benchmark purposes.
#' Oftentimes, ignoring this correlation will lead to better sampling.
#' @param update A \code{numeric}. Specifies the number of iterations before printing the current iteration number to the console.
#' @examples
#' \dontrun{
#' library(derp)
#' #### Simulate data
#' data = list(response = rnorm(100,2,1))
#' #### Model parameters
#' pars = list(
#' 'mu' =    list('init'=c(1,10)),
#' 'sd' =    list('init'=c(.1,5))
#' )
#' #### Define model
#' model = function() {
#'     normal_lpdf(mu, 0, 3)
#'     gamma_lpdf(sd, 1, 1)
#'     normal_lpdf(response, mu, sd)
#' }
#' #### Run the sampler
#' samples = run_mcmc(model, pars, data,
#'                   migration_start = 500, migration_end = 700,
#'                   migration_freq = 10, num_samples=3000)
#' mcmc_snip = snip(samples,burnin=1000,thin=1)
#' summary(mcmc_snip)
#' }
#' @export
run_mcmc = function(model, pars=NULL, data,
                    sampler='de', num_samples=NULL, num_chains=NULL, migrate=FALSE,
                    migration_start=NULL, migration_end=NULL,
                    migration_freq=NULL, randomize_phi=TRUE, update=100,
                    init_theta=NULL, init_phi=NULL, return_as_mcmc = TRUE,
                    parallel_backend='none',n_cores=NULL,benchmark=FALSE)
{

     if (!is.null(pars))
     {
          model = Model(model,pars)
     }

     if (is.null(num_samples))
     {
          warning('Number of samples not specified. Setting num_samples to 1000.',
                  call. = FALSE,immediate. = FALSE)
          num_samples = 1000
     }

     if (is.null(num_chains) & class(model)[1] == 'Model.Hierarchical')
     {
          n_pars = length(model$theta)
          n_hpars = max(sapply(model$blocks,function(x) length(x)))
          num_chains = 2*max(n_pars,n_hpars)

          if (num_chains < 4)
          {
               num_chains = 4
          }

     }

     if (is.null(num_chains) & class(model)[1] == 'Model.Individual')
     {
          n_pars = length(model$theta)
          num_chains = 2*n_pars

          if (num_chains < 4)
          {
               num_chains = 4
          }
     }



     if (is.null(migration_start))
     {
          migration_start = Inf
     }

     if (is.null(migration_end))
     {
          migration_end = -Inf
     }

     if (is.null(migration_freq))
     {
          migration_freq = 0
     }

     if(class(model)[1] == 'Model.Hierarchical')
     {

          if (is.character(sampler) & length(sampler)==1)
          {
               sampler_list = lapply(1:(length(model$blocks)+1), function(x) list(sampler=sampler,block=x-1,iter=c(1,num_samples)))
               sampler_funs = list(package.match.fun(sampler))
               names(sampler_funs) = sampler
          }


          if (is.character(sampler) & length(sampler) > 1)
          {
               sampler_list = lapply(1:length(sampler), function(x) list(sampler=sampler,block=x-1,iter=c(1,num_samples)))
               sampler_funs = lapply(unique(sampler), function(x) package.match.fun(x))
               names(sampler_funs) = unique(sampler)
          }


          if (is.list(sampler))
          {
               sampler_list = sampler
               sampler_names = sapply(sampler,function(x) x$sampler)
               sampler_funs = lapply(unique(sampler_names), function(x) package.match.fun(x))
               names(sampler_funs) = unique(sampler_names)
          }

          sampler_matrix = get_sampler_matrix(sampler_list)

          out = de.sample(model, data, sampler_funs, sampler_matrix, num_samples,
                          num_chains, migrate, migration_start, migration_end,
                          migration_freq, randomize_phi, update, init_theta, init_phi, return_as_mcmc,
                          parallel_backend,n_cores,benchmark)

     } else {

          if (is.character(sampler) & length(sampler) == 1)
          {
               sampler_list = list(list(sampler=sampler,iter=c(1,num_samples)))
               sampler_funs = list(package.match.fun(sampler))
               names(sampler_funs) = sampler
          }


          if (is.list(sampler))
          {
               sampler_list = sampler
               sampler_names = sapply(sampler,function(x) x$sampler)
               sampler_funs = lapply(unique(sampler_names), function(x) package.match.fun(x))
               names(sampler_funs) = unique(sampler_names)
          }

          sampler_vector = get_sampler_vector(sampler_list)

          out = de.sample_ind(model, data, sampler_funs, sampler_vector, num_samples,
                              num_chains, migration_start, migration_end,
                              migration_freq, randomize_phi, update)
     }
     return(out)
}

package.match.fun = function(x)
{
     if (x == 'de')
     {
          return(de)
     } else {
          return(match.fun(x))
     }
}

get_sampler_matrix = function(sampler_list)
{
     block_idx = lapply(sampler_list,function(x) x$block)
     iter_idx = lapply(sampler_list,function(x) x$iter)
     max_block = max(unlist(block_idx)+1)
     max_iter = max(unlist(iter_idx))
     mat = matrix('',max_iter,max_block)
     for (i in 1:length(sampler_list)) {
          block = sampler_list[[i]]$block
          iter = sampler_list[[i]]$iter
          sampler = sampler_list[[i]]$sampler
          mat[min(iter):max(iter),min(block+1):max(block+1)] = sampler
     }
     return(mat)
}

get_sampler_vector = function(sampler_list)
{
     iter_idx = lapply(sampler_list,function(x) x$iter)
     max_iter = max(unlist(iter_idx))
     mat = NULL
     for (i in 1:length(sampler_list))
     {
          iter = sampler_list[[i]]$iter
          sampler = sampler_list[[i]]$sampler
          mat[min(iter):max(iter)] = sampler
     }
     return(mat)
}
