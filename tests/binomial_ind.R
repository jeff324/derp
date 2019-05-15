rm(list=ls())
library(derp)


#### Simulate data
size = 50
data = list(response=rbinom(100,size=size,prob = .8),size=size)

#### Model parameters
pars = list(
'phi' =   list('init'=c(-1,1))
)


#### Define model
model = function() {

     theta = pnorm(phi)

     normal_lpdf(phi,0,3)
     binomial_lpdf(response, size, theta)

}

#### Define a new sampler
mh = function(current_chain, current_iter, num_chains, pars, par_names)
{
     x = pars[current_chain,] + rnorm(1,0,.1)
}


sampler = list(
list(sampler='de',iter=c(1,2000)),
list(sampler='mh',iter=c(2001,3000))
)


#### Run the sampler
samples = run_mcmc(model, pars, data, sampler=sampler,
                   migration_start = 500, migration_end = 700,
                   migration_freq = 10, num_samples=3000)

# run_mcmc returns an mcmc.list object from the coda package
# this makes it easy to run dignostics, summary statistics, and plots using coda
mcmc_snip = snip(samples,burnin=1000,thin=2)
summary(mcmc_snip)

# it is also easy to convert to matrix
mcmc_mat = as.matrix(mcmc_snip,iters=TRUE,chains=TRUE)

# or to a data frame
mcmc_df = as.data.frame(mcmc_mat)

coda::traceplot(mcmc_snip)
coda::densplot(mcmc_snip)
