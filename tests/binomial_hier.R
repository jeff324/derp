rm(list=ls())
library(derp)
#
# z = function()
# {
#      x = 1
#      y = get('y',envir = parent.env(parent.frame()))
#      print(paste0('current ',y))
#      print(ls(envir = parent.env(parent.frame())))
#      new_y = y+1
#      assign('y',new_y,envir= parent.env(parent.frame()))
#      print(paste0('new ',new_y))
# }
# f = function()
# {
#      z()
# }
# e = new.env()
# environment(f) = e
# assign('y',1,envir=e)
#
# f()


#### Simulate data
data = list()
binom_subj = function(size=50) {
     #expected mean of theta will be .2
     #expected mean of mu ~= -.84
     mu = rnorm(1,qnorm(.2),.1)
     sigma = rgamma(1,1,scale=1)
     phi = rnorm(1,mu,sigma)
     theta = pnorm(phi)
     sim_dat = rbinom(100,size=size,prob = theta)
     data = list('response'=sim_dat,'size'=size)
}

data = lapply(1:30,function(x) binom_subj())


#### Model parameters
pars = list(
'phi' =   list('init'=c(-1,1),'block'=0),
'mu'  =   list('init'=c(-1,1),'block'=1),
'sigma' = list('init'=c(.5,5),'block'=1)
)


#### Define model
model = function() {

     theta = pnorm(phi)

     normal_lpdf(mu,0,3)
     gamma_lpdf(sigma,1,1)
     normal_lpdf(phi,mu,sigma)
     binomial_lpdf(response, size, theta)

}

#### Run the sampler
samples = run_mcmc(model, pars, data, sampler='de',
                   migration_start = 500, migration_end = 700,
                   migration_freq = 10, num_samples=3000)

# run_mcmc returns an mcmc.list object from the coda package
# this makes it easy to run dignostics, summary statistics, and plots using coda
mcmc_snip = snip(samples,burnin=1000,thin=1)
summary(mcmc_snip)

coda::traceplot(mcmc_snip[,'mu',])

# it is also easy to convert to matrix
mcmc_mat = as.matrix(mcmc_snip,iters=TRUE,chains=TRUE)

# or to a data frame
mcmc_df = as.data.frame(mcmc_mat)

