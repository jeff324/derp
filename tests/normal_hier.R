rm(list=ls())
library(derp)


#### Simulate data
data = list()
norm_subj = function() {
     sd_shape = rgamma(1,1,scale=1)
     sd_scale = rgamma(1,1,scale=1)
     mu_mu = rnorm(1,10,.1)
     mu_sd = rgamma(1,1,scale=1)
     mu = rnorm(1,mu_mu,mu_sd) #expected mean of 10, var of 1
     sd = rgamma(1,sd_shape,scale=sd_scale) #expected mean of 1, var of 1
     sim_dat = rnorm(200,mu,sd)
     data = list('response'=sim_dat)
}

for (i in 1:10) {
     data[[i]] = norm_subj()
}


#### Model parameters
pars = list(
'mu' =    list('init'=c(1,10),'block'=0),
'sd' =    list('init'=c(.1,5),'block'=0),
'mu_mu' = list('init'=c(5,15),'block'=1),
'mu_sd' = list('init'=c(.1,5),'block'=1),
'sd_shape' = list('init'=c(.1,5),'block'=2),
'sd_scale' = list('init'=c(.1,5),'block'=2)
)


#### Define model
model = function() {

     normal_lpdf(mu_mu, 0, 5)
     gamma_lpdf(mu_sd, 1, 1)
     gamma_lpdf(sd_shape, 1, 1)
     gamma_lpdf(sd_scale, 1, 1)

     normal_lpdf(mu, mu_mu, mu_sd)
     gamma_lpdf(sd, sd_shape, sd_scale)

     normal_lpdf(response, mu, sd)

}

#### Define a new sampler
mh = function(current_chain, current_iter, num_chains, pars, par_names)
{
     if (current_iter < 2000) {
          x = pars[current_chain,] + rnorm(1,0,1)
     } else {
          x = pars[current_chain,] + rnorm(1,0,.1)
     }
}

##### Define a complicated sampling scheme.
##### Note: This is likely to produce bad samples
##### and is just used an illustration
sampler = list(
list(sampler='de', block=0, iter=c(2,1000)),
list(sampler='mh', block=0, iter=c(1001,3000)),
list(sampler='de', block=c(1,2), iter=c(2,1000)),
list(sampler='mh', block=c(1,2), iter=c(1001,2000)),
list(sampler='de', block=1, iter=c(2001,3000)),
list(sampler='mh', block=2, iter=c(2001,3000))
)

#### Run the sampler
samples = run_mcmc(model, pars, data, sampler=sampler,
                   migration_start = 500, migration_end = 700,
                   migration_freq = 10, num_samples=3000)

# discard burnin and summarize
# Note the samples will probably not be good
# because of our poor sampling scheme
mcmc_snip = snip(samples,burnin=1000,thin=1)
summary(mcmc_snip)

coda::traceplot(mcmc_snip[,'mu_mu',])
mat = as.matrix(mcmc_snip)
hist(mat[,'mu_mu'],breaks = 50)
