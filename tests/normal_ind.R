rm(list=ls())
library(derp)

#### Simulate data
data = list(response = rnorm(100,2,1))

#### Model parameters
pars = list(
'mu' =    list('init'=c(1,10)),
'sd' =    list('init'=c(.1,5))
)

#### Define model
model = function() {
     normal_lpdf(mu, 0, 3)
     gamma_lpdf(sd, 1, 1)
     normal_lpdf(response, mu, sd)
}

#### Run the sampler
samples = run_mcmc(model, pars, data,
                   migration_start = 500, migration_end = 700,
                   migration_freq = 10, num_samples=3000)

mcmc_snip = snip(samples,burnin=1000,thin=1)
summary(mcmc_snip)


