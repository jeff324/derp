library(derp)
library(rtdists)

data_cond_1 = rLBA(n=300,A=1,b=1.5,t0=.3,distribution = 'norm',mean_v=c(1,2),sd_v=c(1,1),silent=TRUE)
data_cond_2 = rLBA(n=300,A=1,b=1.5,t0=.3,distribution = 'norm',mean_v=c(1,4),sd_v=c(1,1),silent=TRUE)
cond = rep(1:2,each=300)
data = rbind(data_cond_1,data_cond_2)
data = list(response=data$response,rt=data$rt,cond=cond)

pars = list(
'A' = list('init'=c(0,2)),
'b' = list('init'=c(0,2)),
't0' = list('init'=c(0,1)),
'v1' = list('init'=c(.1,3)),
'v2_c1' = list('init'=c(.1,3)),
'v2_c2' = list('init'=c(.1,3)),
'sd1' = list('init'=c(.5,1.5))
)

lba = custom_lpdf(dLBA)

model = function()
{
     #concatenate for ease of use
     v2 = to_list(v2_c1,v2_c2)

     for (i in 1:2)
     {
          tnormal_lpdf(to_subset(v2,i),2,2,0,Inf)
     }

     tnormal_lpdf(A,1,1,0,Inf)
     tnormal_lpdf(b,1,1,0,Inf)
     tnormal_lpdf(v1,1,1,0,Inf)
     tnormal_lpdf(sd1,.5,.5,0,Inf)
     tnormal_lpdf(t0,.3,.3,0,Inf)

     for (i in 1:2)
     {
          lba(rt=to_subset(rt,cond==i),
              response=to_subset(response,cond==i),
              A=c(A,A),
              b=c(A+b,A+b),
              t0=c(t0,t0),
              distribution='norm',
              mean_v=c(v1,to_subset(v2,i)),
              sd_v=c(sd1,1),
              silent=TRUE)
     }
}

samples = run_mcmc(model,pars,data,num_samples = 2000,
                   migration_start=500,migration_end = 700,migration_freq = 10)

mcmc_snip = snip(samples,burnin=1000,thin=1)
summary(mcmc_snip)

coda::traceplot(mcmc_snip)
