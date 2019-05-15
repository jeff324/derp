library(derp)

lba_subj = function(n_sims)
{
     require(msm)
     require(rtdists)
     A_mu = rtnorm(1,1,.1,0,Inf)
     A_sd = rtnorm(1,.5,.1,0,Inf)
     b_mu = rtnorm(1,.5,.1,0,Inf)
     b_sd = rtnorm(1,.5,.1,0,Inf)
     t0_mu = rtnorm(1,.3,.1,0,Inf)
     t0_sd = rtnorm(1,.3,.1,0,Inf)
     v1_mu = rtnorm(1,1,.1,0,Inf)
     v1_sd = rtnorm(1,.5,.1,0,Inf)
     v2_c1_mu = rtnorm(1,2,.1,0,Inf)
     v2_c1_sd = rtnorm(1,.5,.1,0,Inf)
     v2_c2_mu = rtnorm(1,4,.1,0,Inf)
     v2_c2_sd = rtnorm(1,1,.1,0,Inf)
     sd1_mu = rtnorm(1,1,.1,0,Inf)
     sd1_sd = rtnorm(1,.5,.1,0,Inf)
     A = rtnorm(1,A_mu,A_sd,0,Inf)
     b = rtnorm(1,b_mu,b_sd,0,Inf)
     t0 = rtnorm(1,t0_mu,t0_sd,0,Inf)
     v1 = rtnorm(1,v1_mu,v1_sd,0,Inf)
     v2_c1 = rtnorm(1,v2_c1_mu,v2_c1_sd,0,Inf)
     v2_c2 = rtnorm(1,v2_c2_mu,v2_c2_sd,0,Inf)
     sd1 = rtnorm(1,sd1_mu,sd1_sd,0,Inf)
     data_cond_1 = rLBA(n=n_sims,A=A,b=A+b,t0=t0,distribution = 'norm',mean_v=c(v1,v2_c1),sd_v=c(sd1,1),silent=TRUE)
     data_cond_2 = rLBA(n=n_sims,A=A,b=A+b,t0=t0,distribution = 'norm',mean_v=c(v1,v2_c2),sd_v=c(sd1,1),silent=TRUE)
     cond = rep(1:2,each=n_sims)
     data = rbind(data_cond_1,data_cond_2)
     data = list(response=data$response,rt=data$rt,cond=cond)
}

data = lapply(1:10,function(x) lba_subj(100))



pars = list(
'A' = list('init'=c(0,2),'block'=0),
'b' = list('init'=c(0,2),'block'=0),
't0' = list('init'=c(0,1),'block'=0),
'v1' = list('init'=c(.1,3),'block'=0),
'v2_c1' = list('init'=c(.1,3),'block'=0),
'v2_c2' = list('init'=c(.1,3),'block'=0),
'sd1' = list('init'=c(.5,1.5),'block'=0),
'A_mu' = list('init'=c(0,2),'block'=1),
'b_mu' = list('init'=c(0,2),'block'=2),
't0_mu' = list('init'=c(0,1),'block'=3),
'v1_mu' = list('init'=c(.1,3),'block'=4),
'v2_c1_mu' = list('init'=c(.1,3),'block'=5),
'v2_c2_mu' = list('init'=c(.1,3),'block'=6),
'sd1_mu' = list('init'=c(.5,1.5),'block'=7),
'A_sd' = list('init'=c(0,2),'block'=1),
'b_sd' = list('init'=c(0,2),'block'=2),
't0_sd' = list('init'=c(0,1),'block'=3),
'v1_sd' = list('init'=c(.1,3),'block'=4),
'v2_c1_sd' = list('init'=c(.1,3),'block'=5),
'v2_c2_sd' = list('init'=c(.1,3),'block'=6),
'sd1_sd' = list('init'=c(.5,1.5),'block'=7)
)

lba = custom_lpdf(dLBA)

model = function()
{
     #concatenate for ease of use
     v2 = to_list(v2_c1,v2_c2)
     v2_mu = to_list(v2_c1_mu,v2_c2_mu)
     v2_sd = to_list(v2_c1_sd,v2_c2_sd)

     tnormal_lpdf(A_mu,1,3,0,Inf)
     tnormal_lpdf(A_sd,1,3,0,Inf)
     tnormal_lpdf(b_mu,.5,.5,0,Inf)
     tnormal_lpdf(b_sd,.5,.5,0,Inf)
     tnormal_lpdf(v1_mu,1,3,0,Inf)
     tnormal_lpdf(v1_sd,1,3,0,Inf)
     tnormal_lpdf(sd1_mu,1,1,0,Inf)
     tnormal_lpdf(sd1_sd,1,1,0,Inf)
     tnormal_lpdf(t0_mu,.3,.3,0,Inf)
     tnormal_lpdf(t0_sd,.3,.3,0,Inf)


     for (i in 1:2)
     {
          tnormal_lpdf(to_subset(v2_mu,i),3,2,0,Inf)
          tnormal_lpdf(to_subset(v2_sd,i),1,1,0,Inf)
          tnormal_lpdf(to_subset(v2,i),to_subset(v2_mu,i),to_subset(v2_sd,i),0,Inf)
     }

     tnormal_lpdf(A,A_mu,A_sd,0,Inf)
     tnormal_lpdf(b,b_mu,b_sd,0,Inf)
     tnormal_lpdf(v1,v1_mu,v1_sd,0,Inf)
     tnormal_lpdf(sd1,sd1_mu,sd1_sd,0,Inf)
     tnormal_lpdf(t0,t0_mu,t0_sd,0,Inf)

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

samples = run_mcmc(model,pars,data,num_samples = 1000,
                   migration_start=500,migration_end = 700,migration_freq = 10,update=10)

mcmc_snip = snip(samples,burnin=1,thin=1)
summary(mcmc_snip)

coda::traceplot(mcmc_snip)
