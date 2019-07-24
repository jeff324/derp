#' @import foreach

de = function(current_chain, current_iter, num_chains, pars, par_names)
{
     gamma = 2.38/sqrt(2*length(pars[current_chain,]))
     pop = c(1:num_chains)[-current_chain]
     index = sample(pop, size=2)
     x_update = pars[current_chain,] + gamma*(pars[index[1],]-pars[index[2],]) + runif(1,-.001,.001)
     return(x_update)
}

proposal = function(f,current_chain,current_iter,num_chains,pars,par_names)
{
     if (is.vector(pars))
     {
          pars = matrix(pars,num_chains,2)
          out = f(current_chain,current_iter,num_chains,pars,par_names)[1]
          names(out) = par_names
          return(out)
     } else {
          f(current_chain,current_iter,num_chains,pars,par_names)
     }
}

accept = function(cur_weight,weight)
{
     alpha = exp(weight-cur_weight)
     if (is.nan(alpha))
     {
          return(FALSE)
     } else if (is.na(alpha))
     {
          return(FALSE)
     } else if (runif(1) < alpha)
     {
          return(TRUE)
     } else
     {
          return(FALSE)
     }
}

migrate = function(x,cur_weight)
{

     new_weight = cur_weight
     n_chains = length(cur_weight)
     num_migration_chains = sample(1:n_chains,size=1)
     use_chains = sample(1:n_chains,size=num_migration_chains,replace=FALSE)

     for (i in 1:num_migration_chains) {

          #get weight of random chain
          migration_cur_weight = cur_weight[use_chains[i]]

          #go to the previous chain in use_chains
          new_chain = i - 1
          #if there is no previous chain, then go to the last chain
          if (new_chain == 0){
               new_chain = num_migration_chains
          }

          migration_weight = cur_weight[use_chains[new_chain]]

          #move the chain given acceptance probability
          if (accept(migration_cur_weight,migration_weight)) {
               if (is.vector(x))
               {
                    x[use_chains[i]] = x[use_chains[new_chain]]
               } else {
                    x[use_chains[i],] = x[use_chains[new_chain],]
               }
               new_weight[use_chains[i]] = new_weight[use_chains[new_chain]]

          }
     }

     return(list(x,new_weight))
}


de.sample_ind = function(model, data, sampler, sampler_vector,
                         num_samples, n_chains,
                         migrate_start, migrate_end,
                         migrate_freq, rand_phi, update)
{
     n_pars = length(model$theta)
     theta = array(NA,c(n_chains,n_pars,num_samples))
     weight_theta = array(-Inf,c(n_chains,num_samples))

     e_lp = new.env()
     log_prob = model$lp
     environment(log_prob) = e_lp

     colnames(theta) = model$theta
     par_names = model$theta

     cat('\n','Initializing parameters')

     for (k in 1:n_chains)
     {
          while (weight_theta[k,1] == -Inf)
          {
               for (p in 1:n_pars)
               {
                    theta[k,p,1] = model$initializer(model$theta[p])
               }
               x_theta = theta[k, ,1]
               list2env(c(set_eval_true(data),set_eval_true(x_theta),c('lp__'=0)),e_lp)
               log_prob()
               weight_theta[k,1] = e_lp$lp__ #likelihood + prior
          }
     }

     cat('\n','Sampling')
     for (i in 2:num_samples)
     {
          if (i %% update == 0)
          {
               cat('\n', i, '/', num_samples)
          }

          if ((i > migrate_start) & (i < migrate_end) & (i %% migrate_freq == 0))
          {
               m_out = migrate(theta[,,i-1],weight_theta[,i-1])
               theta[,,i] = m_out[[1]]
               weight_theta[,i] = m_out[[2]]
          } else {
               for (k in 1:n_chains)
               {
                    x_theta = proposal(sampler[[sampler_vector[i]]],k,i,n_chains,theta[,,i-1],par_names)
                    list2env(c(set_eval_true(data),set_eval_true(x_theta),c('lp__'=0)),e_lp)
                    log_prob()
                    weight = e_lp$lp__ #likelihood + prior
                    if (accept(weight_theta[k,i-1],weight))
                    {
                         theta[k,,i] = x_theta
                         weight_theta[k,i] = weight
                    } else {
                         theta[k,,i] = theta[k,,i-1]
                         weight_theta[k,i] = weight_theta[k,i-1]
                    }
               }
          }
     }

     samples = list('level_1'=theta)
     mcmc_list = utility.as_mcmc_list_ind(samples)
     return(mcmc_list)


}

set_eval_true = function(x){
     if (is.list(x))
     {
          rapply(x,function(x){ attr(x,'eval')=TRUE; x},how='replace')
     } else {
          lapply(x,function(x){ if (!is.na(x)) attr(x,'eval')=TRUE; x})
     }

}



de.sample = function(model, data, sampler, sampler_matrix,
                     num_samples, n_chains, migrate,
                     migrate_start, migrate_end,
                     migrate_step, rand_phi, update, init_theta, init_phi, return_as_mcmc,
                     parallel_backend,n_cores,benchmark)
{

     n_pars = length(model$theta)
     n_hpars = length(model$phi)
     n_subj = length(data)
     n_blocks = length(model$blocks)


     theta = array(NA,c(n_chains,n_pars,num_samples,n_subj))
     phi = array(NA,c(n_chains,n_hpars,num_samples))
     weight_theta = array(-Inf,c(n_chains,num_samples,n_subj))
     weight_phi = array(-Inf,c(n_chains,num_samples))

     e_lp = new.env()
     log_prob = model$lp
     environment(log_prob) = e_lp

     theta_names = model$theta
     phi_names = model$phi
     colnames(theta) = theta_names
     colnames(phi) = phi_names

     # initialize theta chains

     cat('\n','Initializing level-1 parameters')
     if (!is.null(init_theta))
     {
          cat('\n','Initializing with user-supplied values')
     }
     for (s in 1:n_subj) {
          cat('\n','Subject:',s,'/',n_subj)
          for (k in 1:n_chains) {
               while (weight_theta[k,1,s] == -Inf) {
                    for (p in 1:n_pars) {
                         if (is.null(init_theta))
                         {
                              theta[k,p,1,s] = model$initializer(model$theta[p])
                         } else {
                              init_try = as.numeric(init_theta[[s]][theta_names[p]]) + model$initializer(model$theta[p])
                              init_try = pmax(init_try,model$pars[[p]]$init[3])
                              init_try = pmin(init_try,model$pars[[p]]$init[4])
                              theta[k,p,1,s] = init_try
                         }

                    }
                    for (p in 1:n_hpars) {
                         if (is.null(init_phi))
                         {
                              phi[k,p,1] = model$initializer(model$phi[p])
                         } else {
                              phi[k,p,1] = as.numeric(init_phi[phi_names[p]]) + model$initializer(model$phi[p])
                         }
                    }
                    x_phi = phi[k,,1]
                    x_theta = theta[k, ,1,s]
                    list2env(c(set_eval_true(data[[s]]),x_theta,x_phi,c('lp__'=0)),e_lp)
                    log_prob()
                    weight_theta[k,1,s] = e_lp$lp__ #likelihood
               }
          }
     }


     #initialize phi chains
     cat('\n','Initializing level-2 parameters')
     if (!is.null(init_phi))
     {
          cat('\n','Initializing with user-supplied values')
     }
     cat('\n','Chain: ')
     for (k in 1:n_chains) {
          cat(k,' ')
          while (weight_phi[k,1] == -Inf) {
               for (p in 1:n_hpars) {
                    if (is.null(init_phi))
                    {
                         phi[k,p,1] = model$initializer(model$phi[p])
                    } else {
                         phi[k,p,1] = as.numeric(init_phi[phi_names[p]]) + model$initializer(model$phi[p])
                    }
               }
               x_phi = phi[k,,1]
               lp = 0
               for (s in 1:n_subj) {
                    x_theta = theta[k,,1,s]
                    list2env(c(data[[s]],set_eval_true(x_theta),x_phi,c('lp__'=0)),e_lp)
                    log_prob()
                    lp = e_lp$lp__ + lp #prior
               }
               list2env(c(data[[1]],x_theta,set_eval_true(x_phi)),e_lp)
               log_prob()
               weight_phi[k,1] = e_lp$lp__ + lp #hyperprior
          }
     }


     if (parallel_backend == 'doParallel')
     {
          if (is.null(n_cores))
          {
               n_cores = parallel::detectCores(all.tests = FALSE, logical = TRUE)
          }
          cat('\n','Running on', n_cores, 'cores')
          cl = parallel::makeCluster(n_cores)
          doParallel::registerDoParallel(cl)

     }

     # run DE-MCMC
     cat('\n','Sampling')
     chain_idx = 1:n_chains
     for (i in 2:num_samples) {
          if (i %% update == 0) {
               cat('\n', i, '/', num_samples)
          }

          #sample phi
          if (rand_phi) {
               chain_idx = sample(1:n_chains, size=n_chains, replace=FALSE)
          }

          if (benchmark)
          {
               start_time = Sys.time()
          }
          phi[,,i] = phi[,,i-1]
          for (p in 1:n_blocks) {
               par_range = model$blocks[[p]]
               #migration step
               if ((i > migrate_start) & (i < migrate_end) & (i %% migrate_step == 0) & migrate) {
                    if(p == 1) cat('\n','Migration Step')
                    phi_constant = phi[,,i]
                    #gets weights corresponding to parameter block, holding other parameters constant
                    weight_constant = NULL
                    for (k in 1:n_chains) {
                         #fix all parameters except current parameter block across all chains
                         all_pars = 1:length(phi[k,,i])
                         anti_par_block = all_pars[-par_range]
                         phi_constant[k,anti_par_block] = phi_constant[1,anti_par_block]
                         #update weight for each chain
                         x_phi = phi_constant[k,]
                         lp = 0
                         for (s in 1:n_subj) {
                              x_theta = theta[chain_idx[k],,i-1,s]
                              list2env(c(data[[1]],set_eval_true(x_theta),x_phi,c('lp__'=0)),e_lp)
                              log_prob()
                              lp = lp + e_lp$lp__ #prior
                         }
                         list2env(c(data[[1]],x_theta,set_eval_true(x_phi),c('lp__'=0)),e_lp)
                         log_prob()
                         lp = lp + e_lp$lp__ #prior + hyperprior
                         weight_constant = c(weight_constant, lp)
                    }
                    phi[,par_range,i] = migrate(phi[,par_range,i],weight_constant)[[1]]
                    #if we have updated all blocks then recompute weight for updated phi
                    if (p == n_blocks) {
                         for (k in 1:n_chains) {
                              #update weight for each chain
                              x_phi = phi[k,,i]
                              lp = 0
                              for (s in 1:n_subj) {
                                   x_theta = theta[chain_idx[k],,i-1,s]
                                   list2env(c(data[[1]],set_eval_true(x_theta),x_phi,c('lp__'=0)),e_lp)
                                   log_prob()
                                   lp = lp + e_lp$lp__ #prior
                              }
                              list2env(c(data[[1]],x_theta,set_eval_true(x_phi),c('lp__'=0)),e_lp)
                              log_prob()
                              lp = lp + e_lp$lp__
                              weight_phi[k,i] = lp #prior + hyperprior
                         }
                    }
               } else {
                    #crossover step

                    for (k in 1:n_chains) {
                         temp = phi[k,,i]
                         temp[par_range] = proposal(sampler[[sampler_matrix[i,p+1]]],k,i,n_chains,phi[,par_range,i],phi_names[par_range])
                         x_phi = temp

                         lp = 0
                         for (s in 1:n_subj)
                         {
                              x_theta = theta[chain_idx[k],,i-1,s]
                              list2env(c(data[[1]],set_eval_true(x_theta),x_phi,c('lp__'=0)),e_lp)
                              log_prob()
                              lp = lp + e_lp$lp__ #prior
                         }

                         list2env(c(data[[1]],x_theta,set_eval_true(x_phi),c('lp__'=0)),e_lp)
                         log_prob()
                         weight = lp + e_lp$lp__  #prior + hyperprior
                         if (accept(weight_phi[k,i-1],weight)) {
                              phi[k,par_range,i] = temp[par_range]
                              weight_phi[k,i] = weight
                         } else {
                              weight_phi[k,i] = weight_phi[k,i-1]
                         }
                    }
               }
          }

          if (benchmark)
          {
               end_time = Sys.time()
               total_time_l2 = difftime(end_time, start_time)
               cat('\n','Level-2 time',total_time_l2,units(total_time_l2))
          }

          #sample theta
          if (rand_phi) {
               chain_idx = sample(1:n_chains, size=n_chains, replace=FALSE)
          }

          if (benchmark)
          {
               start_time = Sys.time()
          }

          if (parallel_backend != 'none')
          {
               if ((i > migrate_start) & (i < migrate_end) & (i %% migrate_step == 0) & migrate) {
                    cat('\n','Migration Step')
                    for (s in 1:n_subj) {
                         m_out = migrate(theta[,,i-1,s],weight_theta[,i-1,s])
                         theta[,,i,s] = m_out[[1]]
                         weight_theta[,i,s] = m_out[[2]]
                    }
               } else {
                    theta_i = theta[,,(i-1):i,]
                    weight_theta_i = weight_theta[,(i-1):i,]
                    out = foreach(s = 1:n_subj) %dopar% {
                         for (k in 1:n_chains) {
                              temp = proposal(sampler[[sampler_matrix[i,1]]],k,i,n_chains,theta_i[,,1,s],theta_names)
                              x_theta = temp
                              x_phi = phi[chain_idx[k],,i]
                              list2env(c(set_eval_true(data[[s]]),set_eval_true(x_theta),x_phi,c('lp__'=0)),e_lp)
                              log_prob()
                              weight = e_lp$lp__ #likelihood + prior
                              if (accept(weight_theta_i[k,1,s],weight)) {
                                   theta_i[k,,2,s] = temp
                                   weight_theta_i[k,2,s] = weight
                              } else {
                                   theta_i[k,,2,s] = theta_i[k,,1,s]
                                   weight_theta_i[k,2,s] = weight_theta_i[k,1,s]
                              }
                         }
                         list(theta_i[,,,s],weight_theta_i[,,s])
                    }
                    for (s in 1:n_subj)
                    {
                         theta[,,i,s] = out[[s]][[1]][,,2]
                         weight_theta[,i,s] = out[[s]][[2]][,2]
                    }
               }

          } else {

               for (s in 1:n_subj) {
                    if ((i > migrate_start) & (i < migrate_end) & (i %% migrate_step == 0) & migrate) {
                         m_out = migrate(theta[,,i-1,s],weight_theta[,i-1,s])
                         theta[,,i,s] = m_out[[1]]
                         weight_theta[,i,s] = m_out[[2]]
                    } else {
                         for (k in 1:n_chains) {
                              temp = proposal(sampler[[sampler_matrix[i,1]]],k,i,n_chains,theta[,,i-1,s],theta_names)
                              x_theta = temp
                              x_phi = phi[chain_idx[k],,i]
                              list2env(c(set_eval_true(data[[s]]),set_eval_true(x_theta),x_phi,c('lp__'=0)),e_lp)
                              log_prob()
                              weight = e_lp$lp__ #likelihood + prior
                              if (accept(weight_theta[k,i-1,s],weight)) {
                                   theta[k,,i,s] = temp
                                   weight_theta[k,i,s] = weight
                              } else {
                                   theta[k,,i,s] = theta[k,,i-1,s]
                                   weight_theta[k,i,s] = weight_theta[k,i-1,s]
                              }
                         }
                    }
               }
          }
          if (benchmark)
          {
               end_time = Sys.time()
               total_time_l1 = difftime(end_time, start_time)
               total_time = total_time_l1 + total_time_l2
               projected_time = total_time * (num_samples - i)
               cat('\n','Level-1 time',total_time_l1, units(total_time_l1))
               cat('\n','Total iteration time',total_time, units(total_time))
               if (units(projected_time) == 'secs' & as.numeric(projected_time,units='secs') > 60)
               {
                    #compute projected time with minutes
                    projected_time = as.numeric(projected_time,units='mins')
                    if (projected_time > 60)
                    {
                         cat('\n','Projected time', projected_time/60, 'hours')
                    } else {
                         cat('\n','Projected time', projected_time, 'mins')
                    }
               } else if (units(projected_time) == 'mins' & as.numeric(projected_time,units='mins') > 60)
               {
                    cat('\n','Projected time', projected_time/60, 'hours')
               } else {
                    cat('\n','Projected time', projected_time, units(projected_time))
               }
          }
     }

     if (parallel_backend == 'doParallel')
     {
          doParallel::stopImplicitCluster()
     }

     samples = list('level_1'=theta,'level_2'=phi)

     if (return_as_mcmc)
     {
          samples = utility.as_mcmc_list(samples)
     }

     return(samples)

}




