## Program: DALEC_DEMCMC.R
## Author: John W Smith Jr
## Date: Feb 19th, 2019
## Description: Runs an MCMC routine in order to estimate parameters in an ecosystem state space model, in 
## particular the DALEC model. Utilizes library Rcpp for fast likelihood evaluations for MH steps
## Sampling is done through DEMCMC method (Der Braak 2007)

## sources MCMC preparation (see file for details)
setwd('/home/john/Sampler')
source('prepare_MCMC.R')
#pupper[1] = .02
source('make_u.R')
#library(truncnorm)
g <- .2
l.g <- .2
#plower[9] <- 9e-6
## ----- Performing MCMC ------ ##
time1 <- Sys.time()
#library(profvis)
#profvis({
#chain_length <- 8000
eps <- sqrt(.Machine$double.eps)
for (i in 2:chain_length){
  
  for (gen in 1:num_gen){
    p[i,,gen] <- p[i-1,,gen]
    buildCpredi <- Cpred(p[i,,gen], C_samples[, , i-1,gen], LAI_ = pmax(.1, C_samples[,1,i-1,gen]/LMA_e), maxt, mint, Ca, lat_e, yearday, Nit_e, rad, 730)
    Cpredi[,,gen] <- cbind(buildCpredi$Cfpred, buildCpredi$Cwpred, buildCpredi$Crpred,
                    buildCpredi$Clitpred, buildCpredi$Csompred)
    Gpred[,gen] <- buildCpredi$G
    rm(buildCpredi)
    m <- sample(1:num_gen, 2)
    u <- make_u(p[i-1,,gen], plower, pupper)
    if (any(u > 1 - eps)){
      u[which(u > 1-eps)] <- u[which(u > 1-eps)]  - runif(length(which(u > 1-eps)), eps, 1e-6)
    }
    if (any(u < eps)){
      u[which(u < eps)] <- u[which(u < eps)] + runif(length(which(u < eps)), eps, 1e-6) 
    }
    l <- log(u / (1-u))
    
    um1 <- make_u(p[i-1,,m[1]], plower, pupper)
    um2 <- make_u(p[i-1,,m[2]], plower, pupper)
    l1 <- log(um1 / (1-um1))
    l2 <- log(um2 / (1-um2))
    l_star <- l + (g*rbinom(11, 1, .85))*(l1 - l2) + runif(11, 1e-4, 1e-4)
                                                             #.2*(c(.2, .05, .1, .25, .25, .2, .1, .3, .2, .15)))
    #l_star 
    #exp(l_star) / (1 + exp(l_star))
    #(exp(l_star) / (1 + exp(l_star)))*(pupper - plower) + plower
    #if (any(is.nan(l_star)) || any(is.infinite(l_star))){
    #  if (any(is.nan(l_star))){
    #    l_star[which(is.nan(l_star))] <- l[which(is.nan(l_star))]
    #  } 
    #  if (any(is.infinite(l_star))){
    #    l_star[which(is.infinite(l_star))] <- l[which(is.infinite(l_star))]
    #  }
    #}
    u_star <-  exp(l_star) / (1 + exp(l_star)) #pcauchy(l_star)
    p_star <- u_star*(pupper - plower) + plower
    #if (!(is.null(inf.ind.1))){
    #  p_star[inf.ind.1] <- plower[inf.ind.1] + runif(length(inf.ind.1), 0, 1e-6)
    #}
    #if (!(is.null(inf.ind.2))){
    #  p_star[inf.ind.2] <- pupper[inf.ind.2] - runif(length(inf.ind.2), 0, 1e-6)
    #}
    buildCpredstar <- Cpred(p_star, C_samples[, , i-1,gen], LAI_ = pmax(.1, C_samples[,1,i-1,gen]/LMA_e), maxt, mint, Ca, lat_e, yearday, Nit_e, rad, 730)
    Cpred_star <- cbind(buildCpredstar$Cfpred, buildCpredstar$Cwpred, buildCpredstar$Crpred,
                        buildCpredstar$Clitpred, buildCpredstar$Csompred)
    Gpred_star <- buildCpredstar$G
    rm(buildCpredstar)
    
    threshold <- cpp_ss_likelihood(C_obs, C_samples[, , i-1, gen], 1 / prec_obs, Cpred_star, 1 / prec_add_sample[i-1,,gen], 
                                   c(100,9200,100,20,11000), init_var, 730, G_pred = Gpred_star, G_e = G_obs, G = G[i-1,,gen], ind = missing_data_ind - 1) - 
      cpp_ss_likelihood(C_obs, C_samples[, , i-1,gen], 1 / prec_obs, Cpredi[,,gen], 1 / prec_add_sample[i-1,,gen], 
                        c(100,9200,100,20,11000), init_var, 730, G_pred <- Gpred[,gen], G_e = G_obs, G = G[i-1,,gen], ind = missing_data_ind - 1)
    if (log(runif(1,0,1)) < threshold){
      #acceptance2[i,,gen] <- 1
      acceptance2[i,gen] <- 1
      p[i,,gen] <- p_star
      Cpredi[,,gen] <- Cpred_star
      Gpred[,gen] <- Gpred_star
    }
    
    prec_add_sample[i,,gen] <- prec_add_sample[i-1,,gen]
    m.prec <- sample(1:num_gen, 2)
    l.prec <- log(prec_add_sample[i,,gen])
    prec_star <- exp(l.prec + l.g * rbinom(5, 1, .8) * (log(prec_add_sample[i-1,,m.prec[1]]) - log(prec_add_sample[i-1,,m.prec[2]])) + rnorm(5,0,.005) ) 
    threshold <- cpp_ss_likelihood(C_obs, C_samples[, , i-1,gen], 1 / prec_obs, Cpredi[,,gen], 1 / prec_star, 
                                   c(100,9200,100,20,11000), init_var, 730, G_pred = Gpred[,gen], G_e = G_obs, G=G[i-1,,gen], ind = missing_data_ind - 1) - 
      cpp_ss_likelihood(C_obs, C_samples[, , i-1,gen], 1 / prec_obs, Cpredi[,,gen], 1 / prec_add_sample[i-1,,gen], 
                        c(100,9200,100,20,11000), init_var, 730, G_pred = Gpred[,gen], G_e = G_obs, G=G[i-1,,gen], ind = missing_data_ind - 1)
    if (log(runif(1,0,1)) < threshold){
      prec_add_sample[i,,gen] <- prec_star
      #prec_Accept2[i,,gen] <- 1
      prec_Accept2[i,gen] <- 1
    }
    
    C.0[i,1,gen] <- get_latent_sample(atp = 1 - p[i,5,gen], btp = G[i-1,1,gen]*(1-p[i,2,gen])*p[i,3,gen], initial_prec = prec_init[1], add_prec = prec_add_sample[i-1,1,gen],
                                  xtp1 = C_samples[1,1,i-1,gen], data = C_obs[1,1], type = 'first', init_mean = init_mean[1])
    
    C.0[i,2,gen] <- get_latent_sample(atp = 1 - p[i,6,gen], btp = G[i-1,1,gen]*(1-p[i,2,gen])*p[i,3,gen], initial_prec = prec_init[2], add_prec = prec_add_sample[i-1,2,gen],
                                  xtp1 = C_samples[1,2,i-1,gen], data = C_obs[1,2], type = 'first', init_mean = init_mean[2])
    
    C.0[i,3,gen] <- get_latent_sample(atp = 1 - p[i,7,gen], btp = G[i-1,1,gen]*(1-p[i,2,gen])*(1-p[i,3,gen])*p[i,4,gen],
                                  initial_prec = prec_init[3], add_prec = prec_add_sample[i-1,3,gen],
                                  xtp1 = C_samples[1,3,i-1,gen], data = C_obs[1,3], type = 'first', init_mean = init_mean[3])
    
    C.0[i,4,gen] <- get_latent_sample(atp = (1 - .5*exp(.5*p[i,10,gen]*(mint[1] + maxt[1]))*(p[i,8,gen]+p[i,1,gen])), 
                                  btp = p[i,5,gen]*C_samples[1,1,i-1,gen] + p[i,7,gen]*C_samples[1,3,i-1,gen],
                                  initial_prec = prec_init[4], add_prec = prec_add_sample[i-1,4,gen],
                                  xtp1 = C_samples[1,4,i-1,gen], data = C_obs[1,4], type = 'first', init_mean = init_mean[4])
    
    C.0[i,5,gen] <- get_latent_sample(atp = 1 - .5*p[i,9,gen]*exp(.5*p[i,10,gen]*(maxt[1] + mint[1])), 
                                  btp = p[i,6,gen]*C_samples[1,2,i-1,gen] + .5*p[i,1,gen]*exp(.5*p[i,10,gen]*(mint[1] + maxt[1]))*C_samples[1,4,i-1,gen],
                                  initial_prec = prec_init[5], add_prec = prec_add_sample[i-1,5,gen],
                                  xtp1 = C_samples[1,5,i-1,gen], data = C_obs[1,5], type = 'first', init_mean = init_mean[5])
    
    C_samples[,,i,gen] <- C_samples[,,i-1,gen]
    
    C_samples[1,1,i,gen] <- get_latent_sample(at = 1 - p[i,5,gen], atp = 1 - p[i,5,gen], bt = G[i-1,1,gen]*(1-p[i,2,gen])*p[i,3,gen], 
                                          btp = G[i-1,2,gen]*(1-p[i,2,gen])*p[i,3,gen], 
                                          add_prec = prec_add_sample[i-1,1,gen],
                                          obs_prec = prec_obs[1], xtp1 = C_samples[2, 1,i,gen], xtm1 = C.0[i,1,gen],
                                          data = C_obs[1,1], type = 'middle')
    
    C_samples[1,2,i,gen] <- get_latent_sample(at = 1 - p[i,6,gen], atp = 1 - p[i,6,gen], bt = G[i-1,1,gen]*(1 - p[i,2,gen])*(1 - p[i,3,gen])*(1 - p[i,4,gen]), 
                                          btp = G[i-1,2,gen]*(1 - p[i,2,gen])*(1 - p[i,3,gen])*(1 - p[i,4,gen]), 
                                          add_prec = prec_add_sample[i-1,2,gen],
                                          obs_prec = prec_obs[2], xtp1 = C_samples[2, 2,i,gen], xtm1 = C.0[i,2,gen],
                                          data = C_obs[1,2], type = 'middle')
    
    C_samples[1,3,i,gen] <- get_latent_sample(at = 1 - p[i,7,gen], atp = 1 - p[i,7,gen], bt = G[i-1,1,gen]*(1-p[i,2,gen])*(1-p[i,3,gen])*p[i,4,gen], 
                                          btp = G[i-1,2,gen]*(1-p[i,2,gen])*(1-p[i,3,gen])*p[i,4,gen], 
                                          add_prec = prec_add_sample[i-1,3,gen],
                                          obs_prec = prec_obs[3], xtp1 = C_samples[2, 3,i,gen], xtm1 = C.0[i,3,gen],
                                          data = C_obs[1,3], type = 'middle')
    
    C_samples[1,4,i,gen] <- get_latent_sample(at = (1 - .5*exp(.5*p[i,10,gen]*(mint[1] + maxt[1]))*(p[i,8,gen]+p[i,1,gen])), 
                                          atp = (1 - .5*exp(.5*p[i,10,gen]*(mint[2] + maxt[2]))*(p[i,8,gen]+p[i,1,gen])), 
                                          bt = p[i,5,gen]*C_samples[1,1,i,gen] + p[i,7,gen]*C_samples[1,3,i,gen], 
                                          btp = p[i,5,gen]*C_samples[2,1,i,gen] + p[i,7,gen]*C_samples[2,3,i,gen], 
                                          add_prec = prec_add_sample[i-1,4,gen],
                                          obs_prec = prec_obs[4], xtp1 = C_samples[2, 4,i,gen], xtm1 = C.0[i,4,gen],
                                          data = C_obs[1,4], type = 'middle')
    
    C_samples[1,5,i,gen] <- get_latent_sample(at = 1 - .5*p[i,9,gen]*exp(.5*p[i,10,gen]*(maxt[1] + mint[1])), 
                                          atp = 1 - .5*p[i,9,gen]*exp(.5*p[i,10,gen]*(maxt[2] + mint[2])), 
                                          bt = p[i,6,gen]*C_samples[1,2,i,gen] + .5*p[i,1,gen]*exp(.5*p[i,10,gen]*(mint[1] + maxt[1]))*C_samples[1,4,i,gen], 
                                          btp = p[i,6,gen]*C_samples[2,2,i,gen] + .5*p[i,1,gen]*exp(.5*p[i,10,gen]*(mint[2] + maxt[2]))*C_samples[2,4,i,gen], 
                                          add_prec = prec_add_sample[i-1,5,gen],
                                          obs_prec = prec_obs[5], xtp1 = C_samples[2, 5,i,gen], xtm1 = C.0[i,5,gen],
                                          data = C_obs[1,5], type = 'middle')
    for (k in 2:729){
      ## Gibbs steps for "middle" latent states
      
      C_samples[k,1,i,gen] <- get_latent_sample(at = 1 - p[i,5,gen], atp = 1 - p[i,5,gen], bt = G[i-1,k,gen]*(1-p[i,2,gen])*p[i,3,gen], btp = G[i-1,k+1,gen]*(1-p[i,2,gen])*p[i,3,gen],
                                            obs_prec = prec_obs[1], add_prec = prec_add_sample[i-1,1,gen], xtm1 = C_samples[k-1,1,i,gen], 
                                            xtp1 = C_samples[k+1,1,i,gen], data = C_obs[k,1], type = 'middle')
      
      C_samples[k,2,i,gen] <- get_latent_sample(at = 1 - p[i,6,gen], atp = 1 - p[i,6,gen], bt = G[i-1,k,gen]*(1 - p[i,2,gen])*(1 - p[i,3,gen])*(1 - p[i,4,gen]), 
                                            btp = G[i-1,k+1,gen]*(1 - p[i,2,gen])*(1 - p[i,3,gen])*(1 - p[i,4,gen]),
                                            obs_prec = prec_obs[2], add_prec = prec_add_sample[i-1,2,gen], xtm1 = C_samples[k-1,2,i,gen], 
                                            xtp1 = C_samples[k+1,2,i,gen], data = C_obs[k,2], type = 'middle')
      
      C_samples[k,3,i,gen] <- get_latent_sample(at = 1 - p[i,7,gen], atp = 1 - p[i,7,gen], bt = G[i-1,k,gen]*(1-p[i,2,gen])*(1-p[i,3,gen])*p[i,4,gen], 
                                            btp = G[i-1,k+1,gen]*(1-p[i,2,gen])*(1-p[i,3,gen])*p[i,4,gen],
                                            obs_prec = prec_obs[3], add_prec = prec_add_sample[i-1,3,gen], xtm1 = C_samples[k-1,3,i,gen], 
                                            xtp1 = C_samples[k+1,3,i,gen], data = C_obs[k,3], type = 'middle')
      
      C_samples[k,4,i,gen] <- get_latent_sample(at = (1 - .5*exp(.5*p[i,10,gen]*(mint[k] + maxt[k]))*(p[i,8,gen]+p[i,1,gen])),
                                            atp = (1 - .5*exp(.5*p[i,10,gen]*(mint[k+1] + maxt[k+1]))*(p[i,8,gen]+p[i,1,gen])), 
                                            bt = p[i,5,gen]*C_samples[k-1,1,i,gen] + p[i,7,gen]*C_samples[k-1,3,i,gen], 
                                            btp = p[i,5,gen]*C_samples[k,1,i,gen] + p[i,7,gen]*C_samples[k,3,i,gen],
                                            obs_prec = prec_obs[4], add_prec = prec_add_sample[i-1,4,gen], xtm1 = C_samples[k-1,4,i,gen], 
                                            xtp1 = C_samples[k+1,4,i,gen], data = C_obs[k,4], type = 'middle')
      
      C_samples[k,5,i,gen] <- get_latent_sample(at = 1 - .5*p[i,9,gen]*exp(.5*p[i,10,gen]*(maxt[k] + mint[k])),
                                            atp = 1 - .5*p[i,9,gen]*exp(.5*p[i,10,gen]*(maxt[k+1] + mint[k+1])), 
                                            bt = p[i,6,gen]*C_samples[k-1,2,i,gen] + .5*p[i,1,gen]*exp(.5*p[i,10,gen]*(mint[k] + maxt[k]))*C_samples[k-1,4,i,gen], 
                                            btp = p[i,6,gen]*C_samples[k,2,i,gen] + .5*p[i,1,gen]*exp(.5*p[i,10,gen]*(mint[k+1] + maxt[k+1]))*C_samples[k,4,i,gen],
                                            obs_prec = prec_obs[5], add_prec = prec_add_sample[i-1,5,gen], xtm1 = C_samples[k-1,5,i,gen], 
                                            xtp1 = C_samples[k+1,5,i,gen], data = C_obs[k,5], type = 'middle')
    }
    ## Gibbs steps for final state
    C_samples[730, 1, i,gen] <- get_latent_sample(at = 1 - p[i,5,gen], bt = G[i-1,k,gen]*(1-p[i,2,gen])*p[i,3,gen], obs_prec = prec_obs[1],
                                              add_prec = prec_add_sample[i-1,1,gen], xtm1 = C_samples[729, 1, i,gen], data = C_obs[730, 1], type = 'last')
    
    C_samples[730, 2, i,gen] <- get_latent_sample(at = 1 - p[i,6,gen], bt = G[i-1,k,gen]*(1 - p[i,2,gen])*(1 - p[i,3,gen])*(1 - p[i,4,gen]), obs_prec = prec_obs[2],
                                              add_prec = prec_add_sample[i-1,2,gen], xtm1 = C_samples[729, 2, i,gen], data = C_obs[730, 2], type = 'last')
    
    C_samples[730, 3, i,gen] <- get_latent_sample(at = 1 - p[i,7,gen], bt = G[i-1,k,gen]*(1-p[i,2,gen])*(1-p[i,3,gen])*p[i,4,gen], obs_prec = prec_obs[3],
                                              add_prec = prec_add_sample[i-1,3,gen], xtm1 = C_samples[729, 3, i,gen], data = C_obs[730, 3], type = 'last')
    
    C_samples[730, 4, i,gen] <- get_latent_sample(at = (1 - .5*exp(.5*p[i,10,gen]*(mint[730] + maxt[730]))*(p[i,8,gen]+p[i,1,gen])), 
                                              bt = p[i,5,gen]*C_samples[729,1,i,gen] + p[i,7,gen]*C_samples[729,3,i,gen], obs_prec = prec_obs[4],
                                              add_prec = prec_add_sample[i-1,4,gen], xtm1 = C_samples[729, 4, i,gen], data = C_obs[730, 4], type = 'last')
    
    C_samples[730, 5, i,gen] <- get_latent_sample(at = 1 - .5*p[i,9,gen]*exp(.5*p[i,10,gen]*(maxt[730] + mint[730])), 
                                              bt = p[i,6,gen]*C_samples[729,2,i,gen] + .5*p[i,1,gen]*exp(.5*p[i,10,gen]*(mint[730] + maxt[730]))*C_samples[729,4,i,gen], 
                                              obs_prec = prec_obs[5], add_prec = prec_add_sample[i-1,5,gen], xtm1 = C_samples[729, 5, i,gen], data = C_obs[730, 5], 
                                              type = 'last')
    Gmeans <- (phi_g*Gpred[,gen] + tau_g*G_obs) / (phi_g + tau_g)
    G_sd <- 1 / sqrt(phi_g + tau_g)
    G[i,,gen] <- rnorm(730, mean = Gmeans, sd = G_sd)
  }
  
  if (i %% 500 == 1){
    cat(paste(i, 'iterations done \n'))
    #print(sd)
    if (i >= (500)){
      print(paste0('Process parameter acceptance over last 500 iterations is ',mean(acceptance2[((i-500):i),])))
      print(paste0('Precision parameter acceptance over last 500 iterations is ', mean(prec_Accept2[(i-500) : i,])))
    }
  }
}
## ADD GPP PRECISION
#})
time2 <- Sys.time()
time2 - time1

for (i in 1:11){
  for (j in 1:num_gen){
    if (j == 1) {
      plot(p[1:chain_length, i, j], type = 'l', col = j, ylim = c(plower[i], pupper[i]))
    } else{
      points(p[1:chain_length, i, j], type = 'l', col = j)
    }
    
  }
  readline('')
}

for (i in 1:5){
  for (j in 1:num_gen){
    if (j == 1) {
      plot(prec_add_sample[burn:chain_length, i, j], type = 'l', col = j)
    } else{
      points(prec_add_sample[burn:chain_length, i, j], type = 'l', col = j)
    }
    
  }
  readline('')
}
burn <- 6000
p_hists <- p[burn:chain_length, ,1]
for (z in 2:num_gen){
  p_hists <- rbind(p_hists, p[burn:chain_length, ,z])
}

for (i in 1:11){
  hist(p_hists[,i])
  print(mean(p_hists[,i]))
  print(paste0('True mean ',p_e[i]))
  readline('')
}
