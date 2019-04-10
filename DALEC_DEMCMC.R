## Program: sampler.R
## Author: John W Smith Jr
## Date: Feb 19th, 2019
## Description: Runs an MCMC routine in order to estimate parameters in an ecosystem state space model, in 
## particular the DALEC model. Utilizes library Rcpp for fast likelihood evaluations for MH steps

## Changes working directory 
setwd('~/Reflex/Sampler')
## Create vector for initial condition mean + stdev, as well as variances and precisions for obs and add
init_var <- c(1,4,.5,.5,4)
init_mean <- c(100,9200,100,20,11000)
var_add <- c(2, 9, (2)^2, .25, 16)
var_obs <- c(4, 16, 1.5^2, 1^2, 25)
prec_add <- 1 / var_add
prec_obs <- 1 / var_obs
prec_init <- 1 / init_var
## Sources R scripts to generate data
source('./DALEC_data_gen.R')
#source('./DALEC_nn.R')

# Change directory back
setwd('~/Reflex/Sampler')
# Load c++ functions through Rcpp
if (!("Rcpp" %in% installed.packages())) {
  install.packages('Rcpp')
}
library(Rcpp)
## Load matrixStats library for colMeans
if (!('matrixStats' %in% installed.packages())) {
  install.packages('matrixStats')
}
library(matrixStats)
## Load laGP for GP initializations of latent states
if (!('laGP' %in% installed.packages())) {
  install.packages('laGP')
}
library(laGP)

sourceCpp('./src/lmvnd_ll.cpp')
sourceCpp('./src/likelihood_SS.cpp')
sourceCpp('./src/lgamma.cpp')
sourceCpp('./src/Cpred.cpp')
source('cpp_ss_likelihood.R')
source('get_latent_sample.R')
source('latent_initialization.R')
source('DEMCMC_proposal.R')

## ------------ Sampler Test Run -----------##
# Set carbon stock variables from data generation scripts
Cf <- Cf_e
Cw <- Cw_e
Cr <- Cr_e
Clit <- Clit_e
Csom <- Csom_e
# Set driver variables
mint <- drivers_ev$mint[1:730] 
maxt <- drivers_ev$maxt[1:730] 
lat <- lat_e 
Nit <- Nit_e 
LMA <- LMA_e
Ca <- drivers_ev$ca[1:730] 
yearday <- drivers_ev$yearday[1:730] 
rad <- drivers_ev$rad[1:730]
# Set chain length and burn in for MCMC
chain_length <- 20000
burn <- 2000


## ----- Preparing MCMC ------ ##
# Compute LAI estimate to be fed into p_11 optimization routine
LAI <- pmax(.1, Cf.obs/LMA_e)
# Column bind matrices and build initial Cpred matrix
C_obs <- cbind(Cf.obs, Cw.obs, Cr.obs, Clit.obs, Csom.obs)
missing_data_ind = seq(from = 1, to = 730, by = 3)
C_obs[-missing_data_ind,] <- NA
C <- cbind(Cf, Cw, Cr, Clit, Csom)

# Initalize matrix for p draws to be stored
p <- matrix(NA, nrow = chain_length, ncol = 11)
# Set vector of upper and lower bounds for accept reject sampling
plower <- c(10e-6, .2, .01, .01, 10e-4, 10e-6, 10e-6, 10e-5, 10e-6, .05, 5)
pupper <- c(.01, .7, .5, .5, .1, .01, .01, .1, .01, .2, 20)
# Fix p11 and initalize first guesses at the midpoint
p[1,1:11] <- .5*(pupper[1:11] + plower[1:11])
buildCpred <- Cpred(p[1,], C, LAI, maxt, mint, Ca, lat_e, yearday, Nit_e, rad, 730)
G = buildCpred$G
Cpredi <- cbind(buildCpred$Cfpred, buildCpred$Cwpred, buildCpred$Crpred,
                buildCpred$Clitpred, buildCpred$Csompred)
rm(buildCpred)

## ----- Performing MCMC ------ ##
## Accept Reject Normal for p_i ##

# create matrix for acceptance for p_i 
acceptance <- matrix(0, nrow = chain_length, ncol = 11)
acceptance2 <- rep(0, chain_length)
prec_Accept <- matrix(0, nrow = chain_length, ncol = 5)
prec_Accept2 <- rep(0, chain_length)
prec_sd <- rep(.2, 5)

# Initialize arrays for C_samples
C_samples <- array(NA, dim = c(730, 5, chain_length))

# Initialize first entry to observed data
C_samples[, , 1] <- latent_initialization(C_obs, missing_data_ind, 1, 730, type = c('gp', 'gp', 'pwl', 'pwl', 'gp'))

prec_add_sample = matrix(1, nrow = chain_length, ncol = 5)
sd = (log(pupper) - log(plower)) / 10
num_init <- 500


time1 <- Sys.time()
#library(profvis)
#profvis({
for (i in 2:chain_length){
  ## Build Cpred matrix
  p[i,] <- p[i-1,]
  buildCpredi <- Cpred(p[i,], C_samples[, , i-1], LAI_ = pmax(.1, C_samples[,1,i-1]/LMA_e), maxt, mint, Ca, lat_e, yearday, Nit_e, rad, 730)
  Cpredi <- cbind(buildCpredi$Cfpred, buildCpredi$Cwpred, buildCpredi$Crpred,
                  buildCpredi$Clitpred, buildCpredi$Csompred)
  Gpred <- buildCpredi$G
  rm(buildCpredi)
  ## Get initial samples to use in DEMCMC
  if (i <= num_init){
    ## Loop over p_i
    for (j in 1:11){
      p_star <- p[i,]
        tmp <- 0
        ## Accept Reject Sampling step in log space
        while (tmp < log(plower[j]) | tmp > log(pupper[j])){
          tmp <- rnorm(1, log(p[i,j]), sd[j])
        }
        p_star[j] <- exp(tmp)
        #print(p_star[j])
        ## Build new Cpred matrix
        buildCpredstar <- Cpred(p_star, C_samples[, , i-1], LAI_ = pmax(.1, C_samples[,1,i-1]/LMA_e), maxt, mint, Ca, lat_e, yearday, Nit_e, rad, 730)
        Cpred_star <- cbind(buildCpredstar$Cfpred, buildCpredstar$Cwpred, buildCpredstar$Crpred,
                            buildCpredstar$Clitpred, buildCpredstar$Csompred)
        Gpred_star <- buildCpredstar$G
        rm(buildCpredstar)
        ## Compute threshold, perform Metropolis step
        threshold <- cpp_ss_likelihood(C_obs, C_samples[, , i-1], 1 / prec_obs, Cpred_star, 1 / prec_add_sample[i-1,], 
                                       c(100,9200,100,20,11000), init_var, 730, Gpred_star, G_e, ind = missing_data_ind - 1) - 
          cpp_ss_likelihood(C_obs, C_samples[, , i-1], 1 / prec_obs, Cpredi, 1 / prec_add_sample[i-1,], 
                            c(100,9200,100,20,11000), init_var, 730, Gpred, G_e, ind = missing_data_ind - 1)
        if (log(runif(1,0,1)) < threshold){
          p[i,j] <- p_star[j]
          Cpredi <- Cpred_star
          G <- Gpred_star
        }
      } 
     
  } else if (i>num_init) {
      ## DEMCMC step 
      p[i,] <- p[i-1,]
      buildCpredi <- Cpred(p[i,], C_samples[, , i-1], LAI_ = pmax(.1, C_samples[,1,i-1]/LMA_e), maxt, mint, Ca, lat_e, yearday, Nit_e, rad, 730)
      Cpredi <- cbind(buildCpredi$Cfpred, buildCpredi$Cwpred, buildCpredi$Crpred,
                    buildCpredi$Clitpred, buildCpredi$Csompred)
      Gpred <- buildCpredi$G
      rm(buildCpredi)
      p_star <- DEMCMC_proposal(i, p, lower = plower, upper = pupper, jitter.max = 1e-6, gamma = .2)
      
      buildCpredstar <- Cpred(p_star, C_samples[, , i-1], LAI_ = pmax(.1, C_samples[,1,i-1]/LMA_e), maxt, mint, Ca, lat_e, yearday, Nit_e, rad, 730)
      Cpred_star <- cbind(buildCpredstar$Cfpred, buildCpredstar$Cwpred, buildCpredstar$Crpred,
                        buildCpredstar$Clitpred, buildCpredstar$Csompred)
      Gpred_star <- buildCpredstar$G
      rm(buildCpredstar)
      
      threshold <- cpp_ss_likelihood(C_obs, C_samples[, , i-1], 1 / prec_obs, Cpred_star, 1 / prec_add_sample[i-1,], 
                                   c(100,9200,100,20,11000), init_var, 730, Gpred_star, G_e, ind = missing_data_ind - 1) - 
        cpp_ss_likelihood(C_obs, C_samples[, , i-1], 1 / prec_obs, Cpredi, 1 / prec_add_sample[i-1,], 
                        c(100,9200,100,20,11000), init_var, 730, Gpred, G_e, ind = missing_data_ind - 1)
      if (log(runif(1,0,1)) < threshold){
        acceptance2[i] <- 1
        p[i,] <- p_star
        Cpredi <- Cpred_star
        G <- Gpred_star
      }
    } 
  prec_add_sample[i,] <- prec_add_sample[i-1,]
  if (i <= num_init){
    ## Get initial estimates to use in DEMCMC step
    for (k in 1:5){
      prec_star <- prec_add_sample[i,]
      tmp <- -1
      while (tmp < 0){
        tmp <- rnorm(1, prec_star[k], sd = prec_sd[k])
      }
      prec_star[k] <- tmp
      threshold <- cpp_ss_likelihood(C_obs, C_samples[, , i-1], 1 / prec_obs, Cpredi, 1 / prec_star, 
                                   c(100,9200,100,20,11000), init_var, 730, Gpred, G_e, ind = missing_data_ind - 1) - 
        cpp_ss_likelihood(C_obs, C_samples[, , i-1], 1 / prec_obs, Cpredi, 1 / prec_add_sample[i-1,], 
                        c(100,9200,100,20,11000), init_var, 730, Gpred, G_e, ind = missing_data_ind - 1)
      if (log(runif(1,0,1)) < threshold){
        prec_add_sample[i,] <- prec_star
      }
    }
  } else{
      ## DEMCMC step for precisions
      prec_star <- DEMCMC_proposal(i, prec_add_sample, lower = rep(0,dim(prec_add_sample)[2]), jitter.max = 1e-6, gamma = .2)
      threshold <- cpp_ss_likelihood(C_obs, C_samples[, , i-1], 1 / prec_obs, Cpredi, 1 / prec_star, 
                                   c(100,9200,100,20,11000), init_var, 730, Gpred, G_e, ind = missing_data_ind - 1) - 
        cpp_ss_likelihood(C_obs, C_samples[, , i-1], 1 / prec_obs, Cpredi, 1 / prec_add_sample[i-1,], 
                        c(100,9200,100,20,11000), init_var, 730, Gpred, G_e, ind = missing_data_ind - 1)
      if (log(runif(1,0,1)) < threshold){
        prec_add_sample[i,] <- prec_star
        prec_Accept2[i] <- 1
      }
    }
  ## Gibbs Step for C_i for the initial state
  C_samples[,,i] <- C_samples[,,i-1]
  
  C_samples[1,1,i] <- get_latent_sample(at = 1 - p[i,5], bt = G[1]*(1-p[i,2])*p[i,3], initial_prec = prec_init[1], add_prec = prec_add_sample[i-1,1],
                                        obs_prec = prec_obs[1], init_mean = init_mean[1], xtp1 = C_samples[2, 1,i], 
                                        data = C_obs[1,1], type = 'first')
  
  C_samples[1,2,i] <- get_latent_sample(at = 1 - p[i,6], bt = G[1]*(1 - p[i,2])*(1 - p[i,3])*(1 - p[i,4]), initial_prec = prec_init[2], 
                                        add_prec = prec_add_sample[i-1,2], obs_prec = prec_obs[2], init_mean = init_mean[2], 
                                        xtp1 = C_samples[2, 2,i], data = C_obs[1,2], type = 'first')
  
  C_samples[1,3,i] <- get_latent_sample(at = 1 - p[i,7], bt = G[1]*(1-p[i,2])*(1-p[i,3])*p[i,4], initial_prec = prec_init[3], 
                                        add_prec = prec_add_sample[i-1,3], obs_prec = prec_obs[3], init_mean = init_mean[3], 
                                        xtp1 = C_samples[2, 3,i], data = C_obs[1,3], type = 'first')
  
  C_samples[1,4,i] <- get_latent_sample(at = (1 - .5*exp(.5*p[i,10]*(mint[1] + maxt[1]))*(p[i,8]+p[i,1])), 
                                        bt = p[i,5]*C_samples[1,1,i] + p[i,7]*C_samples[1,3,i], initial_prec = prec_init[4], 
                                        add_prec = prec_add_sample[i-1,4], obs_prec = prec_obs[4], init_mean = init_mean[4], 
                                        xtp1 = C_samples[2, 4,i], data = C_obs[1,4], type = 'first')
  
  C_samples[1,5,i] <- get_latent_sample(at = 1 - .5*p[i,9]*exp(.5*p[i,10]*(maxt[1] + mint[1])), 
                                        bt = p[i,6]*C_samples[1,2,i] + .5*p[i,1]*exp(.5*p[i,10]*(mint[1] + maxt[1]))*C_samples[1,4,i], 
                                        initial_prec = prec_init[5], add_prec = prec_add_sample[i-1,5], obs_prec = prec_obs[5], 
                                        init_mean = init_mean[5], xtp1 = C_samples[2, 5,i], data = C_obs[1,5], type = 'first')
  for (k in 2:729){
    ## Gibbs steps for "middle" latent states
    C_samples[k,1,i] <- get_latent_sample(at = 1 - p[i,5], atp = 1 - p[i,5], bt = G[k]*(1-p[i,2])*p[i,3], btp = G[k+1]*(1-p[i,2])*p[i,3],
                                          obs_prec = prec_obs[1], add_prec = prec_add_sample[i-1,1], xtm1 = C_samples[k-1,1,i], 
                                          xtp1 = C_samples[k+1,1,i], data = C_obs[k,1], type = 'middle')
    
    C_samples[k,2,i] <- get_latent_sample(at = 1 - p[i,6], atp = 1 - p[i,6], bt = G[k]*(1 - p[i,2])*(1 - p[i,3])*(1 - p[i,4]), 
                                          btp = G[k+1]*(1 - p[i,2])*(1 - p[i,3])*(1 - p[i,4]),
                                          obs_prec = prec_obs[2], add_prec = prec_add_sample[i-1,2], xtm1 = C_samples[k-1,2,i], 
                                          xtp1 = C_samples[k+1,2,i], data = C_obs[k,2], type = 'middle')
    
    C_samples[k,3,i] <- get_latent_sample(at = 1 - p[i,7], atp = 1 - p[i,7], bt = G[k]*(1-p[i,2])*(1-p[i,3])*p[i,4], 
                                          btp = G[k+1]*(1-p[i,2])*(1-p[i,3])*p[i,4],
                                          obs_prec = prec_obs[3], add_prec = prec_add_sample[i-1,3], xtm1 = C_samples[k-1,3,i], 
                                          xtp1 = C_samples[k+1,3,i], data = C_obs[k,3], type = 'middle')
    
    C_samples[k,4,i] <- get_latent_sample(at = (1 - .5*exp(.5*p[i,10]*(mint[k] + maxt[k]))*(p[i,8]+p[i,1])),
                                          atp = (1 - .5*exp(.5*p[i,10]*(mint[k+1] + maxt[k+1]))*(p[i,8]+p[i,1])), 
                                          bt = p[i,5]*C_samples[k-1,1,i] + p[i,7]*C_samples[k-1,3,i], 
                                          btp = p[i,5]*C_samples[k,1,i] + p[i,7]*C_samples[k,3,i],
                                          obs_prec = prec_obs[4], add_prec = prec_add_sample[i-1,4], xtm1 = C_samples[k-1,4,i], 
                                          xtp1 = C_samples[k+1,4,i], data = C_obs[k,4], type = 'middle')
    
    C_samples[k,5,i] <- get_latent_sample(at = 1 - .5*p[i,9]*exp(.5*p[i,10]*(maxt[k] + mint[k])),
                                          atp = 1 - .5*p[i,9]*exp(.5*p[i,10]*(maxt[k+1] + mint[k+1])), 
                                          bt = p[i,6]*C_samples[k-1,2,i] + .5*p[i,1]*exp(.5*p[i,10]*(mint[k] + maxt[k]))*C_samples[k-1,4,i], 
                                          btp = p[i,6]*C_samples[k,2,i] + .5*p[i,1]*exp(.5*p[i,10]*(mint[k+1] + maxt[k+1]))*C_samples[k,4,i],
                                          obs_prec = prec_obs[5], add_prec = prec_add_sample[i-1,5], xtm1 = C_samples[k-1,5,i], 
                                          xtp1 = C_samples[k+1,5,i], data = C_obs[k,5], type = 'middle')
  }
  ## Gibbs steps for final state
  C_samples[730, 1, i] <- get_latent_sample(at = 1 - p[i,5], bt = G[k]*(1-p[i,2])*p[i,3], obs_prec = prec_obs[1],
                                            add_prec = prec_add_sample[i-1,1], xtm1 = C_samples[729, 1, i], data = C_obs[730, 1], type = 'last')
  
  C_samples[730, 2, i] <- get_latent_sample(at = 1 - p[i,6], bt = G[k]*(1 - p[i,2])*(1 - p[i,3])*(1 - p[i,4]), obs_prec = prec_obs[2],
                                            add_prec = prec_add_sample[i-1,2], xtm1 = C_samples[729, 2, i], data = C_obs[730, 2], type = 'last')
  
  C_samples[730, 3, i] <- get_latent_sample(at = 1 - p[i,7], bt = G[k]*(1-p[i,2])*(1-p[i,3])*p[i,4], obs_prec = prec_obs[3],
                                            add_prec = prec_add_sample[i-1,3], xtm1 = C_samples[729, 3, i], data = C_obs[730, 3], type = 'last')
  
  C_samples[730, 4, i] <- get_latent_sample(at = (1 - .5*exp(.5*p[i,10]*(mint[730] + maxt[730]))*(p[i,8]+p[i,1])), 
                                            bt = p[i,5]*C_samples[729,1,i] + p[i,7]*C_samples[729,3,i], obs_prec = prec_obs[4],
                                            add_prec = prec_add_sample[i-1,4], xtm1 = C_samples[729, 4, i], data = C_obs[730, 4], type = 'last')
  
  C_samples[730, 5, i] <- get_latent_sample(at = 1 - .5*p[i,9]*exp(.5*p[i,10]*(maxt[730] + mint[730])), 
                                            bt = p[i,6]*C_samples[729,2,i] + .5*p[i,1]*exp(.5*p[i,10]*(mint[730] + maxt[730]))*C_samples[729,4,i], 
                                            obs_prec = prec_obs[5], add_prec = prec_add_sample[i-1,5], xtm1 = C_samples[729, 5, i], data = C_obs[730, 5], 
                                            type = 'last')
  if (i %% 500 == 0){
    cat(paste(i, 'iterations done \n'))
    #print(sd)
    if (i >= (burn + 500)){
      print(mean(acceptance2[((i-500):i)]))
    }
  }
}
## ADD GPP PRECISION
#})
time2 <- Sys.time()
time2 - time1

sim_directory = './Simulations/MissingEveryThird/'
burn = 10000
for (i in 1:11){
  #png(filename = paste0(sim_directory, paste0(paste0('/p',i),'.png')))
  hist(p[burn:chain_length, i], breaks = 25,  xlab = paste0('p',i),
       main = paste0(paste0(paste0(paste0('Histogram of p',i)), ', true value = '), p_e[i]))
  #dev.off()
  cat(paste('The mean is', mean(p[burn:chain_length, i])))
  cat('\n')
  cat(paste('The true value is', p_e[i]))
  cat('\n')
  readline('Press enter to continue')
  plot(p[1:chain_length, i], type = 'l')
  readline('Press enter to continue')
}

for (i in 1:11){
  #png(filename = paste0(sim_directory, paste0(paste0('/p_dens/p',i),'density.png')))
  plot(density(p[burn:chain_length, i]),  main = paste0(paste0(paste0(paste0('Estimated Density of p',i)), ', true value = '), 
                                                        p_e[i]), xlab = paste0('p',i))
  readline('\n')
  #abline(v = mean(p[burn:chain_length, i]), lty = 2)
  #dev.off()
}

for (j in 1:5){
  #png(filename = paste0(paste0('./Plots/sd_obs',j),'.png'))
  #png(filename = paste0(sim_directory, paste0(paste0('/prec_',j),'.png')))
  hist(prec_add_sample[burn:chain_length, j], breaks = 25, main = paste0('True Value is ', prec_add[j]))
  #dev.off()
  cat(paste('The mean is', mean(prec_add_sample[burn:chain_length,j])))
  cat('\n')
  readline('Press Enter')
}

for (j in 1:5){
  png(paste0(sim_directory, paste0(paste0('/prec_dens/Prec_Obs_',colnames(C)[j]),'density.png')))
  plot(density(prec_obs_sample[burn:chain_length, j]), xlab = paste0('Observation Precision for ', colnames(C)[j]),
       main = paste0(paste0(paste0('Estimated Density for Prec Obs ', colnames(C)[j]), '\nTrue Value = '), prec_obs[j]))
  abline(v = mean(prec_obs_sample[burn:chain_length, j]), lty = 2)
  dev.off()
  #readline('')
}

png(filename = './Plots/pairs_postburn.png', width = 1100, height = 1100)
pairs(p[burn:chain_length, 1:11])
dev.off()

png(filename = paste0(sim_directory, 'latent_plots/latent_Cf.png'), width = 600, height = 600)
plot(C[,1], type = 'l', ylab = colnames(C)[1], xlab = 'Time (Days)', main = 'Comparison of Estimated vs\n True Latent States')
points(rowMeans(C_samples[,1,burn:chain_length]), type = 'l', col = 'red')
points(rowQuantiles(C_samples[,1,burn:chain_length], prob = .05), type = 'l', col = 'red', lty = 2)
points(rowQuantiles(C_samples[,1,burn:chain_length], prob = .95), type = 'l', col = 'red', lty = 2)
legend('topleft', legend = c('True Latent State', 'Estimated Latent State', 'Credible Bounds'), col = c('black', 'red', 'red'), lty = c(1,1,2))
dev.off()

png(filename = paste0(sim_directory, 'latent_plots/latent_Cw.png'), width = 600, height = 600)
plot(C[,2], type = 'l', ylab = colnames(C)[2], xlab = 'Time (Days)', main = 'Comparison of Estimated vs\n True Latent States')
points(rowMeans(C_samples[,2,burn:chain_length]), type = 'l', col = 'red')
points(rowQuantiles(C_samples[,2,burn:chain_length], prob = .05), type = 'l', col = 'red', lty = 2)
points(rowQuantiles(C_samples[,2,burn:chain_length], prob = .95), type = 'l', col = 'red', lty = 2)
legend('topleft', legend = c('True Latent State', 'Estimated Latent State', 'Credible Bounds'), col = c('black', 'red', 'red'), lty = c(1,1,2))
dev.off()
#points(C_samples[,2,10000], type = 'l', col = 'red')

png(filename = paste0(sim_directory, 'latent_plots/latent_Cr.png'), width = 600, height = 600)
plot(C[,3], type = 'l', ylab = colnames(C)[3], xlab = 'Time (Days)', main = 'Comparison of Estimated vs\n True Latent States')
points(rowMeans(C_samples[,3,burn:chain_length]), type = 'l', col = 'red')
points(rowQuantiles(C_samples[,3,burn:chain_length], prob = .05), type = 'l', col = 'red', lty = 2)
points(rowQuantiles(C_samples[,3,burn:chain_length], prob = .95), type = 'l', col = 'red', lty = 2)
legend('topleft', legend = c('True Latent State', 'Estimated Latent State', 'Credible Bounds'), col = c('black', 'red', 'red'), lty = c(1,1,2))
dev.off()
#points(C_samples[,3,10000], type = 'l', col = 'red')

png(filename = paste0(sim_directory, '/latent_Clit.png'), width = 1100, height = 1100)
plot(C[,4], type = 'l')
points(rowMeans(C_samples[,4,burn:chain_length]), type = 'l', col = 'red')
points(rowQuantiles(C_samples[,4,burn:chain_length], prob = .05), type = 'l', col = 'red', lty = 2)
points(rowQuantiles(C_samples[,4,burn:chain_length], prob = .95), type = 'l', col = 'red', lty = 2)
dev.off()

png(filename = paste0(sim_directory, 'latent_plots/latent_Csom.png'), width = 600, height = 600)
plot(C[,5], type = 'l', ylab = colnames(C)[5], xlab = 'Time (Days)', main = 'Comparison of Estimated vs\n True Latent States', lwd = 2)
points(rowMeans(C_samples[,5,burn:chain_length]), type = 'l', col = 'red', lwd = 2)
points(rowQuantiles(C_samples[,5,burn:chain_length], prob = .05), type = 'l', col = 'red', lty = 2, lwd = 2)
points(rowQuantiles(C_samples[,5,burn:chain_length], prob = .95), type = 'l', col = 'red', lty = 2, lwd = 2)
legend('topleft', legend = c('True Latent State', 'Estimated Latent State', 'Credible Bounds'), col = c('black', 'red', 'red'), lty = c(1,1,2))
dev.off()

## Next:
## Different process models
## Generalize to multiple plots?
## Larger missing data intervals
## Larger process model noise