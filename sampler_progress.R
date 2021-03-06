## Program: sampler_progress.R
## Author: John W Smith Jr
## Date: Feb 19th, 2019
## Description: Runs an MCMC routine in order to estimate parameters in an ecosystem state space model, in 
## particular the REFLEX model. Utilizes library Rcpp for fast likelihood evaluations for MH steps

## Changes working directory 
setwd('~/Reflex/Sampler')
## Create vector for initial condition mean + stdev, as well as variances and precisions for obs and add
init_var <- c(4,4,4,4,4)
init_mean <- c(100,9200,100,20,11000)
#var_add <- c(.25^2, 1, (.03125/4)^2, .25, 1)
var_add <- c(2, 9, (.03125/4)^2, .25, 16)
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
sourceCpp('./src/lmvnd_ll.cpp')
sourceCpp('./src/likelihood_SS.cpp')
sourceCpp('./src/lgamma.cpp')
sourceCpp('./src/Cpred.cpp')
source('cpp_ss_likelihood.R')
source('optim_p11.R')

## ------------ Sampler Test Run -----------##
# Set carbon stock variables from data generation scripts
Cf <- Cf_e
Cw <- Cw_e
Cr <- Cr_e
Clit <- Clit_e
Csom <- Csom_e
# Set standard deviations 
sd_cf <- .25
sd_cw <- 1
sd_cr <- (.03125/4)
sd_clit <- .25
sd_csom <- 1
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
chain_length <- 8000
burn <- 2000

## ----- Preparing MCMC ------ ##
# Compute LAI estimate to be fed into p_11 optimization routine
LAI <- pmax(.1, Cf.obs/LMA_e)
# Obtain estimate for p_11 based off of GPP flux data
#p_11_est <- optim_p11(1001, G_e)
# Column bind matrices and build initial Cpred matrix
C_obs <- cbind(Cf.obs, Cw.obs, Cr.obs, Clit.obs, Csom.obs)
C <- cbind(Cf, Cw, Cr, Clit, Csom)
buildCpred <- Cpred(p[1,], C, LAI, maxt, mint, Ca, lat_e, yearday, Nit_e, rad, 730)
G = buildCpred$G
Cpredi <- cbind(buildCpred$Cfpred, buildCpred$Cwpred, buildCpred$Crpred,
               buildCpred$Clitpred, buildCpred$Csompred)
rm(buildCpred)

## ----- Performing MCMC ------ ##
## Accept Reject Normal for p_i ##

# create matrix for acceptance for p_i 
acceptance <- matrix(0, nrow = chain_length, ncol = 11)

# set proposal standard deviations
sd <- rep(.5,11)
sd[c(1,6,9)] <- .3
sd[c(2,3,4,5,7,8,10)] <- .1
# Set point to start adapting proposal SD
#start.adapt = 2500

# Initialize arrays for C_samples
C_samples <- array(NA, dim = c(730, 5, chain_length))

# Initialize first entry to observed data
C_samples[, , 1] <- C_obs
# Initalize matrix for p draws to be stored
p <- matrix(NA, nrow = chain_length, ncol = 11)
# Set vector of upper and lower bounds for accept reject sampling
plower <- c(10e-6, .2, .01, .01, 10e-4, 10e-6, 10e-6, 10e-5, 10e-6, .05, 5)
pupper <- c(.01, .7, .5, .5, .1, .01, .01, .1, .01, .2, 20)
# Fix p11 and initalize first guesses at the midpoint
p[,11] <- rep(p_11_est, chain_length)
p[1,1:11] <- .5*(pupper[1:11] + plower[1:11])
p[,11] = rep(18, chain_length)
sd[11] = .05
sd[2] = .05
sd[3] = .03

prec_obs_sample = matrix(1, nrow = chain_length, ncol = 5)

time1 <- Sys.time()
#library(profvis)
#profvis({
for (i in 2:chain_length){
  ## ----- Metropolis Step for p_i ----- ##
  p[i,] <- p[i-1,]
  buildCpredi <- Cpred(p[i,], C_samples[, , i-1], LAI_ = pmax(.1, C_samples[,1,i-1]/LMA_e), maxt, mint, Ca, lat_e, yearday, Nit_e, rad, 730)
  Cpredi <- cbind(buildCpredi$Cfpred, buildCpredi$Cwpred, buildCpredi$Crpred,
                 buildCpredi$Clitpred, buildCpredi$Csompred)
  Gpred <- buildCpredi$G
  rm(buildCpredi)
  for (j in 1:11){
    p_star <- p[i,]
    tmp <- 0
    while (tmp < log(plower[j]) | tmp > log(pupper[j])){
      tmp <- rnorm(1, log(p[i,j]), sd[j])
    }
    p_star[j] <- exp(tmp)
    
    buildCpredstar <- Cpred(p_star, C_samples[, , i-1], LAI_ = pmax(.1, C_samples[,1,i-1]/LMA_e), maxt, mint, Ca, lat_e, yearday, Nit_e, rad, 730)
    Cpred_star <- cbind(buildCpredstar$Cfpred, buildCpredstar$Cwpred, buildCpredstar$Crpred,
                        buildCpredstar$Clitpred, buildCpredstar$Csompred)
    Gpred_star <- buildCpredstar$G
    rm(buildCpredstar)
    
    threshold <- cpp_ss_likelihood(C_obs, C_samples[, , i-1], 1 / prec_obs_sample[i-1,], Cpred_star, var_add, 
                                  c(100,9200,100,20,11000), c(4,4,4,4,4), 730, Gpred_star, G_e) - 
      cpp_ss_likelihood(C_obs, C_samples[, , i-1], 1 / prec_obs_sample[i-1,], Cpredi, var_add, 
                        c(100,9200,100,20,11000), c(4,4,4,4,4), 730, Gpred, G_e)
    if (log(runif(1,0,1)) < threshold){
      p[i,j] <- p_star[j]
      Cpredi <- Cpred_star
      G <- Gpred_star
      acceptance[i,j] <- 1
    }
    #if (i%%10 == 0 & i > start.adapt & mean(acceptance[(burn:i), j]) < .25){
    #  sd[j] <- sd[j] / 1.05
    #} else if (i%%50 == 0 & i > start.adapt & mean(acceptance[(burn:i), j]) > .3){
    # sd[j] <- sd[j] * 1.05
    #}
  }
  ## Gibbs Step for C_i for the initial state
  C_samples[,,i] <- C_samples[,,i-1]
  a_t <- 1 - p[i,5]
  b_t <- G[1]*(p[i,2])
  C_samples[1,1,i] <- rnorm(1, 
                            mean = (prec_obs_sample[i-1,1]*C_obs[1,1] + prec_add[1]*a_t*C_samples[2, 1,i] - 
                                      prec_add[1]*a_t*b_t + prec_init[1]*init_mean[1] ) / 
                              (prec_obs_sample[i-1,1] + (a_t^2)*prec_add[1] + prec_init[1] ), 
                            sd = 1 / sqrt((prec_obs_sample[i-1,1] + (a_t^2)*prec_add[1] + prec_init[1] ) ))
  a_t <- 1 - p[i,6]
  b_t <- G[1]*(p[i,3])
  C_samples[1,2,i] <- rnorm(1, 
                            mean = (prec_obs_sample[i-1,2]*C_obs[1,2] + prec_add[2]*a_t*C_samples[2, 2,i] 
                                    - prec_add[2]*a_t*b_t + prec_init[2]*init_mean[2] ) / 
                              (prec_obs_sample[i-1,2] + (a_t^2)*prec_add[2] + prec_init[2] ), 
                            sd = 1 / sqrt((prec_obs_sample[i-1,2] + (a_t^2)*prec_add[2] + prec_init[2] ) ))
  
  a_t <- 1-p[i,7]
  b_t <- G[1]*p[i,4]
  C_samples[1,3,i] <- rnorm(1, 
                            mean = (prec_obs_sample[i-1,3]*C_obs[1,3] + prec_add[3]*a_t*C_samples[2, 3,i] 
                                    - prec_add[3]*a_t*b_t + prec_init[3]*init_mean[3] ) / 
                              (prec_obs_sample[i-1,3] + (a_t^2)*prec_add[3] + prec_init[3] ), 
                            sd = 1 / sqrt((prec_obs_sample[i-1,3] + (a_t^2)*prec_add[3] + prec_init[3] ) ))
  
  a_t <- (1 - .5*exp(.5*p[i,10]*(mint[1] + maxt[1]))*(p[i,8]+p[i,1]))
  b_t <- p[i,5]*C_samples[1,1,i] + p[i,7]*C_samples[1,3,i]
  C_samples[1,4,i] <- rnorm(1, 
                            mean = (prec_obs_sample[i-1,4]*C_obs[1,4] + prec_add[4]*a_t*C_samples[2, 4,i] 
                                    - prec_add[4]*a_t*b_t + prec_init[4]*init_mean[4] ) / 
                              (prec_obs_sample[i-1,4] + (a_t^2)*prec_add[4] + prec_init[4] ), 
                            sd = 1 / sqrt((prec_obs_sample[i-1,4] + (a_t^2)*prec_add[4] + prec_init[4] ) ))
  
  a_t <- 1 - .5*p[i,9]*exp(.5*p[i,10]*(maxt[1] + mint[1]))
  b_t <- p[i,6]*C_samples[1,2,i] + .5*p[i,1]*exp(.5*p[i,10]*(mint[1] + maxt[1]))*C_samples[1,4,i]
  C_samples[1,5,i] <- rnorm(1, 
                            mean = (prec_obs_sample[i-1,5]*C_obs[1,5] + prec_add[5]*a_t*C_samples[2, 5,i] 
                                    - prec_add[5]*a_t*b_t + prec_init[5]*init_mean[5] ) / 
                              (prec_obs_sample[i-1,5] + (a_t^2)*prec_add[5] + prec_init[5] ), 
                            sd = 1 / sqrt((prec_obs_sample[i-1,5] + (a_t^2)*prec_add[5] + prec_init[5] ) ))
  ## Gibbs Step for Cobs_i for the inital state
  
  for (k in 2:729){
    ## Gibbs steps for "middle" latent states
    a_t <- 1 - p[i,5]
    a_tp <- 1 - p[i,5]
    b_t <- G[k]*(1-p[i,2])*p[i,3]
    b_tp <- G[k+1]*(1-p[i,2])*p[i,3]
    
    C_samples[k,1,i] <- rnorm(n = 1, 
                      mean = (prec_add[1]*(a_t*C_samples[k-1,1,i] +b_t ) + prec_add[1]*(a_tp*C_samples[k+1,1,i]-a_tp*b_tp) 
                      + prec_obs_sample[i-1,1]*C_obs[k,1]) / (prec_add[1]*(1 + (a_tp)^2 ) + prec_obs_sample[i-1,1]),
                      sd = 1 / sqrt((prec_add[1]*(1 + (a_tp)^2 ) + prec_obs_sample[i-1,1])))
    
    a_t <- 1 - p[i,6]
    a_tp <- 1 - p[i,6]
    b_t <- G[k]*(1 - p[i,2])*(1 - p[i,3])*(1 - p[i,4])
    b_tp <- G[k+1]*(1 - p[i,2])*(1 - p[i,3])*(1 - p[i,4])
    
    C_samples[k,2,i] <- rnorm(n = 1, 
                      mean = (prec_add[2]*(a_t*C_samples[k-1,2,i] +b_t ) + prec_add[2]*(a_tp*C_samples[k+1,2,i]-a_tp*b_tp) 
                      + prec_obs_sample[i-1,2]*C_obs[k,2]) / (prec_add[2]*(1 + (a_tp)^2 ) + prec_obs_sample[i-1,2]),
                      sd = 1 / sqrt((prec_add[2]*(1 + (a_tp)^2 ) + prec_obs_sample[i-1,2])))
    
    a_t <- 1-p[i,7]
    a_tp <- 1-p[i,7]
    b_t <- G[k]*(1-p[i,2])*(1-p[i,3])*p[i,4]
    b_tp <- G[k+1]*(1-p[i,2])*(1-p[i,3])*p[i,4]
    
    C_samples[k,3,i] <- rnorm(n = 1, 
                      mean = (prec_add[3]*(a_t*C_samples[k-1,3,i] +b_t ) + prec_add[3]*(a_tp*C_samples[k+1,3,i]-a_tp*b_tp) 
                      + prec_obs_sample[i-1,3]*C_obs[k,3]) / (prec_add[3]*(1 + (a_tp)^2 ) + prec_obs_sample[i-1,3]),
                      sd = 1 / sqrt((prec_add[3]*(1 + (a_tp)^2 ) + prec_obs_sample[i-1,3])))
    
    a_t <- (1 - .5*exp(.5*p[i,10]*(mint[k] + maxt[k]))*(p[i,8]+p[i,1]))
    a_tp <- (1 - .5*exp(.5*p[i,10]*(mint[k+1] + maxt[k+1]))*(p[i,8]+p[i,1]))
      
    b_t <- p[i,5]*C_samples[k-1,1,i] + p[i,7]*C_samples[k-1,3,i]
    b_tp <- p[i,5]*C_samples[k,1,i] + p[i,7]*C_samples[k,3,i]
    
    C_samples[k,4,i] <- rnorm(n = 1, 
                      mean = (prec_add[4]*(a_t*C_samples[k-1,4,i] +b_t ) + prec_add[4]*(a_tp*C_samples[k+1,4,i]-a_tp*b_tp) 
                      + prec_obs_sample[i-1,4]*C_obs[k,4]) / (prec_add[4]*(1 + (a_tp)^2 ) + prec_obs_sample[i-1,4]),
                      sd = 1 / sqrt((prec_add[4]*(1 + (a_tp)^2 ) + prec_obs_sample[i-1,4])))
    
    a_t <- 1 - .5*p[i,9]*exp(.5*p[i,10]*(maxt[k] + mint[k]))
    a_tp <- 1 - .5*p[i,9]*exp(.5*p[i,10]*(maxt[k+1] + mint[k+1]))
    
    b_t <- p[i,6]*C_samples[k-1,2,i] + .5*p[i,1]*exp(.5*p[i,10]*(mint[k] + maxt[k]))*C_samples[k-1,4,i]
    b_tp <- p[i,6]*C_samples[k,2,i] + .5*p[i,1]*exp(.5*p[i,10]*(mint[k+1] + maxt[k+1]))*C_samples[k,4,i]
    
    
    C_samples[k,5,i] <- rnorm(n = 1, 
                      mean = (prec_add[5]*(a_t*C_samples[k-1,5,i] +b_t ) + prec_add[5]*(a_tp*C_samples[k+1,5,i]-a_tp*b_tp) 
                      + prec_obs_sample[i-1,5]*C_obs[k,5]) / (prec_add[5]*(1 + (a_tp)^2 ) + prec_obs_sample[i-1,5]),
                      sd = 1 / sqrt((prec_add[5]*(1 + (a_tp)^2 ) + prec_obs_sample[i-1,5])))
  }
  ## Gibbs steps for final state
  C_samples[730, 1, i] <- rnorm(n = 1, 
                               mean = (prec_add[1]*C_samples[729, 1, i] + prec_obs_sample[i-1,1]*C_obs[730, 1]) / (prec_add[1] + prec_obs_sample[i-1,1]),
                               sd = 1 / sqrt(prec_add[1] + prec_obs_sample[i-1,1]))
  C_samples[730, 2, i] <- rnorm(n = 1, 
                               mean = (prec_add[2]*C_samples[729, 2, i] + prec_obs_sample[i-1,2]*C_obs[730, 2]) / (prec_add[2] + prec_obs_sample[i-1,2]),
                               sd = 1 / sqrt(prec_add[2] + prec_obs_sample[i-1,2]))
  C_samples[730, 3, i] <- rnorm(n = 1, 
                               mean = (prec_add[3]*C_samples[729, 3, i] + prec_obs_sample[i-1,3]*C_obs[730, 3]) / (prec_add[3] + prec_obs_sample[i-1,3]),
                               sd = 1 / sqrt(prec_add[3] + prec_obs_sample[i-1,3]))
  C_samples[730, 4, i] <- rnorm(n = 1, 
                               mean = (prec_add[4]*C_samples[729, 4, i] + prec_obs_sample[i-1,4]*C_obs[730, 4]) / (prec_add[4] + prec_obs_sample[i-1,4]),
                               sd = 1 / sqrt(prec_add[4] + prec_obs_sample[i-1,4]))
  C_samples[730, 5, i] <- rnorm(n = 1, 
                               mean = (prec_add[5]*C_samples[729, 5, i] + prec_obs_sample[i-1,5]*C_obs[730, 5]) / (prec_add[5] + prec_obs_sample[i-1,5]),
                               sd = 1 / sqrt(prec_add[5] + prec_obs_sample[i-1,5]))
  
  if (i %% 100 == 0){
    cat(paste(i, 'iterations done \n'))
    #print(sd)
    print(colMeans(acceptance[(1:i),]))
  }
  for (j in 1:5){
    SSCobsi = sum((C_samples[,j,i] - C_obs[,j])^2)
    prec_obs_sample[i, j] = rgamma(1,374, .5*SSCobsi + 1)
  }
}
#})
time2 <- Sys.time()
time2 - time1

for (i in 1:11){
  #png(filename = paste0(paste0('./Plots/p',i),'.png'), )
  hist(p[burn:chain_length, i], breaks = 50,  xlab = paste0('p',i),
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

for (j in 1:5){
  #png(filename = paste0(paste0('./Plots/sd_obs',j),'.png'))
  hist(prec_obs_sample[burn:chain_length, j], breaks = 50, main = paste0('True Value is ', prec_obs[j]))
  #dev.off()
  #cat(paste('The mean is', mean(prec_obs_sample[burn:chain_length,j])))
  #cat('\n')
  readline('Press Enter')
}

png(filename = './Plots/pairs_postburn.png', width = 1100, height = 1100)
pairs(p[burn:chain_length, 1:10])
dev.off()

library(matrixStats)
plot(C[,1], type = 'l')
points(rowMeans(C_samples[,1,burn:chain_length]), type = 'l', col = 'red')
points(rowQuantiles(C_samples[,1,burn:chain_length], prob = .05), type = 'l', col = 'red')
points(rowQuantiles(C_samples[,1,burn:chain_length], prob = .95), type = 'l', col = 'red')


plot(C[,2], type = 'l')
points(rowMeans(C_samples[,2,burn:chain_length]), type = 'l', col = 'red')
points(rowQuantiles(C_samples[,2,burn:chain_length], prob = .05), type = 'l', col = 'red')
points(rowQuantiles(C_samples[,2,burn:chain_length], prob = .95), type = 'l', col = 'red')
#points(C_samples[,2,10000], type = 'l', col = 'red')

plot(C[,3], type = 'l')
points(rowMeans(C_samples[,3,burn:chain_length]), type = 'l', col = 'red')
points(rowQuantiles(C_samples[,3,burn:chain_length], prob = .05), type = 'l', col = 'red')
points(rowQuantiles(C_samples[,3,burn:chain_length], prob = .95), type = 'l', col = 'red')
#points(C_samples[,3,10000], type = 'l', col = 'red')

plot(C[,4], type = 'l')
points(rowMeans(C_samples[,4,burn:chain_length]), type = 'l', col = 'red')
points(rowQuantiles(C_samples[,4,burn:chain_length], prob = .05), type = 'l', col = 'red')
points(rowQuantiles(C_samples[,4,burn:chain_length], prob = .95), type = 'l', col = 'red')

plot(C[,5], type = 'l')
points(rowMeans(C_samples[,5,burn:chain_length]), type = 'l', col = 'red')
points(rowQuantiles(C_samples[,5,burn:chain_length], prob = .05), type = 'l', col = 'red')
points(rowQuantiles(C_samples[,5,burn:chain_length], prob = .95), type = 'l', col = 'red')

plot(Cf_e_nn, type = 'l')
points(rowMeans(C_samples[,1,burn:chain_length]), type = 'l', col = 'red')

plot(Cw_e_nn, type = 'l')
points(rowMeans(C_samples[,2,burn:chain_length]), type = 'l', col = 'red')
#points(C_samples[,2,10000], type = 'l', col = 'red')

plot(Cr_e_nn, type = 'l')
points(rowMeans(C_samples[,3,burn:chain_length]), type = 'l', col = 'red')
#points(C_samples[,3,10000], type = 'l', col = 'red')

plot(Clit_e_nn, type = 'l')
points(rowMeans(C_samples[,4,burn:chain_length]), type = 'l', col = 'red')

plot(Csom_e_nn, type = 'l')
points(rowMeans(C_samples[,5,burn:chain_length]), type = 'l', col = 'red')

## Comparison of fits via likelihood
buildCpred_true <- Cpred(p_e, C, LAI, maxt, mint, Ca, lat_e, yearday, Nit_e, rad, 730)
Cpred_true <- cbind(buildCpred_true$Cfpred, buildCpred_true$Cwpred, buildCpred_true$Crpred,
      buildCpred_true$Clitpred, buildCpred_true$Csompred)
LL1 <- cpp_ss_likelihood(C_obs, C, c(.25,16,2.25, .25, 9), Cpred_true, c(sd_cf^2,sd_cw^2,sd_cr^2, sd_clit^2, sd_csom^2), 
                  c(100,9200,100,20,11000), c(4,4,4,4,4), 730)

p_algo <- rep(0, 11)
for (i in 1:11){
  p_algo[i] = mean(p[burn:chain_length, i])
}
C_guess = cbind(rowMeans(C_samples[,1,]), rowMeans(C_samples[,2,]), rowMeans(C_samples[,3,]),
                rowMeans(C_samples[,4,]), rowMeans(C_samples[,5,]))
Cobs_guess = cbind(rowMeans(C_obs[,1,]), rowMeans(C_obs[,2,]), rowMeans(C_obs[,3,]),
                   rowMeans(C_obs[,4,]), rowMeans(C_obs[,5,]))
buildCpred_guess <- Cpred(p_algo, C_guess, LAI, maxt, mint, Ca, lat_e, yearday, Nit_e, rad, 730)
Cpred_guess <- cbind(buildCpred_guess$Cfpred, buildCpred_guess$Cwpred, buildCpred_guess$Crpred,
                    buildCpred_guess$Clitpred, buildCpred_guess$Csompred)
LL2 <- cpp_ss_likelihood(C_obs, C_guess, c(.25,16,2.25, .25, 9), Cpred_guess, c(sd_cf^2,sd_cw^2,sd_cr^2, sd_clit^2, sd_csom^2), 
                         c(100,9200,100,20,11000), c(4,4,4,4,4), 730)

## turn off adaptive part
## pairs plots
## sample these in log space

## set all parameters to be fixed 
## fix p_2 at .40889, try fitting
## try fitting just the latent states with ALL parameters fixed, see if there is a problem there
## multiple synthetic data sets?

## start with predicted GPP
## add g fluxes into the likelihood
## check how robust we are to process noise and observational noise
## robust to thinning data / missing data?

## try transforming variables at the end rather than inside the loop