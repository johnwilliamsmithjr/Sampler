# Changes working directory and sources R script to generate DALEC data
setwd('~/Reflex/Sampler')
source('./DALEC_data_gen.R')
# Change directory back
setwd('~/Reflex/Sampler')
# Load c++ functions through Rcpp
if (!("Rcpp" %in% installed.packages())) {
  install.packages('Rcpp')
}
library(Rcpp)
sourceCpp('./src/likelihood_SS.cpp')
sourceCpp('./src/lgamma.cpp')
sourceCpp('./src/debugcpred2.cpp')
source('cpp_ss_likelihood.R')

## ------------ Sampler Test Run -----------##
# Set carbon stock variables from data generation scripts
Cf <- Cf_e
Cw <- Cw_e
Cr <- Cr_e
Clit <- Clit_e
Csom <- Csom_e
# Set standard deviations 
sd_cf <- .125^2
sd_cw <- 1
sd_cr <- .25^2
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
chain_length <- 10000
burn <- 2000
# Initalize matrix for p draws to be stored
p <- matrix(NA, nrow = chain_length, ncol = 11)
# Set vector of upper and lower bounds for accept reject sampling
plower <- c(10e-6, .2, .01, .01, 10e-4, 10e-6, 10e-6, 10e-5, 10e-6, .05, 5)
pupper <- c(.01, .7, .5, .5, .1, .01, .01, .1, .01, .2, 20)
# Set inital guesses for p
p[1,] <- (plower + pupper) / 2

## ----- Preparing MCMC ------ ##
LAI <- pmax(.1, Cf/LMA_e)
C_obs <- cbind(Cf.obs, Cw.obs, Cr.obs, Clit.obs, Csom.obs)
C <- cbind(Cf, Cw, Cr, Clit, Csom)
buildCpred <- Cpred(p[1,], C, LAI, maxt, mint, Ca, lat_e, yearday, Nit_e, rad, 730)
Cpredi <- cbind(buildCpred$Cfpred, buildCpred$Cwpred, buildCpred$Crpred,
               buildCpred$Clitpred, buildCpred$Csompred)
rm(buildCpred)

## ----- Performing MCMC ------ ##
## Accept Reject Normal for p_i ##
#chain_length = 30000
LAI <- pmax(.1, Cf/LMA_e)
C_obs <- cbind(Cf.obs, Cw.obs, Cr.obs, Clit.obs, Csom.obs)
C <- cbind(Cf, Cw, Cr, Clit, Csom)
buildCpred <- Cpred(p[1,], C, LAI, maxt, mint, Ca, lat_e, yearday, Nit_e, rad, 730)
G = buildCpred$G
Cpredi <- cbind(buildCpred$Cfpred, buildCpred$Cwpred, buildCpred$Crpred,
               buildCpred$Clitpred, buildCpred$Csompred)
rm(buildCpred)

# Create vector for initial condition mean + stdev, as well as variances and precisions for obs and add
init_var <- c(4,4,4,4,4)
init_mean <- c(100,9200,100,20,11000)
var_add <- c(.125^2, 1, .25^2, .25, 1)
var_obs <- c(.25, 16, 1.5^2, .5^2, 9)
prec_add <- 1 / var_add
prec_obs <- 1 / var_obs

# create matrix for acceptance for p_i 
acceptance <- matrix(0, nrow = chain_length, ncol = 11)

# set proposal standard deviations
sd <- (pupper - plower) / 30
sd[9:10] <- 5*sd[2:11]
# Set point to start adapting proposal SD
start.adapt = 500

# Initialize arrays for C_samples and Cobs_samples
C_samples <- array(NA, dim = c(730, 5, chain_length))
Cobs_samples <- array(NA, dim = c(730, 5, chain_length))
# Initialize first entry to true values (for test runs, will not be done in practice)
Cobs_samples[, , 1] <- C_obs
C_samples[, , 1] <- C

for (i in 2:chain_length){
  ## Metropolis Step for p_i ##
  #p[i,] <- p[i-1,]
  p[i,] <- p_e
  buildCpredi <- Cpred(p[i,], C_samples[, , i-1], LAI, maxt, mint, Ca, lat_e, yearday, Nit_e, rad, 730)
  Cpredi <- cbind(buildCpredi$Cfpred, buildCpredi$Cwpred, buildCpredi$Crpred,
                 buildCpredi$Clitpred, buildCpredi$Csompred)
  Gpred <- buildCpredi$G
  rm(buildCpredi)
  for (j in 1:11){
    p_star <- p[i,]
    tmp <- 0
    while (tmp < plower[j] | tmp > pupper[j]){
      tmp <- rnorm(1, p[i,j], sd[j])
    }
    p_star[j] <- tmp
    buildCpredstar <- Cpred(p_star, C_samples[, , i-1], LAI, maxt, mint, Ca, lat_e, yearday, Nit_e, rad, 730)
    Cpred_star <- cbind(buildCpredstar$Cfpred, buildCpredstar$Cwpred, buildCpredstar$Crpred,
                       buildCpredstar$Clitpred, buildCpredstar$Csompred)
    rm(buildCpredstar)
    
    threshold <- cpp_ss_likelihood(Cobs_samples[, , i-1], C_samples[, , i-1], c(.25,16,2.25, .25, 9), Cpred_star, c(sd_cf^2,sd_cw^2,sd_cr^2, sd_clit^2, sd_csom^2), 
                                  c(100,9200,100,20,11000), c(4,4,4,4,4), 730) - 
      cpp_ss_likelihood(Cobs_samples[, , i-1], C_samples[, , i-1], c(.25,16,2.25, .25, 9), Cpredi, c(sd_cf^2,sd_cw^2,sd_cr^2, sd_clit^2, sd_csom^2), 
                        c(100,9200,100,20,11000), c(4,4,4,4,4), 730)
    if (log(runif(1,0,1)) < threshold){
      p[i,j] <- p_star[j]
      Cpredi <- Cpred_star
      G <- Gpred
      acceptance[i,j] <- 1
    }
    #if (i%%10 == 0 & i > start.adapt & mean(acceptance[(burn:i), j]) < .25 ){
    #  sd[j] <- sd[j] / 1.02
    #} else if (i%%50 == 0 & i > start.adapt & mean(acceptance[(burn:i), j]) > .35){
    #  sd[j] <- sd[j] * 1.02
    #}
  }
  ## DOUBLE CHECK THESE GIBBS FULL CONDITIONALS
  ## Gibbs Step for C_i for the initial state
  C_samples[,,i] <- C_samples[,,i-1]
  C_samples[1,1,i] <- rnorm(1, init_mean[1], sqrt(init_var[1]))
  C_samples[1,2,i] <- rnorm(1, init_mean[2], sqrt(init_var[2]))
  C_samples[1,3,i] <- rnorm(1, init_mean[3], sqrt(init_var[3]))
  C_samples[1,4,i] <- rnorm(1, init_mean[4], sqrt(init_var[4]))
  C_samples[1,5,i] <- rnorm(1, init_mean[5], sqrt(init_var[5]))
  ## Gibbs Step for Cobs_i for the inital state
  Cobs_samples[,,i] <- Cobs_samples[,,i-1]
  Cobs_samples[1,1,i] <- rnorm(1, C_samples[1,1,i], sqrt(var_obs[1]))
  Cobs_samples[1,2,i] <- rnorm(1, C_samples[1,2,i], sqrt(var_obs[2]))
  Cobs_samples[1,3,i] <- rnorm(1, C_samples[1,3,i], sqrt(var_obs[3]))
  Cobs_samples[1,4,i] <- rnorm(1, C_samples[1,4,i], sqrt(var_obs[4]))
  Cobs_samples[1,5,i] <- rnorm(1, C_samples[1,5,i], sqrt(var_obs[5]))
  
  for (k in 2:729){
    ## Gibbs steps for "middle" latent states
    
    a_t <- 1 - p[i,5]
    a_tp <- 1 - p[i,5]
    b_t <- G[k]*(1-p[i,2])*p[i,3]
    b_tp <- G[k+1]*(1-p[i,2])*p[i,3]
    
    C_samples[k,1,i] <- rnorm(n = 1, 
                      mean = (prec_add[1]*(a_t*C_samples[k-1,1,i] +b_t ) + prec_add[1]*(a_tp*C_samples[k+1,1,i]-a_tp*b_tp) 
                      + prec_obs[1]*Cobs_samples[k,1,i-1]) / (prec_add[1]*(1 + (a_tp)^2 ) + prec_obs[1]),
                      sd = 1 / sqrt((prec_add[1]*(1 + (a_tp)^2 ) + prec_obs[1])))
    
    a_t <- 1 - p[i,6]
    a_tp <- 1 - p[i,6]
    b_t <- G[k]*(1 - p[i,2])*(1 - p[i,3])*(1 - p[i,4])
    b_tp <- G[k+1]*(1 - p[i,2])*(1 - p[i,3])*(1 - p[i,4])
    
    C_samples[k,2,i] <- rnorm(n = 1, 
                      mean = (prec_add[2]*(a_t*C_samples[k-1,2,i] +b_t ) + prec_add[2]*(a_tp*C_samples[k+1,2,i]-a_tp*b_tp) 
                      + prec_obs[2]*Cobs_samples[k,2,i-1]) / (prec_add[2]*(1 + (a_tp)^2 ) + prec_obs[2]),
                      sd = 1 / sqrt((prec_add[2]*(1 + (a_tp)^2 ) + prec_obs[2])))
    
    a_t <- 1-p[i,7]
    a_tp <- 1-p[i,7]
    b_t <- G[k]*(1-p[i,2])*(1-p[i,3])*p[i,4]
    b_tp <- G[k+1]*(1-p[i,2])*(1-p[i,3])*p[i,4]
    
    C_samples[k,3,i] <- rnorm(n = 1, 
                      mean = (prec_add[3]*(a_t*C_samples[k-1,3,i] +b_t ) + prec_add[3]*(a_tp*C_samples[k+1,3,i]-a_tp*b_tp) 
                      + prec_obs[3]*Cobs_samples[k,3,i-1]) / (prec_add[3]*(1 + (a_tp)^2 ) + prec_obs[3]),
                      sd = 1 / sqrt((prec_add[3]*(1 + (a_tp)^2 ) + prec_obs[3])))
    
    a_t <- (1 - .5*exp(.5*p[i,10]*(mint[k] + maxt[k]))*(p[i,8]+p[i,1]))
    a_tp <- (1 - .5*exp(.5*p[i,10]*(mint[k+1] + maxt[k+1]))*(p[i,8]+p[i,1]))
      
    b_t <- p[i,5]*C_samples[k-1,1,i] + p[i,7]*C_samples[k-1,3,i]
    b_tp <- p[i,5]*C_samples[k,1,i] + p[i,7]*C_samples[k,3,i]
    
    C_samples[k,4,i] <- rnorm(n = 1, 
                      mean = (prec_add[4]*(a_t*C_samples[k-1,4,i] +b_t ) + prec_add[4]*(a_tp*C_samples[k+1,4,i]-a_tp*b_tp) 
                      + prec_obs[4]*Cobs_samples[k,4,i-1]) / (prec_add[4]*(1 + (a_tp)^2 ) + prec_obs[4]),
                      sd = 1 / sqrt((prec_add[4]*(1 + (a_tp)^2 ) + prec_obs[4])))
    
    a_t <- 1 - .5*p[i,9]*exp(.5*p[i,10]*(maxt[k] + mint[k]))
    a_tp <- 1 - .5*p[i,9]*exp(.5*p[i,10]*(maxt[k+1] + mint[k+1]))
    
    b_t <- p[i,6]*C_samples[k-1,2,i] + .5*p[i,1]*exp(.5*p[i,10]*(mint[k] + maxt[k]))*C_samples[k-1,4,i]
    b_tp <- p[i,6]*C_samples[k,2,i] + .5*p[i,1]*exp(.5*p[i,10]*(mint[k+1] + maxt[k+1]))*C_samples[k,4,i]
    
    
    C_samples[k,5,i] <- rnorm(n = 1, 
                      mean = (prec_add[5]*(a_t*C_samples[k-1,5,i] +b_t ) + prec_add[5]*(a_tp*C_samples[k+1,5,i]-a_tp*b_tp) 
                      + prec_obs[5]*Cobs_samples[k,5,i-1]) / (prec_add[5]*(1 + (a_tp)^2 ) + prec_obs[5]),
                      sd = 1 / sqrt((prec_add[5]*(1 + (a_tp)^2 ) + prec_obs[5])))
    
    
    Cobs_samples[k,1,i] <- rnorm(1, C_samples[k,1,i], sqrt(var_obs[1]))
    Cobs_samples[k,2,i] <- rnorm(1, C_samples[k,2,i], sqrt(var_obs[2]))
    Cobs_samples[k,3,i] <- rnorm(1, C_samples[k,3,i], sqrt(var_obs[3]))
    Cobs_samples[k,4,i] <- rnorm(1, C_samples[k,4,i], sqrt(var_obs[4]))
    Cobs_samples[k,5,i] <- rnorm(1, C_samples[k,5,i], sqrt(var_obs[5]))
  }
  ## Gibbs steps for final state
  C_samples[730, 1, i] <- rnorm(n = 1, 
                               mean = (prec_add[1]*C_samples[729, 1, i] + prec_obs[1]*Cobs_samples[730, 1, i-1]) / (prec_add[1] + prec_obs[1]),
                               sd = 1 / sqrt(prec_add[1] + prec_obs[1]))
  C_samples[730, 2, i] <- rnorm(n = 1, 
                               mean = (prec_add[2]*C_samples[729, 2, i] + prec_obs[2]*Cobs_samples[730, 2, i-1]) / (prec_add[2] + prec_obs[2]),
                               sd = 1 / sqrt(prec_add[2] + prec_obs[2]))
  C_samples[730, 3, i] <- rnorm(n = 1, 
                               mean = (prec_add[3]*C_samples[729, 3, i] + prec_obs[3]*Cobs_samples[730, 3, i-1]) / (prec_add[3] + prec_obs[3]),
                               sd = 1 / sqrt(prec_add[3] + prec_obs[3]))
  C_samples[730, 4, i] <- rnorm(n = 1, 
                               mean = (prec_add[4]*C_samples[729, 4, i] + prec_obs[4]*Cobs_samples[730, 4, i-1]) / (prec_add[4] + prec_obs[4]),
                               sd = 1 / sqrt(prec_add[4] + prec_obs[4]))
  C_samples[730, 5, i] <- rnorm(n = 1, 
                               mean = (prec_add[5]*C_samples[729, 5, i] + prec_obs[5]*Cobs_samples[730, 5, i-1]) / (prec_add[5] + prec_obs[5]),
                               sd = 1 / sqrt(prec_add[5] + prec_obs[5]))
  Cobs_samples[730,1,i] <- rnorm(1, C_samples[730,1,i], sqrt(var_obs[1]))
  Cobs_samples[730,2,i] <- rnorm(1, C_samples[730,2,i], sqrt(var_obs[2]))
  Cobs_samples[730,3,i] <- rnorm(1, C_samples[730,3,i], sqrt(var_obs[3]))
  Cobs_samples[730,4,i] <- rnorm(1, C_samples[730,4,i], sqrt(var_obs[4]))
  Cobs_samples[730,5,i] <- rnorm(1, C_samples[730,5,i], sqrt(var_obs[5]))
  
  if (i %% 100 == 0){
    cat(paste(i, 'iterations done \n'))
    #print(sd)
    print(mean(acceptance[(1:i),]))
  }
}

for (i in 1:11){
  hist(p[burn:chain_length, i], breaks = 50)
  cat(paste('The mean is', mean(p[burn:chain_length, i])))
  cat('\n')
  cat(paste('The true value is', p_e[i]))
  cat('\n')
  readline('Press enter to continue')
  plot(p[burn:chain_length, i], type = 'l')
  readline('Press enter to continue')
}
plot(C[,1], type = 'l')
points(rowMeans(C_samples[,1,]), type = 'l', col = 'red')

plot(C[,2], type = 'l')
points(rowMeans(C_samples[,2,]), type = 'l', col = 'red')
points(C_samples[,2,10000], type = 'l', col = 'red')

plot(C[,3], type = 'l')
points(rowMeans(C_samples[,3,]), type = 'l', col = 'red')
points(C_samples[,3,10000], type = 'l', col = 'red')

plot(C[,4], type = 'l')
points(rowMeans(C_samples[,4,]), type = 'l', col = 'red')

plot(C[,5], type = 'l')
points(rowMeans(C_samples[,5,]), type = 'l', col = 'red')

