init_var <- c(1,4,.5,.5,4)
init_mean <- c(100,9200,100,20,11000)
#var_add <- c(2, 9, (2)^2, .25, 16)
var_add <- c(3, 16, .5^2, .25, 16)
#var_obs <- c(4, 16, 1.5^2, 1^2, 25)
var_obs <- c(4, 16, .5^2, 1^2, 25)
prec_add <- 1 / var_add
prec_obs <- 1 / var_obs
prec_init <- 1 / init_var
## Sources R scripts to generate data
source('./DALEC_data_gen.R')
#source('./DALEC_nn.R')

# Change directory back
setwd('~/Sampler')
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
library(lhs)
library(HDInterval)

## Compile C++ scripts from src folder 
sourceCpp('./src/lmvnd_ll.cpp')
sourceCpp('./src/likelihood_SS.cpp')
sourceCpp('./src/lgamma.cpp')
sourceCpp('./src/Cpred.cpp')
## Load required R functions
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
chain_length <- 8000
burn <- 2000
num_gen <- 5

## ----- Preparing MCMC ------ ##
# Compute LAI estimate to be fed into p_11 optimization routine
LAI <- pmax(.1, Cf.obs/LMA_e)
# Column bind matrices and build initial Cpred matrix
C_obs <- cbind(Cf.obs, Cw.obs, Cr.obs, Clit.obs, Csom.obs)
#missing_data_ind = seq(from = 1, to = 730, by = 3)
missing_data_ind = seq(from = 1, to = 730, by = 3)
#missing_data_ind = c(missing_data_ind, 730)
C_obs[-missing_data_ind,] <- NA
C <- cbind(Cf, Cw, Cr, Clit, Csom)

# Initalize matrix for p draws to be stored
p <- array(NA, dim = c(chain_length, 11, num_gen))
# Set vector of upper and lower bounds for accept reject sampling
plower <- c(1e-6, .2, .01, .01, 1e-4, 1e-6, 1e-6, 1e-5, 1e-6, .05, 5)
pupper <- c(.01, .7, .5, .5, .1, .01, .01, .1, .01, .2, 20)
# Fix p11 and initalize first guesses at the midpoint

#for (gen in 1:num_gen){
#  p[1, ,gen] <- runif(dim(p)[2], plower, pupper)
#}
#p[1,,] <- (t(randomLHS(num_gen,11)) * (pupper - plower) + plower)

#p[1,1:11] <- .5*(pupper[1:11] + plower[1:11])
Cpredi <- array(NA, dim = c(730, 5, num_gen))
Gpred <- array(NA, dim = c(730, num_gen))
for (gen in 1:num_gen){
  buildCpred <- Cpred(p[1,,gen], C, LAI, maxt, mint, Ca, lat_e, yearday, Nit_e, rad, 730)
  Gpred[,gen] = buildCpred$G
  Gpred[1,gen] <- G_obs[1]
  Cpredi[,,gen] <- cbind(buildCpred$Cfpred, buildCpred$Cwpred, buildCpred$Crpred,
                         buildCpred$Clitpred, buildCpred$Csompred)
  rm(buildCpred)
}

G <- array(NA, dim = c(chain_length, 730, num_gen))
G[1,,] <- G_obs

phi_g <- 1 / .025
tau_g <- 1 / .04
# create matrix for acceptance for p_i 
#acceptance <- matrix(0, nrow = chain_length, ncol = 11)
acceptance2 <- array(0, dim = c(chain_length, 11, num_gen))
#prec_Accept <- matrix(0, nrow = chain_length, ncol = 5)
#prec_Accept2 <- matrix(0, nrow = chain_length, ncol = num_gen)
prec_Accept2 <- array(0, dim = c(chain_length, 5, num_gen))
#prec_sd <- rep(.2, 5)

# Initialize arrays for C_samples
C_samples <- array(NA, dim = c(730, 5, chain_length, num_gen))
C.0 <- array(NA, dim = c(chain_length, 5, num_gen))
C.0[1,,] <- init_mean
# Initialize first entry to observed data
C_samples[, , 1,] <- latent_initialization(C_obs, missing_data_ind, 1, 730, type = c('gp', 'gp', 'pwl', 'pwl', 'gp'))

#G[,,]

prec_add_sample = array(1, dim = c(chain_length, 5, num_gen))

for (i in 1:num_gen){
  prec_add_sample[1,,i] <- runif(5, .1, 10)
}

p_cands <- randomLHS(1000, 11)
i = 2
LL <- rep(NA, 1000)
for (j in 1:1000){
  p_est <- p_cands[j,] * (pupper - plower) + plower
  buildCpredi <- Cpred(p_est, C_samples[, , i-1,gen], LAI_ = pmax(.1, C_samples[,1,i-1,gen]/LMA_e), maxt, mint, Ca, lat_e, yearday, Nit_e, rad, 730)
  Cpredi <- cbind(buildCpredi$Cfpred, buildCpredi$Cwpred, buildCpredi$Crpred,
                  buildCpredi$Clitpred, buildCpredi$Csompred)
  Gpred <- buildCpredi$G
  rm(buildCpredi)
  LL[j] <- cpp_ss_likelihood(C_obs, C_samples[, , i-1, gen], 1 / prec_obs, Cpredi, 1 / prec_add_sample[i-1,,gen], 
                             c(100,9200,100,20,11000), init_var, 730, G_pred = Gpred, G_e = G_obs, G = G[i-1,,gen], ind = missing_data_ind - 1)
}

for (gen in 1:num_gen){
  p[1,,gen] <- p_cands[order(LL, decreasing = TRUE)[gen], ] * (pupper - plower) + plower
}

Cpredi <- array(NA, dim = c(730, 5, num_gen))
Gpred <- array(NA, dim = c(730, num_gen))
for (gen in 1:num_gen){
  buildCpred <- Cpred(p[1,,gen], C, LAI, maxt, mint, Ca, lat_e, yearday, Nit_e, rad, 730)
  Gpred[,gen] = buildCpred$G
  Gpred[1,gen] <- G_obs[1]
  Cpredi[,,gen] <- cbind(buildCpred$Cfpred, buildCpred$Cwpred, buildCpred$Crpred,
                         buildCpred$Clitpred, buildCpred$Csompred)
  rm(buildCpred)
}