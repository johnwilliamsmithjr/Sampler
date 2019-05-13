In order to run the DALEC DEMCMC routine, open the file DALEC_par_DEMCMC.R. From there, run the lines that source 
the files "make_u.R" and "prepare_MCMC.R". The file "make_u.R" allows us to transform the bounded parameters to
an infinite dimensional space so that we can perform the DEMCMC algorithm from Ter Braak without having to worry if
the proposal distribution is symmetric. The file "prepare_MCMC.R" can be edited to change values for: the initial conditions,
process precisions, observation precisions, and others. This file will also load several libraries, including: laGP, lhs,
Rcpp, matrixStats, and others. C++ code that implements the process models is compiled here so that it may be called from R. 
R functions that aggregrate the data to call the likelihood function from C++ as well as the latent state updates are also
run here.
