cpp_ss_likelihood = function(Cobs, C, Sigma, Cpred, sd, init_mean, init_sd, N){
  LL = 0
  for (i in 1:5){
    LL = LL + lgammad(sd[i], 10, 1)
  }
  LL = LL + SSLL(Cobs, C, Sigma, Cpred, sd, init_mean, init_sd, N)
  return(LL)
}