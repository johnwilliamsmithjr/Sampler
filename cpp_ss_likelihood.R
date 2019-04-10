cpp_ss_likelihood = function(Cobs, C, Sigma, Cpred, sd, init_mean, init_sd, N, G, G_e, ind){
  LL = 0
  for (i in 1:5){
    LL = LL + lgammad(sd[i], 10, 1)
  }
  LL = LL + SSLL(Cobs, C, Sigma, Cpred, sd, init_mean, init_sd, N, G, G_e, rep(.025, length(G_e)), ind = ind)
  #LL = LL + dmvnorm(G, G_e, sigma = diag(.025, length(G_e)), log = T)
  return(LL)
}

cpp_ss_likelihood2 = function(Cobs, C, Sigma, Cpred, sd, init_mean, init_sd, N, G, G_e, ind){
  LL = 0
  #for (i in 1:5){
  #  LL = LL + lgammad(sd[i], 10, 1)
  #}
  LL = LL + SSLL(Cobs, C, Sigma, Cpred, sd, init_mean, init_sd, N, G, G_e, rep(.025, length(G_e)), ind = ind)
  #LL = LL + dmvnorm(G, G_e, sigma = diag(.025, length(G_e)), log = T)
  return(LL)
}
