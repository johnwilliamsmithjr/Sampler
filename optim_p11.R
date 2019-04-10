optim_p11 = function(seq.length, G_vec){
p11_seq = seq(from = 5, to = 20, length.out = seq.length)
sum = rep(0, seq.length)
for (j in 1:seq.length){
  for (i in 2:730){
    sum[j] = sum[j] + (GPP(lai = LAI[i], p11 = p11_seq[j], maxt = drivers_ev$maxt[i], mint = drivers_ev$mint[i], Ca = drivers_ev$ca[i], 
        lat = lat_e, yearday = drivers_ev$yearday[i], Nit = Nit_e, rad = drivers_ev$rad[i]) - G_vec[i])^2
  }
}
estimate = p11_seq[which(sum == min(sum))]
return(estimate)
}
