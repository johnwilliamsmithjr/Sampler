## Program: DALEC_ev_noisy.R
## Author: John W Smith Jr
## Date: 10/19/2018
## Description: Simulates data using the Reflex Evergreen model. Adds both process and observational noise into the data.
## After simulating the data, a STAN model is compiled (see model_fixG.stan for details). This model is then fit using STAN, 
## a Hamiltonian Monte Carlo sampler. Plots are created for densities of the underlying parameters with credible intervals, and
## plots of the latent carbon pools are compared to the true (deterministic / noiseless) states. 

## Set directory to be location of the folder containing the following files:
##    - DALEC_ev_noisy.R (this file)
##    - DALEC_ev_R.R
##    - model_fixG.stan
##    - EV1_drivers.csv
##    - EV1_initial.csv
direc = 'Reflex/Sampler'
setwd(paste('~/', direc, sep =''))

GPP = function(lai, p11, maxt, mint, Ca, lat, yearday, Nit, rad){
  ## Creates GPP function in order to compute the GPP for a day based on values input from the drivers
  
  ## Sets parameters
  psid = -2
  rtot = 1
  trange = .5*(maxt - mint)
  a = c(p11, 
        0.0156935, 
        4.22273,
        208.868,
        0.0453194, 
        0.37836,
        7.19298,
        0.011136,
        2.1001, 
        0.789798)
  gs = abs((psid))^a[10] / (a[6]*rtot + trange)
  pp = lai*Nit/gs*a[1]*exp(a[8]*maxt)
  qq = a[3] - a[4]
  ci = .5*(Ca + qq - pp + ((Ca+qq-pp)^2 - 
                             4*(Ca*qq - pp*a[3]))^(.5) )
  e0 = (a[7]*lai**2 / (lai**2 + a[9]) )
  dec = -23.4*cos((360*(yearday + 10)/365)*pi/180)*pi/180
  mult = tan(lat)*tan(dec)
  
  ## Changes dayl variable depending on value of mult
  if (mult >= 1){
    dayl = 24
  } else if(mult <= -1){
    dayl = 0
  } else{
    dayl = 24*acos(-mult)/pi
  }
  cps = e0*rad*gs*(Ca - ci) / (e0*rad + gs*(Ca - ci))
  gpp = cps*(a[2]*dayl + a[5])
  return(gpp)
}
## Sets global constants
nday = 730
psid = -2
rtot = 1

## Reads in drivers and sets column names
if ('EV1_drivers.csv' %in% dir()){
  drivers_ev = read.csv('./EV1_drivers.csv', header = F)
  colnames(drivers_ev) = c('projectday', 'mint', 'maxt', 'rad', 'ca', 'yearday')
} else{
  cat('Please add file EV1_drivers.csv to the directory. \n')
}

## Reads in initial conditions
if ('EV1_initial.csv' %in% dir()){
  init_ev = read.csv('./EV1_initial.csv', header = F)
} else{
  cat('Please add file EV1_initial.csv to the directory. \n')
}
## Sets values for p and initializes relevant variables
#p_e = init_ev[1:11,]
#p_e = rep(NA, 11)
#p_e[1] = 1e-6 + (.01 - 1e-6)*runif(1,0,1)
#p_e[2] = .2 + .5*runif(1,0,1)
#p_e[3] = .01 + .49*runif(1,0,1)
#p_e[4] = .01 + .49*runif(1,0,1)
#p_e[5] = 1e-4 + (.1 - 1e-4)*runif(1,0,1)
#p_e[6] = 1e-6 + (.01 - 1e-6)*runif(1,0,1)
#p_e[7] = 1e-4 + (.01 - 1e-4)*runif(1,0,1)
#p_e[8] = 1e-5 + (.1 - 1e-5)*runif(1,0,1)
#p_e[9] = 1e-6 + (.01 - 1e-6)*runif(1,0,1)
#p_e[10] = .05 + .15*runif(1,0,1)
#p_e[11] = 5 + 15*runif(1,0,1)
p_e = c(0.009000000, 0.408894241, 0.398498134, 0.015000000, 0.005620647,
        0.0001, 0.000150000, 0.073184997, 0.000010000, 0.190000000, 18.000000000)

Cf_e_nn = rep(NA, nday)
Cr_e_nn = rep(NA, nday)
Cw_e_nn = rep(NA, nday)
Clit_e_nn = rep(NA, nday)
Csom_e_nn = rep(NA, nday)
G_e = rep(NA, nday)
Ra_e = rep(NA, nday)
Af_e = rep(NA, nday)
Aw_e = rep(NA, nday)
Ar_e = rep(NA, nday)
Lf_e = rep(NA, nday)
Lw_e = rep(NA, nday)
Lr_e = rep(NA, nday)
Rh1_e = rep(NA, nday)
Rh2_e = rep(NA, nday)
D_e = rep(NA, nday)
NEE_e_nn = rep(NA, nday)
LAI = rep(NA, nday)
Cf_e_nn[1] = init_ev[12,]
Cr_e_nn[1] = init_ev[13,]
Cw_e_nn[1] = init_ev[14,]
Clit_e_nn[1] = init_ev[15,]
Csom_e_nn[1] = init_ev[16,]

## Scales lat, sets Nit and LMA
lat_e = init_ev[17,]
lat_e = lat_e * pi / 180
Nit_e = init_ev[18,]
LMA_e = init_ev[19,]

## Iterates through forloop to compute carbon fluxes
for (i in 1:nday){
  ## j is a dummy variable that was needed to match the results from the original F90 model, which did not use indices in the loop
  ## and instead wrote to the .csv file at every iteration in the loop
  if (i == 1){
    j = 1
  } else{
    j = i-1
  }
  ## Computes LAI, passes LAI, p11, and initial conditions and drivers into GPP function to calculate G_e
  LAI[i] = max(.1, Cf_e_nn[j]/LMA_e)
  G_e[i] = GPP(lai = LAI[i], p11 = p_e[11], maxt = drivers_ev$maxt[i], mint = drivers_ev$mint[i], Ca = drivers_ev$ca[i], 
               lat = lat_e, yearday = drivers_ev$yearday[i], Nit = Nit_e, rad = drivers_ev$rad[i])
  ## Calculates Trate
  Trate = .5*exp(p_e[10]*.5*(drivers_ev$maxt[i] + drivers_ev$mint[i]))
  ## Calculates values needed to compute carbon pools
  Ra_e[i] = p_e[2]*G_e[i]
  Af_e[i] = (G_e[i] - Ra_e[i])*p_e[3]
  Ar_e[i] = (G_e[i] - Ra_e[i] - Af_e[i])*p_e[4]
  Aw_e[i] = G_e[i] - Ra_e[i] - Af_e[i] - Ar_e[i]
  Lf_e[i] = p_e[5]*Cf_e_nn[j]
  Lw_e[i] = p_e[6]*Cw_e_nn[j]
  Lr_e[i] = p_e[7]*Cr_e_nn[j]
  Rh1_e[i] = p_e[8]*Clit_e_nn[j]*Trate
  Rh2_e[i] = p_e[9]*Csom_e_nn[j]*Trate
  D_e[i] = p_e[1]*Clit_e_nn[j]*Trate
  ## Compute Carbon Pools, add process error noise
  Cf_e_nn[i] = Cf_e_nn[j] + Af_e[i] - Lf_e[i] 
  Cw_e_nn[i] = Cw_e_nn[j] + Aw_e[i] - Lw_e[i]
  Cr_e_nn[i] = Cr_e_nn[j] + Ar_e[i] - Lr_e[i] 
  Clit_e_nn[i] = Clit_e_nn[j] + Lf_e[i] + Lr_e[i] - Rh1_e[i] - D_e[i] 
  Csom_e_nn[i] = Csom_e_nn[j] + D_e[i] - Rh2_e[i] + Lw_e[i] 
  NEE_e_nn[i] = Ra_e[i] + Rh1_e[i] + Rh2_e[i] - G_e[i]
}

plot(Cf_e_nn, type = 'l')

## Possible parameters!
## 0.001997662  0.408894241  0.398498134  0.264426926  0.005620647  0.009841923  0.006736792  0.073184997  0.007472635  0.137474255 11.977141728

# p_e
#[1]  0.001997662  0.408894241  0.398498134  0.015000000  0.005620647  0.000100000  0.006736792  0.073184997  0.007472635  0.137474255 11.977141728

# p_e
# [1]  0.001997662  0.408894241  0.398498134  0.015000000  0.005620647  0.000100000  0.000500000  0.073184997  0.007472635  0.137474255 11.977141728

# p_e - WORKS
# > p_e
# [1]  0.009000000  0.408894241  0.398498134  0.015000000  0.005620647  0.000100000  0.000150000  0.073184997  0.000010000  0.190000000 18.000000000