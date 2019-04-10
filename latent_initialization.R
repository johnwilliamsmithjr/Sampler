## Function for initializations of latent states
## Arguments: 
## obs.ind - indices of when data was observed
## data - matrix or dataframe, with observed data for rows corresponding to obs.ind and NAs elsewhere
## min.ind - index where you want to begin estimating latent states
## max.ind - terminal index for estimating latent states
## type - method to use to estimate initial latent states where data is not observed
## currently supported methods: pwl (piecewise linear interpolation), gp (gaussian process)

latent_initialization <- function(data, obs.ind, min.ind, max.ind, type = rep('pwl', dim(data)[2])){
  ## check type argument here
  latent.init <- c()
  for (i in 1:dim(data)[2]){
    if (type[i] == 'pwl'){
      pwl.approx <- approxfun(obs.ind, data[obs.ind,i])
      latent.init <- cbind(latent.init, pwl.approx(min.ind:max.ind))
      rm(pwl.approx)
    }
    if (type[i] == 'gp'){
      da <- darg(list(mle = TRUE, max = 10), randomLHS(1000,1))
      gpi <- newGPsep(matrix(obs.ind/730), data[obs.ind, i], d = da$start, g = 1e-6, dK = TRUE)
      mle <- mleGPsep(gpi, param = 'd', tmin = da$min, tmax = da$max, ab = da$ab)
      p <- predGPsep(gpi, matrix((min.ind:max.ind)/max.ind, ncol = 1))
      latent.init <- cbind(latent.init, p$mean)
      deleteGPsep(gpi)
    }
  }
  return(latent.init)
}
