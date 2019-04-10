## Function for generating proposals in block using DEMCMC
## Arguments:
## i - current iteration of chain
## theta - parameter set to be sampled in block
## lower - vector of lower bounds for variables being updated in block
## upper - vector of upper bounds for variables being updated in block
## jitter.max - maximum jitter to be added to proposals
## gamma - tuneable scale factor used in DEMCMC algorithm 

DEMCMC_proposal <- function(i, theta, lower = NULL, upper = NULL, jitter.max, gamma = .5){
  if (is.null(lower) & is.null(upper)){
    m1 <- sample((1:(i-1)), 1)
    m2 <- sample((1:(i-1)), 1)
    theta.star <- theta[i-1,] + gamma*(theta[m1,] - theta[m2,]) + runif(dim(theta)[2], -jitter.max, jitter.max)
  }
  if (!(is.null(lower)) & !(is.null(upper))){
    tmp <- Inf
    while (!(all(tmp > (lower)) & all(tmp < (upper)))){
      m1 <- sample(1:(i-1), 1)
      m2 <- sample(1:(i-1), 1)
      tmp <- theta[i-1,] + gamma*(theta[m1,]-theta[m2,]) + runif(dim(theta)[2], min = -1e-6, max = 1e-6)
    }
    theta.star <- tmp
  }
  if (is.null(lower) & !(is.null(upper))){
    tmp <- Inf
    while(!(all(tmp < upper))){
      m1 <- sample(1:(i-1), 1)
      m2 <- sample(1:(i-1), 1)
      tmp <- theta[i-1,] + gamma*(theta[m1,]-theta[m2,]) + runif(dim(theta)[2], min = -1e-6, max = 1e-6)
    }
    theta.star <- tmp
  }
  if (is.null(upper) & !(is.null(lower))){
    tmp <- -Inf
    while(!(all(tmp > lower))){
      m1 <- sample(1:(i-1), 1)
      m2 <- sample(1:(i-1), 1)
      tmp <- theta[i-1,] + gamma*(theta[m1,]-theta[m2,]) + runif(dim(theta)[2], min = -1e-6, max = 1e-6)
    }
    theta.star <- tmp
  }
  return(theta.star)
}
