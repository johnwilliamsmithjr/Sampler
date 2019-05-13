make_u = function(p, plower, pupper){
  if (is.null(dim(p))){
    n <- length(p)
    for (i in 1:n){
      p[i] <- (p[i] - plower[i]) / (pupper[i] - plower[i])
    }
  } else{
    n <- ncol(p)
    for (i in 1:n){
      p[,i] <- (p[,i] - plower[i]) / (pupper[i] - plower[i])
    }
  }
  return(p)
}
