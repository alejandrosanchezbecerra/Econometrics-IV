
sample_dataAR <- function(Phi,Sigma,p,m,n,...) {

  epsilon <- MASS::mvrnorm(n+p,mu=rep(0,m),Sigma=Sigma)
  data    <- matrix(0,n+p,m)

  for(t in (p+1):(p+n)) {

    data[t,] <- epsilon[t,]
    order    <- 1

    while(order <= p) {
      data[t,] <- data[t,]+matrix(c(data[t-order,]),1,m)%*%t(Phi[[order]])
      order    <- order+1
    }

  }

  return(data[(p+1):(p+n),])

}

sample_dataMA <- function(Theta,Sigma,q,m,n,...) {

  epsilon <- MASS::mvrnorm(n+q,mu=rep(0,m),Sigma=Sigma)
  data    <- matrix(0,n+q,m)

  for(t in (q+1):(n+q)) {

    data[t,] <- epsilon[t,]
    order    <- 1

    while(order <= q) {
      data[t,] <- data[t,]+matrix(c(epsilon[t-order,]),1,m)%*%t(Theta[[order]])
      order    <- order+1
    }

  }

  return(data[(q+1):(n+q),])

}

## TODO: Function that computes MSE given prediction and data.

## TODO: Function that computes criteria
