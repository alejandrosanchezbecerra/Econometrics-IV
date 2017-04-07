
#' Title
#'
#' @param Phi
#' @param Sigma
#' @param p
#' @param m
#' @param n
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
sample_dataAR <- function(Phi0,Sigma0,p,m,n,...) {

  epsilon <- MASS::mvrnorm(n+p,mu=rep(0,m),Sigma=Sigma0)
  data    <- matrix(0,n+p,m)

  for(t in (p+1):(p+n)) {

    data[t,] <- epsilon[t,]
    order    <- 1

    while(order <= p) {
      data[t,] <- data[t,]+matrix(c(data[t-order,]),1,m)%*%t(Phi0[[order]])
      order    <- order+1
    }

  }

  return(data[(p+1):(p+n),])

}

#' Title
#'
#' @param Theta
#' @param Sigma
#' @param q
#' @param m
#' @param n
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
sample_dataMA <- function(Theta0,Sigma0,q,m,n,...) {

  epsilon <- MASS::mvrnorm(n+q,mu=rep(0,m),Sigma=Sigma0)
  data    <- matrix(0,n+q,m)

  for(t in (q+1):(n+q)) {

    data[t,] <- epsilon[t,]
    order    <- 1

    while(order <= q) {
      data[t,] <- data[t,]+matrix(c(epsilon[t-order,]),1,m)%*%t(Theta0[[order]])
      order    <- order+1
    }

  }

  return(data[(q+1):(n+q),])

}

## TODO: Function that computes MSE given prediction and data.

## TODO: Function that computes criteria
