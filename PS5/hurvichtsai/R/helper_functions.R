
#' Generate a random sample of an VAR(p) m-variate time series.
#'
#' @param Phi Transition Matrix
#' @param Sigma Covariance of error terms
#' @param p Order of VAR(p)
#' @param m Number of dimensions
#' @param n Number of observations
#' @param NumLags Number of pre-data obs set to zero.
#'
#' @return mx(n+NumLags) matrix
#' @export
#'
#' @examples
sample_data <- function(Phi,Sigma,p,m,n,NumLags) {

  epsilon <- MASS::mvrnorm(n+NumLags+1,mu=rep(0,m),Sigma=Sigma)
  data    <- matrix(0,n+NumLags+1,2)

  for(t in (NumLags+1):(n+NumLags+1)) {

    data[t,] <- epsilon[t,]
    order    <- 1

    while(order <= p) {
      data[t,] <- data[t,]+matrix(c(data[t-order,]),1,2)%*%Phi[[order]]
      order    <- order+1
    }

  }

  return(data)

}

#' Negative log-likelihood
#'
#' @param param Vector of parameters
#' @param Y m-variate time series
#' @param X lagged matrix of order (p). Does not include constant.
#' @param p Order of VAR(p)
#' @param m Dimension of time series
#'
#' @return Negative log-likelihood up to a constant
#' @export
#'
#' @examples
lik_VAR <- function(param,Y,X,p,m,n) {

  coeff <- get_parameters(param,p,m)
  beta  <- coeff$beta
  Sigma <- coeff$Sigma
  lik   <- n*log(det(Sigma))+
           psych::tr((Y-X%*%beta)%*%solve(Sigma)%*%t(Y-X%*%beta))

  return(lik)

}


#' Obtain maximum likelihood estimate
#'
#' @param data mx(n+NumLags) matrix
#' @param m Number of dimensions of time series
#' @param p Order of VAR(p)
#' @param NumLags Maximum number of lags set to zero in pre-data observations.
#'
#' @return
#' @export
#'
#' @examples
mle_VAR <- function(data,m,p,n,NumLags) {

  Y <- data[(NumLags+1):(n+NumLags),]
  X <- get_regressorMat(data,NumLags+1,n,p)

  init_guess <- get_param_vec(beta = matrix(0,m^2,p),
                              Sigma = diag(m) )

  model     <- optim(par=init_guess, fn=function(par) {lik_VAR(par,Y,X,p,m,n)})
  est_coef  <- model$par

  # model     <- nlminb(start=init_guess,
  #                     objective=function(par) {lik_VAR(par,Y,X,p,m,n)})
  # est_coef  <- model$par
  param     <- get_parameters(est_coef,p,m)

  return(param)

}

#' Construct stacked regressor matrix of lags of order p
#'
#' @param Z  Original time series
#' @param n0 Initial observation
#' @param n  Number of observations
#' @param p  Order of VAR(p)
#'
#' @return
#' @export
#'
#' @examples
get_regressorMat <- function(Z,n0,n,p) {

  if(p == 1) {
    X <- Z[(n0-1):((n0-1)+(n-1)),]

  } else if (p >= 1) {

    X <- Z[(n0-1):((n0-1)+(n-1)),]

    for( j in 2:p) {
      X <- cbind(X,Z[(n0-j):((n0-j)+(n-1)),])
    }
  }

  return(X)

}

#' Convert a vector of parameters into beta and Sigma matrices
#'
#' @param param First $p*(m^2)$ parameters are beta. The rest
#' are Sigma. The parameters are the cells of the Cholesky decomposition.-
#' @param p Order of VAR(p)
#' @param m Number of dimensions of time series
#'
#' @return
#' @export
#'
#' @examples
get_parameters <- function(param,p,m) {

  beta      <- matrix(param[1:(p*(m^2))],m*p,m)
  cholSigma <- matrix(0,m,m)
  cholSigma[lower.tri(cholSigma,diag=TRUE)] <-  param[(p*(m^2)+1):length(param)]
  Sigma     <- cholSigma%*%t(cholSigma )

  return(list(beta = beta,Sigma=Sigma))

}

#' Convert beta and Sigma matrices into a vector of parameters
#'
#' @param beta  Stacked matrix of transition matrices (Phi)
#' @param Sigma Covariance matrix of error terms
#'
#' @return Vector of parameters. The cells for Sigma correspond to the Cholesky
#' decomposition.
#' @export
#'
#' @examples
get_param_vec <- function(beta,Sigma) {

  Sigma <- Sigma
  cholSigma <- t(chol(Sigma))
  cholSigma <- cholSigma[lower.tri(cholSigma,diag=TRUE)]

  return(c(as.vector(beta),cholSigma))

}
