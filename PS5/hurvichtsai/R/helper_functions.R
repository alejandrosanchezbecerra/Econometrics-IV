
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
sample_data <- function(Phi,Sigma,p,m,n,NumLags,...) {

  epsilon <- MASS::mvrnorm(n+NumLags+1,mu=rep(0,m),Sigma=Sigma)
  data    <- matrix(0,n+NumLags+1,2)

  for(t in (NumLags+1):(n+NumLags+1)) {

    data[t,] <- epsilon[t,]
    order    <- 1

    while(order <= p) {
      data[t,] <- data[t,]+matrix(c(data[t-order,]),1,2)%*%t(Phi[[order]])
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

conc_lik_VAR <- function(Sigma_vec,beta,Y,X,p,m,n) {

  Sigma <- get_Sigma(Sigma_vec,m)
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

  # USE VAR as initial guess to speed up.
  model_VAR  <- vars::VAR(as.data.frame(data[(NumLags-p+1):(n+NumLags),]),
                          p = p,type="none")
  beta       <- t(vars::Bcoef(model_VAR))

  # Bmat <- matrix(0,2,2)
  # Bmat[lower.tri(Bmat,diag=TRUE)] <- NA
  # Sigma <- SVAR(model_VAR,estmethod="direct",Amat=diag(2),Bmat=Bmat)$Sigma.U/100

  # beta_guess <- cbind(model_VAR$varresult$V1$coefficients,
  #                   model_VAR$varresult$V2$coefficients)
  # Sigma_guess <- summary(model_VAR)$covres
  # init_guess <- get_param_vec(beta = beta_guess,Sigma = Sigma_guess )

#  init_guess <- get_param_vec(beta = matrix(0,m^2,p),
#                              Sigma = diag(m) )

  init_guess <- get_Sigma_vec(Sigma=summary(model_VAR)$covres)
  model     <- optim(par=init_guess,
                     fn=function(Sigma) {conc_lik_VAR(Sigma,beta,Y,X,p,m,n)},
                     method="L-BFGS-B")
  est_coef  <- model$par
  Sigma     <- get_Sigma(est_coef,m)

  # ls(summary(vars::VAR(as.data.frame(data[(NumLags-p+2):(n+NumLags+1),]),
  #          p = p,type="none")))

  # summary(vars::VAR(as.data.frame(data[(NumLags-p+1):(n+NumLags),]),
  #                  p = p,type="none"))$corres

  #model     <- optim(par=init_guess,
  #                   fn=function(par) {lik_VAR(par,Y,X,p,m,n)},
  #                   method="L-BFGS-B")
  #est_coef  <- model$par
  #param     <- get_parameters(est_coef,p,m)

  # grad(function(par) {lik_VAR(par,Y,X,p,m,n)},est_coef)

  # return(param)

  return(list(beta=beta,
              Sigma=Sigma))

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

get_Sigma <- function(Sigma_vec,m) {

  cholSigma <- matrix(0,m,m)
  cholSigma[lower.tri(cholSigma,diag=TRUE)] <-  Sigma_vec
  Sigma     <- cholSigma%*%t(cholSigma )

  return(Sigma)

}

get_Sigma_vec <- function(Sigma) {

  Sigma <- Sigma
  cholSigma <- t(chol(Sigma))
  cholSigma <- cholSigma[lower.tri(cholSigma,diag=TRUE)]

  return(cholSigma)

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

#' Title
#'
#' @param Phi
#' @param m
#' @param p
#'
#' @return
#' @export
#'
#' @examples
Phi_to_Beta <- function(Phi,m,p) {

  beta  <- matrix(0,0,m)
  order <- 1

  while ( order <= p ) {
    beta  <- rbind(beta,t(Phi[[order]]))
    order <- order+1
  }

  return(beta)
}

#' Title
#'
#' @param beta
#' @param m
#' @param p
#'
#' @return
#' @export
#'
#' @examples
Beta_to_Phi <- function(beta,m,p) {

  Phi <- list()

  order <- 1
  while ( order <= p ) {
    Phi[[order]] <- t(beta[(m*(order-1)+1):(order*m),1:m])
    order        <- order+1
  }

  return(Phi)
}


#' Title
#'
#' @param Phi
#' @param m
#' @param p
#' @param n
#'
#' @return
#' @export
#'
#' @examples
construct_H <- function(Phi,m,p,n) {

  H <- matrix(0,m*n,m*n)
  stacked_Phi <- matrix(NA,0,m)

  Phi

  for(i in n:1){

    # print(i)

    if( i == n) {
      stacked_Phi <- diag(m)
    } else {
      if( n-i <= p ) {
        stacked_Phi <- rbind(stacked_Phi,-Phi[[n-i]])
      } else {
        stacked_Phi <- rbind(stacked_Phi,matrix(0,m,m))
      }
    }

    H[((i-1)*m+1):(n*m),((i-1)*m+1):(i*m)] <- stacked_Phi

  }

  return(H)

}

# hola <- example2()
# hola$Phi0

# Phi_to_Beta(hola$Phi0,hola$m,hola$p)
