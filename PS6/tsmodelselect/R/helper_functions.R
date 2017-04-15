
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
sample_dataAR <- function(Phi0,Sigma0,p0,m,n,...) {

  epsilon <- MASS::mvrnorm(n+p0,mu=rep(0,m),Sigma=Sigma0)
  data    <- matrix(0,n+p0,m)

  for(t in (p0+1):(p0+n)) {

    data[t,] <- epsilon[t,]

    order    <- 1
    while(order <= p0) {
      data[t,] <- data[t,]+matrix(c(data[t-order,]),1,m)%*%t(Phi0[[order]])
      order    <- order+1
    }

  }

  return(data[(p0+1):(p0+n),])

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
sample_dataMA <- function(Theta0,Sigma0,q0,m,n,...) {

  epsilon <- MASS::mvrnorm(n+q0,mu=rep(0,m),Sigma=Sigma0)
  data    <- matrix(0,n+q0,m)

  for(t in (q0+1):(n+q)) {

    data[t,] <- epsilon[t,]
    order    <- 1

    while(order <= q0) {
      data[t,] <- data[t,]+matrix(c(epsilon[t-order,]),1,m)%*%t(Theta0[[order]])
      order    <- order+1
    }

  }

  return(data[(q0+1):(n+q0),])

}


#' Sample Data from an ARMA model
#'
#' @param Theta0
#' @param Sigma0
#' @param q
#' @param m
#' @param n
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
sample_dataARMA <- function(Phi0,Theta0,Sigma0,p0,q0,m,n,...) {

  epsilon <- MASS::mvrnorm(n+max(p0,q0),mu=rep(0,m),Sigma=Sigma0)
  data    <- matrix(0,n+max(p0,q0),m)
  data[1:max(p0,q0),] <- epsilon[1:max(p0,q0),]

  for(t in (max(p0,q0)+1):(n+max(p0,q0))) {

    data[t,] <- epsilon[t,]

    order    <- 1
    while(order <= q0 & (q0 != 0)) {
      data[t,]  <- data[t,]+matrix(c(epsilon[t-order,]),1,m)%*%t(Theta0[[order]])
      order     <- order + 1
    }

    order    <- 1
    while(order <= p0 & (p0 != 0)) {

      data[t,]  <- data[t,]+matrix(c(data[t-order,]),1,m)%*%t(Phi0[[order]])
      order     <- order +1
    }

  }

  return(data[(max(p0,q0)+1):(n+max(p0,q0)),])

}

#' Title
#'
#' @param X
#' @param Y
#' @param h
#' @param v
#' @param t
#'
#' @return
#' @export
#'
#' @examples
construct_training <- function(X,Y,h,v,t) {

  training_data <- list()

  if(t <= h+v) {
    training_data$X <- X[(h+v+t+1):length(Y),]
    training_data$Y <- Y[(h+v+t+1):length(Y),]
  } else if (t > h+v & t < dim(Y)[1]-h-v ) {
    training_data$X <- X[-((t-h-v):(t+h+v)),]
    training_data$Y <- Y[-((t-h-v):(t+h+v)),]
  } else {
    training_data$X <- X[1:(t-h-v-1),]
    training_data$Y <- Y[1:(t-h-v-1),]
  }

  return(training_data)

}

construct_test   <- function(X,Y,h,v,t) {

  test_data     <- list()
  test_data$X   <- X[(t-v):(t+v),]
  test_data$Y   <- Y[(t-v):(t+v),]

  return(test_data)
}

outofsample_fit <- function(X,Y,h,v,t,p) {

  # print(t)

  training_data <- construct_training(X,Y,h,v,t)
  test_data     <- construct_test(X,Y,h,v,t)

  beta_train_model   <- do.call(function(X,Y)
                                solve(t(as.matrix(X))%*%as.matrix(X),
                                      t(as.matrix(X))%*%as.matrix(Y)),
                                training_data)

  pred_error         <- test_data$Y - as.matrix(test_data$X) %*% beta_train_model

  # train_model   <- arima(training_data,order=c(p,0,0),include.mean=FALSE)
  # test_model    <- arima(training_data,order=c(p,0,0),include.mean=FALSE,
  #                       fixed = train_model$coef)

  MSE_t         <- t(pred_error)%*%(pred_error)/(length(test_data))

  return(MSE_t)

}


#' Title
#'
#' @param data
#' @param h
#' @param v
#' @param p
#'
#' @return
#' @export
#'
#' @examples
cv_AR <- function(data,h,v,p) {

  Y <- data[(p+1):length(data)]
  X <- matrix(NA,length(data)-p,0)

  for(lag in 1:p) {
     X <- cbind(X,data[(p+1-lag):(length(data)-lag)])
  }

  Y <- as.data.frame(Y)
  X <- as.data.frame(X)

  MSE_vec <- sapply((v):(length(data)-v-p),
                    function(t) outofsample_fit(X,Y,h,v,t,p))

  MSE     <- mean(MSE_vec)

  return(MSE)

}
