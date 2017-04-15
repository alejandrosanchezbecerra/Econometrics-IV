#' Title
#'
#' @param Phi0
#' @param Sigma0
#' @param n
#' @param Nsim
#' @param m
#' @param p0
#' @param maxLags
#'
#' @return
#' @export
#'
#' @examples
test_criteria_AR <- function(Phi0,Sigma0,n,Nsim,m,p0,maxLags) {

  results              <- array(NA,dim=c(Nsim,maxLags,3))
  error_prediction     <- matrix(NA,Nsim,maxLags)

  set.seed(12345)

  for(s in 1:Nsim) {

    print(paste('Simulation:',s,'of',Nsim))
    data <- sample_dataAR(Phi=Phi0,Sigma=Sigma0,p=p0,m=m,n=n+1)

    for(p in 1:maxLags) {

      model_AR  <- arima(as.data.frame(data[(maxLags-p+1):(n)]),order=c(p,0,0),
                         include.mean=FALSE)

      Sigma    <- (1/(n-maxLags-p))*t(model_AR$residuals)%*%model_AR$residuals

      results[s,p,1:2] <- unlist(criteria_func(Sigma=Sigma,n=n,p=p))
      # results[s,p,3] <- cross_validation_method(data,p)

      error_prediction[s,p] <- as.numeric(data[n+1]-
                                          predict(model_AR,n.ahead=1)$pred)

    }
  }

  table_data     <- select_model(results,maxLags)
  MSE            <- colMeans(error_prediction^2)

  return(list(table_data = table_data,
              MSE = MSE))

}

#' Title
#'
#' @param Theta0
#' @param Sigma0
#' @param n
#' @param Nsim
#' @param m
#' @param q0
#' @param maxLags
#'
#' @return
#' @export
#'
#' @examples
test_criteria_MA <- function(Theta0,Sigma0,n,Nsim,m,q0,maxLags) {

  results              <- array(NA,dim=c(Nsim,maxLags,2))
  error_prediction     <- matrix(NA,Nsim,maxLags)
  set.seed(12345)

  for(s in 1:Nsim) {

    print(paste('Simulation:',s,'of',Nsim))
    data <- sample_dataMA(Theta=Theta0,Sigma=Sigma0,q=q0,m=m,n=n+1)

    for(p in 1:maxLags) {

      model_AR  <- arima(as.data.frame(data[(maxLags-p+1):(n)]),order=c(p,0,0),
                         include.mean=FALSE)

      Sigma    <- (1/(n-maxLags-p))*t(model_AR$residuals)%*%model_AR$residuals

      results[s,p,1:2] <- unlist(criteria_func(Sigma=Sigma,n=n,p=p))
      # results[s,p,3] <- cross_validation_method(data,p)

      error_prediction[s,p] <- as.numeric(data[n+1]-
                                          predict(model_AR,n.ahead=1)$pred)

    }
  }

  table_data     <- select_model(results,maxLags)
  MSE            <- colMeans(error_prediction^2)

  return(list(table_data = table_data,
              MSE = MSE))

}


#' Title
#'
#' @param Theta0
#' @param Sigma0
#' @param n
#' @param Nsim
#' @param m
#' @param q0
#' @param maxLags
#'
#' @return
#' @export
#'
#' @examples
test_criteria_ARMA <- function(Phi0,Theta0,Sigma0,n,Nsim,m,p0,q0,maxLags) {

  results              <- array(NA,dim=c(Nsim,maxLags,3))
  error_prediction     <- matrix(NA,Nsim,maxLags)

  for(s in 1:Nsim) {

    print(paste('Simulation:',s,'of',Nsim))
    data <- sample_dataARMA(Phi=Phi0,Theta=Theta0,Sigma=Sigma0,
                            p0=p0,q0=q0,m=m,n=n+1)

    for(p in 1:maxLags) {

      model_AR  <- arima(as.data.frame(data[(maxLags-p+1):(n)]),order=c(p,0,0),
                         include.mean=FALSE)

      Sigma    <- (1/(n-maxLags-p))*t(model_AR$residuals)%*%model_AR$residuals

      results[s,p,1:2] <- unlist(criteria_func(Sigma=Sigma,n=n,p=p))
      results[s,p,3]   <- cv_AR(data[(maxLags-p+1):(n)],h=25,v=5,p=p)

      error_prediction[s,p] <- as.numeric(data[n+1]-
                                          predict(model_AR,n.ahead=1)$pred)

    }
  }

  table_data     <- select_model(results,maxLags)
  MSE            <- colMeans(error_prediction^2)

  return(list(table_data = table_data,
              MSE = MSE))

}
