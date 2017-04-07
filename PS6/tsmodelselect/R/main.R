test_criteria_AR <- function(Phi0,Sigma0,n,Nsim,m,p0,maxLags) {

  results <- array(NA,dim=c(Nsim,maxLags,2))
  delta   <- matrix(NA,Nsim,maxLags)
  set.seed(12345)

  for(s in 1:Nsim) {

    print(paste('Simulation:',s,'of',Nsim))
    data <- sample_dataAR(Phi=Phi0,Sigma=Sigma0,p=p0,m=m,n=n)

    for(p in 1:maxLags) {

      model_AR  <- arima(as.data.frame(data[(maxLags-p+1):(n)]),order=c(p,0,0),
                         include.mean=FALSE)

      Sigma    <- (1/(n-maxLags-p))*t(model_AR$residuals)%*%model_AR$residuals

      results[s,p,1:2] <- unlist(criteria_func(Sigma=Sigma,n=n,p=p))
      # results[s,p,3] <- cross_validation_method(data,p)

    }
  }

  table_data     <- select_model(results,maxLags)
  table_figure   <- rbind(apply(results, c(3,2), mean))

  return(list(table_data = table_data,
              figure_data = table_figure))

}


test_criteria_MA <- function(Theta0,Sigma0,n,Nsim,m,q0,maxLags) {

  results <- array(NA,dim=c(Nsim,maxLags,2))
  delta   <- matrix(NA,Nsim,maxLags)
  set.seed(12345)

  for(s in 1:Nsim) {

    print(paste('Simulation:',s,'of',Nsim))
    data <- sample_dataMA(Theta=Theta0,Sigma=Sigma0,q=q0,m=m,n=n)

    for(p in 1:maxLags) {

      model_AR  <- arima(as.data.frame(data[(maxLags-p+1):(n)]),order=c(p,0,0),
                         include.mean=FALSE)

      Sigma    <- (1/(n-maxLags-p))*t(model_AR$residuals)%*%model_AR$residuals

      results[s,p,1:2] <- unlist(criteria_func(Sigma=Sigma,n=n,p=p))
      # results[s,p,3] <- cross_validation_method(data,p)

    }
  }

  table_data     <- select_model(results,maxLags)
  table_figure   <- rbind(apply(results, c(3,2), mean))

  return(list(table_data = table_data,
              figure_data = table_figure))

}
