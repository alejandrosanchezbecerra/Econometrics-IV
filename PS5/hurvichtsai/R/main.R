
#' Replicate the results in Hurvich and Tsai
#'
#' @param Phi0   A list with p transition matrices of a VAR(p)
#' @param Sigma0 Covariance matrix of error terms
#' @param n      Number of observations
#' @param Nsim   Number of simulation draws
#' @param m      Dimensions of multivariate time series
#' @param NumLags Maximum number of lag-order for evaluation
#'
#' @return A list with tables summarizing the results.
#' @export
#'
#' @examples
replicate_results <- function(Phi0,Sigma0,n,Nsim,m,p0,NumLags) {

  results <- array(NA,dim=c(Nsim,NumLags,4))
  delta   <- matrix(NA,Nsim,NumLags)
  set.seed(12345)

  for(s in 1:Nsim) {

    print(paste('Simulation:',s,'of',Nsim))
    data <- sample_data(Phi=Phi0,Sigma=Sigma0,p=p0,m=m,n=n,NumLags=NumLags)

    for(p in 1:NumLags) {

      # var4_model  <- vars::VAR(as.data.frame(data[(NumLags-p+2):(n+NumLags+1),]),
      #                          p = p,type="none")
      # summary_VAR <- summary(var4_model)
      # Sigma       <- summary_VAR$covres
      # Sigma         <- (mle_VAR(data,m,p,n,NumLags))$Sigma
      param   <- (mle_VAR(data,m,p,n,NumLags))
      beta    <- param$beta
      Sigma   <- param$Sigma
      Phi     <- Beta_to_Phi(beta,m,p)

      delta[s,p] <- estimate_Delta(Phi,Sigma,Phi0,Sigma0,m,p,p0,n)

      # model_VAR  <- vars::VAR(as.data.frame(data[(NumLags-p+1):(n+NumLags),]),
      #                        p = p,type="none")
      # Bmat <- matrix(0,2,2)
      # Bmat[lower.tri(Bmat,diag=TRUE)] <- NA
      # Sigma <- SVAR(model_VAR,estmethod="direct",Amat=diag(2),Bmat=Bmat)$Sigma.U/100

      # print(Sigma)

      results[s,p,] <- unlist(criteria_func(Sigma=Sigma,m=m,n=n,p=p))
    }
  }

  table_data     <- select_model(results,NumLags)

  delta_expected <- colMeans(delta)
  table_figure   <- rbind(apply(results, c(3,2), mean),
                          delta_expected)

  return(list(table_data = table_data,
              figue_data = table_figure))

}
