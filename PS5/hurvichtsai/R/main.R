
#' Replicate the results in Hurvich and Tsai
#'
#' @param Phi0   A list with p transition matrices of a VAR(p)
#' @param Sigma0 Covariance matrix of error terms
#' @param n      Number of observations
#' @param Nsim   Number of simulation draws
#' @param m      Dimensions of multivariate time series
#' @param NumLags Maximum number of lag-order for evaluation
#'
#' @return A list of two tables summarizing the results. The first computes the
#' frequency with which each model order is selected. The second computes the
#' average
#' @export
#'
#' @examples
#' do.call(replicate_results,example1(Nsim=100))
#' do.call(replicate_results,example2(Nsim=100))
#'
replicate_results <- function(Phi0,Sigma0,n,Nsim,m,p0,NumLags) {

  results <- array(NA,dim=c(Nsim,NumLags,4))
  delta   <- matrix(NA,Nsim,NumLags)
  set.seed(12345)

  for(s in 1:Nsim) {

    print(paste('Simulation:',s,'of',Nsim))
    data <- sample_data(Phi=Phi0,Sigma=Sigma0,p=p0,m=m,n=n,NumLags=NumLags)

    for(p in 1:NumLags) {

      param   <- (mle_VAR(data,m,p,n,NumLags))
      beta    <- param$beta
      Sigma   <- param$Sigma
      Phi     <- Beta_to_Phi(beta,m,p)

      delta[s,p] <- estimate_Delta(Phi,Sigma,Phi0,Sigma0,m,p,p0,n)
      results[s,p,] <- unlist(criteria_func(Sigma=Sigma,m=m,n=n,p=p))
    }
  }

  table_data     <- select_model(results,NumLags)

  delta_expected <- colMeans(delta)
  table_figure   <- rbind(apply(results, c(3,2), mean),
                          delta_expected)

  return(list(table_data = table_data,
              figure_data = table_figure))

}
