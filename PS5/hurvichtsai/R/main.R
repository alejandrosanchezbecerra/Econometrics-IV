
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
replicate_results <- function(Phi0,Sigma0,n,Nsim,m,NumLags) {

  results <- array(NA,dim=c(Nsim,NumLags,4))
  set.seed(12345)

  for(s in 1:Nsim) {

    print(paste('Simulation:',s,'of',Nsim))
    data <- sample_data(Phi=Phi0,Sigma=Sigma0,p=1,m=m,n=n,NumLags=NumLags)

    for(p in 1:NumLags) {
      Sigma         <- mle_VAR(data,m,p,n,NumLags)$Sigma
      results[s,p,] <- unlist(criteria_func(Sigma=Sigma,m=m,n=n,p=p))
    }
  }

  tables <- select_model(results,NumLags)
  return(tables)

}
