#' Compute the four model-selection criteria indexes
#'
#' @param Sigma MLE estimated covariance matrix of error terms
#' @param m     Dimension of multivariate time series
#' @param n     Number of observations
#' @param p     Order of estimated VAR(p)
#'
#' @return A list with (AIC,AIC_c,AIC_cBD,SIC).
#' @export
#'
#' @examples
criteria_func <- function(Sigma,m,n,p){

  b <- n/(n-(p*m+m+1))

  AIC     <- n*(log(det(Sigma))+m)+2*(p*(m^2)+(m/2)*(m+1))
  AIC_c   <- n*(log(det(Sigma))+m)+2*b*(p*(m^2)+(m/2)*(m+1))
  AIC_cBD <- n*(log(det(Sigma))+m)+2*(p*(m^2)+1)*n/(n-p*(m^2)-2)
  SIC     <- n*(log(det(Sigma))+m)+p*(m^2)*log(n)

  return(list(AIC = AIC,
              AIC_c = AIC_c,
              AIC_cBD = AIC_cBD,
              SIC=SIC))

}

#' Compute frequency of selected criteria
#'
#' @param results (Nsim x NumLags x 4) array with computed criteria
#' @param NumLags Maximum number of lags considered for evaluation
#'
#' @return 4 x NumLags matrix with AIC,AICc, AICcBD and SIC criteria,
#' on each row, respectively.
#' @export
#'
#' @examples
select_model <- function(results,NumLags) {

  select_model  <- matrix(NA,4,NumLags)
  max_criterion <- apply(results, c(3,1), min)
  min_array     <- aperm(sapply(rep(0,NumLags),
                                function(x) max_criterion,simplify='array'),
                         c(2,3,1))
  is_selected   <- results == min_array

  summary_selec <- apply(is_selected, c(3,2), sum)

}

#' Compute E[Delta(hat{beta},hat{Sigma})].
#'
#' @param Phi   Estimated Phi, p-list with mxm transition matrices.
#' @param Sigma Estimated Sigma, mxm matrix with covariance matrix
#' @param Phi0  True Phi, p-list with mxm transition matrices.
#' @param Sigma0 True Sigma, mxm matrix with covariance matrix
#' @param m Number of dimensions of multivariate time series.
#' @param p Order of estimated VAR(p)
#' @param p0 Order of true VAR(p0)
#' @param n Number of observations
#'
#' @return
#' @export
#'
#' @examples
estimate_Delta <- function(Phi,Sigma,Phi0,Sigma0,m,p,p0,n) {

  H  <- construct_H(Phi,m,p,n)
  H0 <- construct_H(Phi0,m,p0,n)

  delta <- n*(log(det(Sigma)))+
           psych::tr((t(H)%*%kronecker(diag(n),solve(Sigma))%*%H)%*%
                     solve(H0)%*%kronecker(diag(n),solve(Sigma0))%*%t(solve(H0)))

  return(delta)

}
