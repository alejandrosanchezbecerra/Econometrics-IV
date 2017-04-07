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
criteria_func <- function(Sigma,n,p){

  AIC     <- (log(det(Sigma)))+2*(p)/n
  SIC     <- (log(det(Sigma)))+p*log(n)/n

  return(list(AIC = AIC,
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

  select_model  <- matrix(NA,2,NumLags)
  max_criterion <- apply(results, c(3,1), min)
  min_array     <- aperm(sapply(rep(0,NumLags),
                                function(x) max_criterion,simplify='array'),
                         c(2,3,1))
  is_selected   <- results == min_array

  summary_selec <- apply(is_selected, c(3,2), sum)

}
