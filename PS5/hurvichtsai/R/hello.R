#' Title
#'
#' @param sample_size
#' @param intercept
#' @param slope
#' @param error_sd Standard deviation of the error term
#'
#' @return
#' @export
#'
#' @examples
dgp <- function(sample_size, intercept, slope, error_sd){
  x <- runif(sample_size)
  epsilon <- rnorm(sample_size, sd = error_sd)
  y <- intercept + slope * x + epsilon
  return(y)
}
