#' Title
#'
#' @return
#' @export
#'
#' @examples
example1 <- function() {

  config           <- list()
  config$Sigma0    <- matrix(c(1,0,0,1),2,2)
  config$Phi0      <- list()
  config$Phi0[[1]] <- matrix(c(-1,-1.5,0.96,1.4),2,2)

  config$n         <- 40
  config$Nsim      <- 100
  config$m         <- 2
  config$NumLags   <- 6

  return(config)

}

#' Title
#'
#' @return
#' @export
#'
#' @examples
example2 <- function() {

  config           <- list()
  config$Sigma0    <- matrix(c(1,-0.08,-0.08,1),2,2)

  config$Phi0      <- list()
  config$Phi0[[1]] <- matrix(c(0.5,0.2,-0.3,0.65),2,2)
  config$Phi0[[2]] <- matrix(c(-0.5,0,0.3,-0.4),2,2)

  config$n         <- 40
  config$Nsim      <- 100
  config$m         <- 2
  config$NumLags   <- 6

  return(config)

}