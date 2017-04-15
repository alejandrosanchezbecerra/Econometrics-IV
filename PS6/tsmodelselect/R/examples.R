#' Parameters required for Example 1
#'
#' @param Nsim Number of simulations
#'
#' @return
#' @export
#'
#' @examples
example1_PS6 <- function(Nsim,n) {

  config           <- list()
  config$Sigma0    <- matrix(1,1,1)

  config$Phi0      <- list()
  config$Phi0[[1]] <- matrix(0.7,1,1)

  config$Theta0      <- list()

  config$p0        <- 1
  config$q0        <- 0
  config$n         <- n
  config$Nsim      <- Nsim
  config$m         <- 1
  config$maxLags   <- 6

  return(config)

}

#' Parameters required for Example 2
#'
#' @param Nsim Number of simulations
#'
#' @return
#' @export
#'
#' @examples
example2_PS6 <- function(Nsim,n) {

  config            <- list()
  config$Sigma      <- matrix(1,1,1)

  config$Phi0      <- list()

  config$Theta0      <- list()
  config$Theta0[[1]] <- matrix(0.6,1,1)

  config$p0         <- 0
  config$q0         <- 1
  config$n          <- 100
  config$Nsim       <- Nsim
  config$m          <- 1
  config$maxLags    <- 6

  return(config)

}
