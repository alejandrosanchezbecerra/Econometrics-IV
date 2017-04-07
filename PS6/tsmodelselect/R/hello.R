# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

sample_dataAR <- function(Phi,Sigma,p,m,n,NumLags,...) {

  epsilon <- MASS::mvrnorm(n+NumLags+1,mu=rep(0,m),Sigma=Sigma)
  data    <- matrix(0,n+NumLags+1,2)

  for(t in (NumLags+1):(n+NumLags+1)) {

    data[t,] <- epsilon[t,]
    order    <- 1

    while(order <= p) {
      data[t,] <- data[t,]+matrix(c(data[t-order,]),1,2)%*%t(Phi[[order]])
      order    <- order+1
    }

  }

  return(data)

}
