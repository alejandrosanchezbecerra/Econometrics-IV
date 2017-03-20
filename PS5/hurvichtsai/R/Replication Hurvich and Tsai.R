#********************************************#
#************* ENVIRONMENT ******************#
#********************************************#

# install.packages('vars')
# install.packages('psych')
# install.packages('numDeriv')
library('MASS')
library('vars')
library('psych')
library('numDeriv')

#********************************************#
#*************** AR(1) MODEL ****************#
#********************************************#

Sigma0  <- matrix(c(1,0,0,1),2,2)
Phi1    <- matrix(c(-1,-1.5,0.96,1.4),2,2)

Phi0      <- list()
Phi0[[1]] <-  Phi1

n       <- 40
Nsim    <- 100
m       <- 2
NumLags <- 6
results <- array(NA,dim=c(Nsim,NumLags,4))

set.seed(12345)

for(s in 1:Nsim) {
  
  print(s)
  
  data <- sample_data(Phi=Phi0,Sigma=Sigma0,p=1,m=m,n=n,NumLags=NumLags)
  
  for(p in 1:NumLags) {
    
    Sigma         <- mle_VAR(data,m,p,NumLags)$Sigma
    results[s,p,] <- unlist(criteria_func(Sigma,m,n,p))
    
  }
  
}

select_model  <- matrix(NA,4,NumLags)
max_criterion <- apply(results, c(3,1), min)
min_array     <- aperm(sapply(rep(0,NumLags),
                              function(x) max_criterion,simplify='array'),
                       c(2,3,1))

is_selected   <- results == min_array
summary_selec <- apply(is_selected, c(3,2), sum)
summary_selec

#********************************************#
#*************** AR(2) MODEL ****************#
#********************************************#

Sigma0  <- matrix(c(1,-0.08,-0.08,1),2,2)
Phi1    <- matrix(c(0.5,0.2,-0.3,0.65),2,2)
Phi2    <- matrix(c(-0.5,0,0.3,-0.4),2,2)

Phi      <- list()
Phi[[1]] <- Phi1
Phi[[2]] <- Phi2

n       <- 40
Nsim    <- 50
m       <- 2
NumLags <- 6
results <- array(NA,dim=c(Nsim,NumLags,4))

set.seed(12345)

for(s in 1:Nsim) {
  
  print(s)

  data <- sample_data(Phi=Phi0,Sigma=Sigma0,p=2,m=m,n=n,NumLags=NumLags)
  
  for(p in 1:NumLags) {
    
    Sigma         <- mle_VAR(data,m,p,NumLags)$Sigma
    results[s,p,] <- unlist(criteria_func(Sigma,m,n,p))
    
  }
  
}

select_model  <- matrix(NA,4,NumLags)
max_criterion <- apply(results, c(3,1), min)
min_array     <- aperm(sapply(rep(0,NumLags),
                              function(x) max_criterion,simplify='array'),
                       c(2,3,1))

is_selected   <- results == min_array
summary_selec <- apply(is_selected, c(3,2), sum)

summary_selec
