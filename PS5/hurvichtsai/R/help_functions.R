
criteria_func <- function(Sigma,m,n,p){
  
  b <- n/(n-(p*m+m+1))
  
  AIC     <- n*(log(det(Sigma))+m)+2*(p*(m^2)+(m/2)*(m+1))
  AIC_c   <- n*(log(det(Sigma))+m)+2*b*(p*(m^2)+(m/2)*(m+1))
  AIC_cBD <- n*(log(det(Sigma))+m)+2*(p*(m^2)+1)*n/(n-p*(m^2)-2)
  SIC     <- n*(log(det(Sigma))+m)+p*(m^2)*log(n)
  
  return(list(AIC = AIC,
              SIC=SIC,
              AIC_c = AIC_c,
              AIC_cBD = AIC_cBD
  ))
  
}

sample_data <- function(Phi,Sigma,p,m,n,NumLags) {
  
  epsilon <- mvrnorm(n+NumLags+1,mu=rep(0,m),Sigma=Sigma)
  data    <- matrix(0,n+NumLags+1,2)
  
  for(t in (NumLags+1):(n+NumLags+1)) {
    
    data[t,] <- epsilon[t,]
    
    order    <- 1
    
    while(order <= p) {
      
      data[t,] <- data[t,]+matrix(c(data[t-order,]),1,2)%*%Phi[[order]]
      order    <- order+1
      
    }
    
  }
  
  return(data)
  
}

lik_VAR <- function(param,Y,X,p,m) {

  # param <- init_guess
  # print(p)
  # print(m)
      
  coeff <- get_parameters(param,p,m)
  beta  <- coeff$beta
  Sigma <- coeff$Sigma
  
  lik   <- n*log(det(Sigma))+
           tr((Y-X%*%beta)%*%solve(Sigma)%*%t(Y-X%*%beta))
  
  return(lik)
  
}

mle_VAR <- function(data,m,p,NumLags) {
  
  Y <- data[(NumLags+1):(n+NumLags),]
  X <- get_regressorMat(data,NumLags+1,n,p)
  
  init_guess <- get_param_vec(beta = matrix(0,m^2,p),
                              Sigma = diag(m) )
  
  model     <- optim(par=init_guess, fn=function(par) {lik_VAR(par,Y,X,p,m)})
  est_coef  <- model$par
  param     <- get_parameters(est_coef,p,m)
  
  return(param)
  
}

get_regressorMat <- function(Z,n0,n,p) {
  
  if(p == 1) {
    
    X <- Z[(n0-1):((n0-1)+(n-1)),]
    
  } else if (p >= 1) {
    
    X <- Z[(n0-1):((n0-1)+(n-1)),]
    
    for( j in 2:p) {
      
      X <- cbind(X,Z[(n0-j):((n0-j)+(n-1)),])
      
    }
    
  }
  
  return(X)
  
}

get_parameters <- function(param,p,m) {
  
  beta      <- matrix(param[1:(p*(m^2))],m*p,m)
  cholSigma <- matrix(0,m,m)
  cholSigma[lower.tri(cholSigma,diag=TRUE)] <-  param[(p*(m^2)+1):length(param)]
  Sigma     <- cholSigma%*%t(cholSigma )
   
  return(list(beta = beta,Sigma=Sigma))

}

get_param_vec <- function(beta,Sigma) {
  
  Sigma <- Sigma0
  cholSigma <- t(chol(Sigma))
  cholSigma <- cholSigma[lower.tri(cholSigma,diag=TRUE)]
  
  return(c(as.vector(beta),cholSigma))
  
}
