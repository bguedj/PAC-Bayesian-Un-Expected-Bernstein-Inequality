## Helper functions
square_diff <- function(a,b){
  return((a-b)^2)
}

bisection <- function(fn, a, b, tol = 1e-8) {
  # If the signs of the function at the evaluated points, a and b, stop the function and return message.
  if (fn(a) * fn(b) > 0)
    stop('signs of f(a) and f(b) are the same')
  p <-  (a + b)/2;
  err = abs(fn(p));
  while (err > tol){
    ifelse(fn(a)*fn(p)<0, b <- p, a <- p)
    p <- (a + b)/2; 
    err <- abs(fn(p));
  }
  return(p)
}

get_sample <- function(type, mean, variance2, n_samples=1){
  return(t(mvrnorm(n = n_samples, mu = mean, Sigma = variance2*diag(d))))
}
## End of helper functions

buildERMsequenceFast <- function(eps = 1e-16){
  ERMsequence <- matrix(nrow = d, ncol = 3, data = 0)
  for(ii in c(1,2,3))   ## change here for the fast version
  {
    if(ii==1) {setPoints <- 1:(ntrain/2)}        ## ERM on first half
    if(ii==2) {setPoints <- (ntrain/2+1):ntrain} ## ERM on second half
    if(ii==3) {setPoints <- 1:ntrain}            ## ERM on full sample
    
    ## Definining the objective and the corresponding gradient
    fn <- function(theta){
      X <- Xtrain[setPoints,]
      Y <- Ytrain[setPoints]
      Z <-  as.numeric(X%*%theta)
      return(-mean(Y*log(sigmoid(Z))+(1-Y)*log(1-sigmoid(Z))) +lambda*dot(theta,theta)/2)
    }
    gr <- function(theta){
      X <- Xtrain[setPoints,]
      Y <- Ytrain[setPoints]
      Z <- as.numeric(X%*%theta)
      return(-as.numeric((Y-sigmoid(Z))%*%X/length(setPoints)) + lambda*theta)
    }
    par <- numeric(d)
    ERMsequence[,ii] <- optim(par, fn, gr = gr, method = "BFGS", hessian = TRUE)$par
  }
  return(ERMsequence)
}

computeExpectation <- function(ERMfull,NMC, sigma2){
  theta_samples <- get_sample(type = distribution, mean=ERMfull, variance2=sigma2, n_samples=NMC)
  result <- loss(Ytrain,predictor(Xtrain,theta_samples))
  return(mean(result))
}

COMP <- function(ERMfull,ERM1,ERM2,sigma2){
  val1 <- sum(square_diff(ERMfull,ERM1))
  val2 <- sum(square_diff(ERMfull,ERM2))
   
  result <-   (val1+val2)/(2*sigma2)  
  return(result)
}

OptimEtaVn <- function(NMC, vnTerm, compTerm){
  result <- numeric(etaGridSize)
  for(etaInd in 1:etaGridSize){
    eta <- etaGrid[etaInd]
    result[etaInd] <- vnTerm*(-log(1-eta*b)/(eta*b*b) - 1/b) +
      (compTerm  + 2*log(sigma2GridSize*etaGridSize/delta))/(eta*ntrain)
  }
  argmin <- which.min(result)
  val <- result[argmin]
  return(list(val=val,etaOpt=argmin))
}

OptimEtaVnPrime <- function(NMC, vnTermPrim){
  result <- numeric(etaGridSize)
  for(etaInd in 1:etaGridSize)
  {
    eta <- etaGrid[etaInd]
    result[etaInd] <- vnTermPrim*(-log(1-eta*b)/(eta*b*b) - 1/b) + 
      log(etaGridSize*sigma2GridSize/delta)/(eta*ntrain)
  }
  argmin <- which.min(result)
  val <- result[argmin]
  return(list(val=val,etaOpt=argmin))
}

VnTerm <- function(ERMfull,ERM1,ERM2,NMC,sigma2){
  theta_samples <- get_sample(type = distribution, mean=ERMfull, variance2=sigma2, n_samples=NMC)
  result  <- matrix(nrow = ntrain, ncol = NMC, data = NA)
  loss1 <- t(matrix(loss(Ytrain[1:(ntrain/2)],predictor(Xtrain[1:(ntrain/2),],ERM1)), nrow=NMC, ncol=ntrain/2, byrow=TRUE))
  loss2 <- t(matrix(loss(Ytrain[(ntrain/2+1):ntrain],predictor(Xtrain[(ntrain/2+1):ntrain,],ERM2)), nrow=NMC, ncol=ntrain/2, byrow=TRUE))
  
  result[1:(ntrain/2),] <- square_diff(loss(Ytrain[1:(ntrain/2)],predictor(Xtrain[1:(ntrain/2),],theta_samples)),loss1)
  result[(ntrain/2+1):ntrain,]<-square_diff(loss(Ytrain[(ntrain/2+1):ntrain],predictor(Xtrain[(ntrain/2+1):ntrain,], theta_samples)),loss2)
  return(mean(result))
}

VnPrimeTerm <- function(ERM1,ERM2,NMC){
  res <- numeric(ntrain)
  res[1:(ntrain/2)] <- loss(Ytrain[1:(ntrain/2)],predictor(Xtrain[1:(ntrain/2),],ERM1))^2
  res[(ntrain/2+1):ntrain] <- loss(Ytrain[(ntrain/2+1):ntrain],predictor(Xtrain[(ntrain/2+1):ntrain,],ERM2))^2
  return(mean(res))
}

## The following bound does not compute the empirical error
mainBoundProba <- function(NMC,sigma2){
  vnTerm <- VnTerm(ERMfull = ERMs[,3],
                          ERM1 = ERMs[,2],
                          ERM2 = ERMs[,1],
                          NMC = NMC, sigma2=sigma2)
  compTerm <- COMP(ERMfull = ERMs[,3],
                   ERM1 = ERMs[,2],
                   ERM2 = ERMs[,1],
                   sigma2 = sigma2)
  tmp1 <- OptimEtaVn(NMC = NMC,
                     vnTerm = vnTerm,
                             compTerm = compTerm)
  etaOpt1 <- tmp1$etaOpt
  val1 <- tmp1$val
  
  vnTermPrim <- VnPrimeTerm(ERM1 = ERMs[,2],
                                ERM2 = ERMs[,1],
                                NMC = NMC)
  tmp2 <- OptimEtaVnPrime(NMC = NMC,
                          vnTermPrim = vnTermPrim)
  etaOpt2 <- tmp2$etaOpt
  val2 <- tmp2$val
  
  val <- val1 + val2
  return(list(val=val,val1=val1,val2=val2,vnTerm=vnTerm,vnTermPrim=vnTermPrim,compTerm=compTerm,etaOpt=c(etaOpt1,etaOpt2)))
}






