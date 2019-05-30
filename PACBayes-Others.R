### Tolstikhin-Seldin bound
boundPBEB <- function(NMC, sigma2){
  c1 <- c2 <- 1.15
  # Computing the KL
  ratio <- initsigma2/(sigma2)
  KL <- d/2 * log(ratio) + d/2*(1/ratio-1) + (1/(2*initsigma2))*dot(ERMs[,3],ERMs[,3]) + log(sigma2GridSize)
    
  v1 <- ceiling((1/log(c1))*log(sqrt((exp(1)-2)*ntrain/(4*log(2/delta))))+1)
  v2 <- ceiling((1/log(c2))*log(.5*sqrt((ntrain-1)/(log(2/delta))+1)+.5))
  
  theta_samplesTS <- get_sample(type = distribution, mean=ERMs[,3], variance2=sigma2, NMC)
  Loss <- loss(Ytrain,predictor(Xtrain,theta_samplesTS))
  Ln <-  mean(Loss)
  VarTS <- mean(apply(X = Loss, MARGIN = 2, FUN = var))
  
  Vn <- VarTS + (1+c2)*sqrt(VarTS*(KL+log(2*v2/delta))/(2*(ntrain-1)))+2*c2*(KL+log(2*v2/delta))/(ntrain-1)
  Vnbar <- min(Vn,1/4)
  
  if(sqrt((KL+log(2*v1/delta))/((exp(1)-2)*Vnbar)) <= sqrt(ntrain)){
    val <- (1+c1)*sqrt((exp(1)-2)*Vnbar*(KL+log(2*v1/delta))/ntrain)
  }else{
    val <- 2*(KL+log(2*v1/delta))/ntrain
  }
  return(list(val=val,VarTS=VarTS,KL=KL))
}

### Maurer's bound
boundPBKL <- function(Ln,sigma2){
  # Computing the KL 
  ratio <- initsigma2/(sigma2)
  KL <- d/2 * log(ratio) + d/2*(1/ratio-1) + (1/(2*initsigma2))*dot(ERMs[,3],ERMs[,3]) + log(sigma2GridSize)
 
  fn <- function(p) Ln*log(Ln/p) + (1-Ln)*log((1-Ln)/(1-p)) - (KL + log(2*sqrt(ntrain)/delta))/ntrain; 
  val <- bisection(fn, Ln, 0.999)
  
  return(list(val=val,KL=KL))
}

boundCatoni <- function(Ln,sigma2){
  alpha <- rho
  # Computing the KL 
  ratio <- initsigma2/(sigma2)
  KL <- d/2 * log(ratio) + d/2*(1/ratio-1) + (1/(2*initsigma2))*dot(ERMs[,3],ERMs[,3]) + log(sigma2GridSize)
  
  D <- KL - log(delta)
  sqrtTerm <- sqrt(2*alpha * D/(ntrain*Ln*(1-Ln)))
  Num <-  1 - exp(-Ln*sqrtTerm - 
                    alpha/ntrain * (D + 2 * log(log(alpha^2*ntrain*sqrtTerm/log(alpha)))))
  Den <- 1 - exp(-sqrtTerm)
  val <- Num/Den
  return(list(val=val,KL=KL))
}


## Bounds on half the data
### Tolstikhin-Seldin bound
boundPBEB_half <- function(NMC,sigma2){
  c1 <- c2 <- 1.15
  
  nhalf <- ntrain/2
  # Computing the KL 
  KL <-  (1/(2*sigma2))*sum(square_diff(ERMs[,1],ERMs[,2]))+ log(sigma2GridSize)

  v1 <- ceiling((1/log(c1))*log(sqrt((exp(1)-2)*nhalf/(4*log(2/delta))))+1)
  v2 <- ceiling((1/log(c2))*log(.5*sqrt((nhalf-1)/(log(2/delta))+1)+.5))
  
  theta_samplesTS <- get_sample(type = distribution, mean=ERMs[,2],variance2=sigma2, NMC)
  Loss <- loss(Ytrain[1:nhalf],predictor(Xtrain[1:nhalf,],theta_samplesTS))
  Ln <-  mean(Loss)
  VarTS <- mean(apply(X = Loss, MARGIN = 2, FUN = var))
  
  Vn <- VarTS + (1+c2)*sqrt(VarTS*(KL+log(2*v2/delta))/(2*(nhalf-1)))+2*c2*(KL+log(2*v2/delta))/(ntrain-1)
  Vnbar <- min(Vn,1/4)
  
  if(sqrt((KL+log(2*v1/delta))/((exp(1)-2)*Vnbar)) <= sqrt(nhalf)){
    val <- (1+c1)*sqrt((exp(1)-2)*Vnbar*(KL+log(2*v1/delta))/nhalf)
  }else{
    val <- 2*(KL+log(2*v1/delta))/nhalf
  }
  return(list(val=val,VarTS=VarTS,KL=KL))
}

### Maurer's bound
boundPBKL_half <- function(NMC,sigma2){
  theta_samplesTS <- get_sample(type = distribution, variance2=sigma2, mean=ERMs[,2], NMC)
  nhalf <- ntrain/2
  Ln <- mean(loss(Ytrain[1:nhalf],predictor(Xtrain[1:nhalf,],theta_samplesTS)))
  # Computing the KL 
  KL <-  (1/(2*sigma2))*sum(square_diff(ERMs[,1],ERMs[,2])) + log(sigma2GridSize)
  
  fn <- function(p) Ln*log(Ln/p) + (1-Ln)*log((1-Ln)/(1-p)) - (KL + log(2*sqrt(nhalf)/delta))/nhalf; 
  val <- bisection(fn, Ln, 0.999)
  
  return(list(val=val,KL=KL))
}

boundCatoni_half <- function(NMC,sigma2){
  alpha <- 2
  theta_samplesTS <- get_sample(type = distribution, mean=ERMs[,2],variance2=sigma2, NMC)
  nhalf <- ntrain/2
  Ln <- mean(loss(Ytrain[1:nhalf],predictor(Xtrain[1:nhalf,],theta_samplesTS)))
  
  # Computing the KL 
  KL <-(1/(2*sigma2))*sum(square_diff(ERMs[,1],ERMs[,2]))+ log(sigma2GridSize)
  D <- KL - log(delta)
  
  sqrtTerm <- sqrt(2*alpha * D/(nhalf*Ln*(1-Ln)))
  Num <-  1 - exp(-Ln*sqrtTerm - alpha/nhalf * (D + 2 * log(log(alpha^2*nhalf*sqrtTerm)/log(alpha))))
  Den <- 1 - exp(-sqrtTerm)
  val <- Num/Den
  return(list(val=val,KL=KL))
}


