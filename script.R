rm(list=ls())
set.seed(123)
library(MASS)
library(tictoc)
library(Matrix)
library(pracma)
library(stringr)
library(caret)
library(R.matlab)

par(mfrow=c(1,1))
path <- "./"

## Experimental setup
data_option = "adults"          # Options are: "sigmoid-synthetic", "haberman", "breast-cancer", 
                                  # "tictactoe", "bank-notes", "kr-vs-kp", "spam", "mushroom", "adults"
problem_type = "classification"   # No other problem type is supported currently  
distribution <- "gaussian"        # No other distribution is supported currently 
lambda <- 0.01                    # Regularization for the logistic regression
                      
initsigma2 <- 0.5                 # Prior variance for the Gaussian-ERM distribution 
half <-  F                        # Set to true so that half the data is used to build a prior for PACBayes-Others (i.e. not our method)

## following block is only relevant for synthetic data
nb.seq <- 10
nmin <- 800
nmax <- 8000
nb.grid <- 2*round(seq(nmin,nmax,length.out = nb.seq)/2)
#theta_star <-  c(3,1,4,1,5,9,2,6,5,3,5,8,9,7,9,3,2,3,8,4,6,2,6,4,3,3,8,3,2,7,9,5,0,2,8,8,4,1,9,7,1,6,9,3,9,9,3,7,5,1)
theta_star <- c(3,1,4,1,5,9,2,6,5,3)
d <- length(theta_star)
noise <- T
proportion <- 0.1

## Global constants
rho <- 2            # This is for the grid-eta discretization
delta <- 0.05       # The probability Threshold for the bound
NMC <- 100          # how many Monte Carlo iterations for the Gaussian expectation?

## Definining the loss and the predictor type
loss <- function(a,b) abs(a-b)
b <- 1    # upper bound on the loss
predictor <- function(x,theta) round(1/(1+exp(-x%*%theta)))

## End of Experimental Setup
source(paste(path, "gen-data.R", sep="/"))
source(paste(path, "PACBayes2.R", sep="/"))
source(paste(path, "PACBayes-Others.R", sep="/"))
str <- paste(c(data_option,"imformedPrior",half),collapse='-')

# Initializing 
ifelse(grepl("synthetic",data_option, fixed=TRUE), nbRepet <- 10, nbRepet <- 5)
results <- array(dim = c(nbRepet,4,nb.seq), data = Inf)
bestSigma2 <- array(dim = c(nbRepet,nb.seq), data = Inf)
LnERMtrain <- LnERMtest <- array(dim = c(nbRepet, nb.seq), data = NA)
Vn <- VnPrim  <- VarTS  <- array(dim = c(nbRepet, nb.seq), data = NA)
comp <- KL <- array(dim = c(nbRepet, nb.seq), data = NA)
val1 <- val2 <- array(dim = c(nbRepet, nb.seq), data = NA)

pb <- txtProgressBar(min = 0, max = nb.seq, style = 3)
for(inb in 1:nb.seq){
  for(irepet in 1:nbRepet){
    ## Generate data
    if(irepet==1 || grepl("synthetic", data_option, fixed=TRUE)){
      # Only regenerate for synthetic data
      tmp <- gendata(option = data_option)
      Xfull <- tmp$X
      Yfull <- tmp$Y
      d <- tmp$d
      nfull <- length(Yfull)
      if(irepet==1)
        print(c("Size of the full data=",nfull))
    }
    # 5-fold train-test split 
    ifelse(irepet %% 5==0,ind <- 4,ind <- irepet %% 5 - 1)
    ifelse((nfull-floor(nfull/5))%% 2==0, ntest<-floor(nfull/5), ntest <-ceil(nfull/5))
    testIND <- c((ind*ntest+1):min((ind+1)*ntest,nfull))
    Xtrain <- Xfull[-testIND,]
    Ytrain <- Yfull[-testIND]
    Xtest <- Xfull[testIND,]
    Ytest <- Yfull[testIND]
  
    ntrain <- length(Ytrain)
    if(irepet==1)
      print(c("Size of the training data=",ntrain))
    
    etaGridSize <- ceil(log(0.5*sqrt(ntrain/log(1/delta)))/log(rho))
    etaGrid <- numeric(etaGridSize)
    for(jj in 1:(etaGridSize))
      etaGrid[jj] <- 1/(b*rho^(jj)) 
    
    # Build grid of posterior variance for the Gaussian-ERM distribution 
    sigma2GridSize <- ceil(log(ntrain)/log(2))
    sigma2Grid <- numeric(sigma2GridSize)
    for(jj in 1:sigma2GridSize)
      sigma2Grid[jj] <- 1/(2^(jj))
    
    ## Compute ERMs
    ERMs <- buildERMsequenceFast()

    ## Loop over the grid of sigma2 and pick the best one for each method
    for(sigma2 in sigma2Grid){    
      ## Compute the shared empirical error
      theta_samplesTS <- get_sample(type = distribution, mean=ERMs[,3],variance2=sigma2, NMC)
      Ln <- mean(loss(Ytrain,predictor(Xtrain,theta_samplesTS)))
      
      ## Our bound 
      tmpBProb <- mainBoundProba(NMC,sigma2)
      if(Ln + tmpBProb$val < results[irepet,1,inb]){
        results[irepet,1,inb] <- Ln + tmpBProb$val
        Vn[irepet,inb] <- tmpBProb$vnTerm
        VnPrim[irepet,inb] <- tmpBProb$vnTermPrim
        comp[irepet,inb] <- tmpBProb$compTerm
        val1[irepet,inb] <- tmpBProb$val1
        val2[irepet,inb] <- tmpBProb$val2
        bestSigma2[irepet,inb] <- sigma2
      }
      
      ## TS bound
      ifelse(half,tmpBEB<- boundPBEB_half(NMC,sigma2),tmpBEB<-boundPBEB(NMC,sigma2))
      if(Ln + tmpBEB$val < results[irepet,2,inb]){
        results[irepet,2,inb] <- Ln + tmpBEB$val
        VarTS[irepet,inb] <- tmpBEB$VarTS
      }
      ## Maurer bound
      ifelse(half,tmpBKL<- boundPBKL_half(NMC, sigma2),tmpBKL<-boundPBKL(Ln, sigma2))
      if(tmpBKL$val < results[irepet,3,inb]){
        results[irepet,3,inb] <-  tmpBKL$val
        KL[irepet,inb] <-  tmpBKL$KL
      }
      
      ## Catoni bound
      ifelse(half, tmpBCT <-boundCatoni_half(NMC,sigma2),tmpBCT <-boundCatoni(Ln,sigma2))
      if(tmpBCT$val < results[irepet,4,inb]){
        results[irepet,4,inb] <- tmpBCT$val
      }
      
      LnERMtrain[irepet,inb] <- mean(loss(Ytrain,predictor(Xtrain,ERMs[,3])))
      LnERMtest[irepet,inb] <- mean(loss(Ytest,predictor(Xtest,ERMs[,3])))
    }
  }
  setTxtProgressBar(pb, inb)
  if(!grepl("synthetic",data_option, fixed=TRUE)){
    results <- results[,,1]
    bestSigma2 <- bestSigma2[,1]
    LnERMtrain <- LnERMtrain[,1]
    LnERMtest <- LnERMtest[,1]
    Vn <- Vn[,1]
    VnPrim <- VnPrim[,1]
    VarTS <- VarTS[,1]
    comp <- comp[,1]
    KL <- KL[,1]
    val1 <- val1[,1]
    val2 <- val2[,1]
    break
  }
}
str <- paste(c(str,d),collapse="-")

if(!grepl("synthetic",data_option, fixed=TRUE)){
  means <- apply(results, 2, mean)
  print(paste(c(data_option, ". Our Bound=",  round(means[1],3),
              ", Maurer bound=", round(means[3],3),
             ", TS bound=", round(means[2],3), 
             ", Catoni bound=", round(means[4],3)),collapse = ""))
  print(paste(c("ERM test error=", mean(LnERMtest)),collapse=""))
  
  print(paste(c(round(mean(LnERMtest),3), " & ",  round(means[1],3),
                " & ", round(means[3],3),
                " & ",round(means[2],3), 
                " & ", round(means[4],3)), collapse = ""))
}else{
  MeanBProb <- apply(X = results[,1,], MARGIN = 2, FUN = mean)
  MeanBTS <- apply(X = results[,2,], MARGIN = 2, FUN = mean)
  MeanBKL <- apply(X = results[,3,], MARGIN = 2, FUN = mean)
  MeanBCatoni <- apply(X = results[,4,], MARGIN = 2, FUN = mean)
  
  ## Saving results 
  writeMat(paste(c(path, "/save/results_",str,".mat"), collapse = ""), labpcexport = results)
  writeMat(paste(c(path, "/save/Vn_",str,".mat"), collapse = ""), labpcexport = Vn)
  writeMat(paste(c(path, "/save/VnPrim_",str,".mat"), collapse = ""), labpcexport = VnPrim)
  writeMat(paste(c(path, "/save/LnERMtrain_",str,".mat"), collapse = ""), labpcexport = LnERMtrain)
}




  
