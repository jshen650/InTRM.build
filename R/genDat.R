#' Generate a simulated data set containing missingness
#'
#' @param size Size of simulated data set
#' @param propMiss Proportion of missingness (30% or 50%)
#' @param seedNum Seed used for random generation
#' @param settingMiss Missingness setting, which could be Original or Modified (to be added later)
#' @import mvtnorm MASS
#'
#' @return List of data frames containing the full data set and the data set with missingness

genDat <- function(size, propMiss, seedNum, settingMiss="Original", settingNumber, splitData=FALSE){
  ## code that will automatically generate a single data set given argument specifications
  ## Original Setting = manuscript settings
  # Additional settings may be added later and have to do with explorations discussed in Supplement

  ## consider 5 covariates
  n <- size
  p <- 5

  set.seed(seedNum)


  X_2.colMean <- c(-1, -0.6, -1.6, 0.5, 0.9)
  X_2.stddev <- c(2, 1.7, 1.6, 2.1, 2.2)

  # correlation matrix
  X_2.corMat <-  matrix(0.1, nrow = p, ncol = p) ## exchangeable correlation of 0.1
  diag(X_2.corMat) <- 1

  # covariance matrix
  X_2.covMat <- X_2.stddev %*% t(X_2.stddev) * X_2.corMat

  X_2.mean <-  mvrnorm(n = n, mu = X_2.colMean, Sigma = X_2.covMat, empirical = FALSE)

  A_2 <- sign(runif(n=n, min=-1, max=1))

  X_2 = cbind (1, X_2.mean)
  sm_gamma2 <- c(6.8, 1.5, 2.0, 1.5, 1.6, 2.5)
  phi2 <- c(0.5, -0.2,0.6, 0.1, 0.5, 0.6 )

  mu_Y <- X_2%*%sm_gamma2 + A_2*(X_2%*%phi2)
  eps <- rnorm(n=n, mean=0, sd=5)
  Y <- mu_Y + eps

  # missingness mechanism
  expit <- function(val){
    exp(val)/(1+exp(val))
  }

  if(settingNumber==1){
    if(propMiss==30){
      if(size==300  & settingMiss=="Original"){
        cm1 <- cbind(1, X_2.mean[,1], A_2*X_2.mean[,2])
        #psi <- c(-1.6, 2.3, -3.2)
      }

      if(size==600  & settingMiss=="Original"){
        cm1 <- cbind(1, X_2.mean[,1], A_2*X_2.mean[,2])
        #psi <- c(-1.6, 2.3, -3.2)
      }

      if(size==1200  & settingMiss=="Original"){
        cm1 <- cbind(1, X_2.mean[,1], A_2*X_2.mean[,2])
        psi <- c(-1.6, 2.3, -3.2)
      }
    }

    if(propMiss==50){
        if(size==300  & settingMiss=="Original"){
          cm1 <- cbind(1, X_2.mean[,1], A_2*X_2.mean[,2])
          #psi <- c(-1.6, 2.3, -3.2)
        }

        if(size==600  & settingMiss=="Original"){
          cm1 <- cbind(1, X_2.mean[,1], A_2*X_2.mean[,2])
          #psi <- c(-1.6, 2.3, -3.2)
        }

        if(size==1200  & settingMiss=="Original"){
          cm1 <- cbind(1, X_2.mean[,1], A_2*X_2.mean[,2])
          psi <- c(-0.3, -0.4, 1.3)
        }
    }

    missBin <- mapply(function(x){ rbinom(1,1,prob=x)}, x=expit(cm1%*%psi))

    fullDat <- data.frame(missBin, Y, A_2, X_2.mean)
    names(fullDat)[-c(1:3)] <- paste0('X', 1:p)

    # prepare data set with missingness
    missDat <- fullDat
    missDat[which(missDat$missBin==1),]$Y <- NA
    ## exclude the missBin variable
    missDat <- missDat[,c(2:ncol(missDat))]
  }

  if(settingNumber==2){
    if(propMiss==30){
      if(size==300  & settingMiss=="Original"){
        cm1 <- cbind(1, X_2.mean[,1], A_2*X_2.mean[,2])
        #psi <- c(-1.6, 2.3, -3.2)
      }

      if(size==600  & settingMiss=="Original"){
        cm1 <- cbind(1, X_2.mean[,1], A_2*X_2.mean[,2])
        #psi <- c(-1.6, 2.3, -3.2)
      }

      if(size==1200  & settingMiss=="Original"){
        cm1 <- cbind(1, X_2.mean[,1], A_2*X_2.mean[,2])
        psi <- c(-3.8, 2.2, -2.8)

        cm1_cov <- cbind(1, X_2.mean[,1])
        psi_cov <- c(-3, -0.5)


      }
    }

    if(propMiss==50){
      if(size==300  & settingMiss=="Original"){
        cm1 <- cbind(1, X_2.mean[,1], A_2*X_2.mean[,2])
        #psi <- c(-1.6, 2.3, -3.2)
      }

      if(size==600  & settingMiss=="Original"){
        cm1 <- cbind(1, X_2.mean[,1], A_2*X_2.mean[,2])
        #psi <- c(-1.6, 2.3, -3.2)
      }

      if(size==1200  & settingMiss=="Original"){
        cm1 <- cbind(1, X_2.mean[,1], A_2*X_2.mean[,2])
        psi <- c(-0.6, 0.5, -0.7)

        cm1_cov <- cbind(1, X_2.mean[,1])
        psi_cov <- c(-2.5, -0.8)
      }
    }

    missBin <- mapply(function(x){ rbinom(1,1,prob=x)}, x=expit(cm1%*%psi))

    missBin_cov <- mapply(function(x){ rbinom(1,1,prob=x)}, x=expit(cm1_cov%*%psi_cov))


    fullDat <- data.frame(missBin, missBin_cov, Y, A_2, X_2.mean)
    names(fullDat)[-c(1:4)] <- paste0('X', 1:p)

    # prepare data set with missingness
    missDat <- fullDat
    missDat[which(missDat$missBin==1),]$Y <- NA
    missDat[which(missDat$missBin_cov==1),]$X4 <- NA ## missingness in covariate

    ## exclude the missBin and missBin_cov variables
    missDat <- missDat[,c(3:ncol(missDat))]


  }


  if(settingNumber==3){
    p_1 <- p-1

    if(propMiss==30){
      if(size==300  & settingMiss=="Original"){
        cm1 <- cbind(1, X_2.mean[,1], A_2*X_2.mean[,2])
        #psi <- c(-1.6, 2.3, -3.2)
      }

      if(size==600  & settingMiss=="Original"){
        cm1 <- cbind(1, X_2.mean[,1], A_2*X_2.mean[,2])
        #psi <- c(-1.6, 2.3, -3.2)
      }

      if(size==1200  & settingMiss=="Original"){
        cm1 <- cbind(1, X_2.mean[,1], A_2*X_2.mean[,2])
        psi <- c(-1.6, 0.1, 0.2)

        cm1_cov <- cbind(1, X_2.mean[,1])
        psi_cov <- c(-3.7, -0.3)

        missBin_covB <- matrix(rbinom(n*p_1, 1, 0.03),n,p_1)


      }
    }

    if(propMiss==50){
      if(size==300  & settingMiss=="Original"){
        cm1 <- cbind(1, X_2.mean[,1], A_2*X_2.mean[,2])
        #psi <- c(-1.6, 2.3, -3.2)
      }

      if(size==600  & settingMiss=="Original"){
        cm1 <- cbind(1, X_2.mean[,1], A_2*X_2.mean[,2])
        #psi <- c(-1.6, 2.3, -3.2)
      }

      if(size==1200  & settingMiss=="Original"){
        cm1 <- cbind(1, X_2.mean[,1], A_2*X_2.mean[,2])
        psi <- c(-3.8, 2.7, -4)

        cm1_cov <- cbind(1, X_2.mean[,1])
        psi_cov <- c(-0.3, 2.5)

        missBin_covB <- matrix(rbinom(n*p_1, 1, 0.05),n,p_1)
      }
    }

    missBin <- mapply(function(x){ rbinom(1,1,prob=x)}, x=expit(cm1%*%psi))

    missBin_cov <- mapply(function(x){ rbinom(1,1,prob=x)}, x=expit(cm1_cov%*%psi_cov))


    fullDat <- data.frame(missBin, missBin_covB[,1:3], missBin_cov, missBin_covB[,4], Y, A_2, X_2.mean)
    names(fullDat)[-c(1:8)] <- paste0('X', 1:p)

    # prepare data set with missingness
    missDat <- fullDat
    missDat[which(missDat$missBin==1),]$Y <- NA

    ## missingness in covariates
    missList <- as.list(fullDat[,2:6])
    targList <- as.list(fullDat[,-c(1:8)])
    missDat[,-c(1:8)] <- as.data.frame(mapply(function(miss, targ) ifelse(miss==1, NA, targ),missList, targList))

    ## exclude the missBin, missBin_cov, and missBin_covB variables
    missDat <- missDat[,c(7:ncol(missDat))]

  }


  dfsList <- list(fullDat, missDat)


  return(dfsList)


}
