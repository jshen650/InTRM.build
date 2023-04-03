#' Estimate the Value and variance of the Value using the AIPW and IPW estimators
#' Based on the formulas provided in the book, "Dynamic Treatment Regimes" (2020) by Tsiatis, Davidian, Holloway, and Laber
#'
#' @param testDf Data frame of interest for estimating the Value and variance of the Value
#' @param coef_dr Vector representing decision rule components in the form
#' of c(intercept, trtmt, main effects, main effects by trtmt interactions)
#' @param isSplit Boolean for whether the Value estimate is obtained on a split data set (rather than full)
#' @param estVar If TRUE, will also estimate the variance of the Value. Otherwise, will just return the AIPW estimate of the Value.
#'
#' @return Either a scalar representing the AIPW estimate of Value or
#' A vector of the estimate of Value and the variance of the Value for the IPW and AIPW estimators
ValEst <- function(testDf, coef_dr, isSplit=TRUE, estVar=TRUE){
  if(isSplit){
    n <- nrow(testDf)/2
  } else{
    n <- nrow(testDf)
  }

  nCov <- ncol(testDf)

  test_df <- as.data.frame(testDf)

  ## assumes Y and A_2 are in the first two columns
  ## with covariates in columns 3 to end
  test_sub <- testDf[,3:ncol(testDf)]

  test_m <- as.matrix(test_sub)
  test_opt <- recTrt_dr(test_m, coef_dr) # get recommended trtmt options under decision rule
  test_C <- which(test_opt==as.numeric(as.character(test_df$A_2))) # get indices for individuals that are consistent
  test_subC <- test_df[test_C,] # subset of individuals who are consistent
  test_Q <- mean(test_subC$Y)

  test_df$C <- ifelse(test_opt==as.numeric(as.character(test_df$A_2)), 1, 0)
  test_ps <- rep(0.5, times=nrow(test_df))
  test_ind <- ifelse(test_opt==-1, 1, 0)
  test_pid <- ifelse(test_ind==1, test_ps, (1-test_ps) )

  ## IPW estimate
  test_ps_ipw <- mean(test_df$C * test_df$Y/test_pid)

  ## IPW variance
  test_An_IPW <- (test_df$C * test_df$Y/test_pid) - test_ps_ipw
  test_A_IPW <- mean(test_An_IPW^2)

  ## AIPW estimate
  test_ps_aipw <- test_ps_ipw - mean( ( (test_df$C - test_pid)/(test_pid) ) * test_Q )

  if(estVar==FALSE){
    return(test_ps_aipw)
  }


  ## AIPW variance
  test_An <- (test_df$C * test_df$Y/test_pid) - (( (test_df$C - test_ps)/(test_pid) ) * test_Q) - test_ps_aipw

  test_A <- mean(test_An^2)


  covariates <- colnames(test_C)[-c(1:2)] # get the Xs

  test_mm <- cbind(1, test_df[,-c(1, ncol(test_df))])
  test_mm$A_2 <- as.numeric(as.character(test_mm$A_2))
  for(i in covariates){
    test_mm[[paste0("A.", i)]] <- with(test_mm, test_mm$A_2*test_mm[,i])
  }

  test_dQd <- as.matrix(test_mm)
  test_part <- as.vector((test_df$C - test_pid)/(test_pid))
  test_D <- colMeans(-1*test_part * test_dQd)
  test_G <-as.vector(test_df$C*{test_df$Y - test_Q}/test_pid^2*{-1}^{test_opt})

  mainCov <- colnames(test_mm)[-c(1:2)]

  test_ps_mm <- model.matrix( glm(as.formula(paste("A_2~", paste(mainCov, collapse="+"))), data = test_df, family = "binomial"))

  test_G <- colMeans(x=test_G*test_ps*(1-test_ps)*test_ps_mm)

  test_C_IPW <- colMeans(x=as.vector( test_df$C*test_df$Y / test_pid * (-1*test_opt + test_ps)) * test_ps_mm)

  test_p <- ncol(x=test_dQd)

  test_Bpart <- (test_df$Y - test_Q)*test_dQd

  test_BBt <- sapply(X=1L:n,
                     FUN= function(i){
                       test_Bpart[i,] %o% test_Bpart[i,]},
                     simplify="array")
  test_B <- rowMeans(test_BBt, dim=2L)
  test_En <- rowMeans(sapply(X=1L:n,
                             FUN= function(i){
                               test_dQd[i,] %o% test_dQd[i,]},
                             simplify="array"), dim=2L)
  test_inv_En <- solve(test_En)
  test_Cpart <- ((as.numeric(as.character(test_df$A_2))- test_ps)/(test_ps *(1-test_ps))) * test_ps_mm

  test_C <- rowMeans( sapply(X = 1L:n,
                             FUN = function(i){ test_Cpart[i,] %o% test_Cpart[i,] },
                             simplify = "array"), dim=2L)

  test_J <-  rowMeans(sapply(X=1L:n,
                             FUN= function(i){
                               (test_ps_mm[i,] %o% test_ps_mm[i,])/(test_ps[i] * (1-test_ps[i])) },
                             simplify="array"), dim=2L)

  test_inv_J <- solve(test_J)


  test_Dbread <- t(test_D)%*%(test_inv_En)
  test_DBD <- test_Dbread%*%test_B %*%t(test_Dbread)

  test_varHat = ((test_A + test_DBD) / nrow(test_df))
  test_varHat_IPW = ((test_A_IPW)/nrow(test_df))

  ## return IPW, AIPW, Var (IPW), Var (AIPW)
  return(c(test_ps_ipw, test_ps_aipw, test_varHat_IPW, test_varHat))

}
