#' Estimate Value using the data splitting approach on data sets with missingness
#'
#' @param missDat Data set containing missingness
#' @param numImp Number of imputations for imputing with MICE
#'
#' @import mice
#' @return Estimates of the Value and variance of the Value under the IPW and AIPW estimators
#'
split_Val <- function(missDat, numImp, seedNum){
  n <- nrow(missDat)/2

  ## assumes training data set comprises the first n rows of the data set
  train_ind <- c(1:n)

  train_dat <- missDat[train_ind, ]
  test_dat <- missDat[-train_ind,]

  ##############################
  ## Obtain imputed data sets ##
  ##############################

  covariates <- colnames(missDat)[-c(1:2)] # get the Xs

  for(i in covariates){
    train_dat[[paste0("A.", i)]] <- with(train_dat, train_dat$A_2*train_dat[,i])
    test_dat[[paste0("A.", i)]] <- with(test_dat, test_dat$A_2*test_dat[,i])

  }

  train_dat[[paste0("Y.A")]] <- with(train_dat, train_dat$Y*train_dat$A_2)
  test_dat[[paste0("Y.A")]] <- with(test_dat, test_dat$Y*test_dat$A_2)

  meth <- make.method(train_dat)
  pred <- make.predictorMatrix(train_dat)
  for(i in covariates){
    interTerm <- paste0("A.", i)
    meth[interTerm] <- paste("~I(A_2*", i, ")", sep="")
    pred[c('A_2', i),interTerm] <- 0
  }

  pred[c('Y', 'A_2'), "Y.A"] <- 0

  imputed_res_train <- mice(train_dat, m=numImp, meth=meth, pred=pred, seed=seedNum, print=FALSE)
  imputed_res_test <- mice(test_dat, m=numImp, meth=meth, pred=pred, seed=seedNum, print=FALSE)

  train_list <- lapply(1:numImp, function(i){
    df <- complete(imputed_res_train, i)[,-c( (ncol(train_dat)-5):(ncol(train_dat)))] # exclude the derived columns so can continue w older code
    df$A_2 <- as.factor(df$A_2)
    df
  })

  test_list <- lapply(1:numImp, function(i){
    df <- complete(imputed_res_test, i)[,-c( (ncol(test_dat)-5):(ncol(test_dat)))] # exclude the derived columns so can continue w older code
    df$A_2 <- as.factor(df$A_2)
    df
  })

  ##############################
  ## Training ##
  ##############################

  train_mod <- lapply(train_list, function(x){
    mod <- lm(Y ~ .*A_2, data=x)
    coefmod <- as.data.frame(coef(mod))
    coefmod <- tibble::rownames_to_column(coefmod, "Covariate")
    return(coefmod)
  })

  ##############################
  ## Model-averaged decision rule ##
  ##############################

  ## gather all coefficients and average across
  coef_dr <- rowMeans(sapply(train_mod, "[[", 2))

  ## gather components of decision rule
  intercept_dr <- coef_dr[1]
  trtcoef_dr <- coef_dr[2]
  covcoef_dr <- coef_dr[3:7]
  intcoef_dr <- coef_dr[8:12]

  ## return vector of recommended trtmts under decision rule
  giveOpt_dhat.mod <- function(matCov){
    assign1 <- 1
    assign2 <- -1

    # consider outcome with trtmt assignment == 1
    Y1 <- intercept_dr + matCov%*%covcoef_dr + assign1*(trtcoef_dr + matCov%*%intcoef_dr)

    # consider outcome with trtmt assignment == -1
    Y2 <- intercept_dr + matCov%*%covcoef_dr

    out <- ifelse(Y1 > Y2, assign1, assign2)

    return(out)
  }

  #######################################
  ## Get Value and Variance estimate of Value for IPW and AIPW estimators ##
  ## Reference: Dynamic Treatment Regimes - Tsiatis, Davidian, Holloway, and Laber
  #######################################

  ValEst <- function(testDf){
    test_df <- as.data.frame(testDf)
    test_sub <- subset.data.frame(testDf, select=c("X1", "X2", "X3", "X4", "X5"))
    test_m <- as.matrix(test_sub)
    test_opt <- giveOpt_dhat.mod(test_m) # get recommended trtmt options under decision rule
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
    test_ps_mm <- model.matrix( glm(A_2 ~ X1 + X2 + X3 + X4 + X5, data = test_df, family = "binomial"))

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

  ValEst_res <- lapply(test_list, ValEst)

  ## IPW Value, averaged acros imputations
  Val_proIPW <- mean(unlist(lapply(ValEst_res, function(x) x[1])))

  ## AIPW Value, averaged across imputations
  Val_AIPW <- mean(unlist(lapply(ValEst_res, function(x) x[2])))

  ## AIPW Variance, across imputations

  # average within-imputation Variance
  w_Var_AIPW_imp <- mean(unlist(lapply(ValEst_res, function(x) x[4])))

  # between-imputation variance
  Val_AIPWimp_sqDev <- (unlist(lapply(ValEst_res, function(x) x[2])) - Val_AIPW)^2
  bt_Var_AIPW <- (1/(numImp-1) )*sum(Val_AIPWimp_sqDev)

  ## combine via Rubin's Rules
  Val_AIPW_var <- w_Var_AIPW_imp + bt_Var_AIPW + bt_Var_AIPW/numImp


  ## IPW Variance (following propensity score model estimation), across imputations
  # average within-imputation Variance
  w_Var_v2IPW_imp <- mean(unlist(lapply(ValEst_res, function(x) x[3])))

  # between-imputation variance
  Val_v2IPWimp_sqDev <- (unlist(lapply(ValEst_res, function(x) x[1])) - Val_proIPW)^2
  bt_Var_v2IPW <- (1/(numImp-1) )*sum(Val_v2IPWimp_sqDev)

  ## combine via Rubin's Rules
  Val_v2IPW_var <- w_Var_v2IPW_imp + bt_Var_v2IPW + bt_Var_v2IPW/numImp

  ##############################
  ## Return Output ##
  ##############################

  out1 <- data.frame(Val_proIPW, Val_v2IPW_var , Val_AIPW, Val_AIPW_var)
  colnames(out1) <- c("Value_IPW","Var_IPW" , "Value_AIPW", "Var_AIPW")
  return(out1)

}
