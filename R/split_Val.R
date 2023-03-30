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

  ValEst_res <- lapply(test_list, function(x){
    ValEst(x, coef_dr)
  })

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
