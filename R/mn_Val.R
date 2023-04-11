#' Run the m-out-of-n bootstrap
#'
#' @param missDat Data frame containing missingness
#' @param numImp Scalar for number of imputations
#' @param seedNum Scalar for seed
#' @param reps Scalar for number of m-out-of-n bootstrap repetitions to perform
#' @param eta Scalar for the (1-eta)\% confidence set
#' @param alpha Scalar for the alpha value to use when estimating m
#'
#' @import mice
#'
#' @return List containing the 95% CI formed with the m-out-of-n bootstraps from R imputations,
#' the choice of 'm', and the final 95% CI with m-out-of-n bootstrap
mn_Val <- function(missDat, numImp, seedNum, reps, eta=0.05, alpha=0.025){

  ##############################
  ## Obtain imputed data sets ##
  ##############################
  covariates <- colnames(missDat)[-c(1:2)] # get the Xs

  for(i in covariates){
    missDat[[paste0("A.", i)]] <- with(missDat, missDat$A_2*missDat[,i])

  }

  missDat[[paste0("Y.A")]] <- with(missDat, missDat$Y*missDat$A_2)

  meth <- make.method(missDat)
  pred <- make.predictorMatrix(missDat)
  for(i in covariates){
    interTerm <- paste0("A.", i)
    meth[interTerm] <- paste("~I(A_2*", i, ")", sep="")
    pred[c('A_2', i),interTerm] <- 0
  }

  pred[c('Y', 'A_2'), "Y.A"] <- 0

  imputed_res <- mice(missDat, m=numImp, meth=meth, pred=pred, seed=seedNum, print=FALSE)


  # form list of imputations then form training and testing lists
  ind_imp <- lapply(1:numImp, function(i){
    df <- complete(imputed_res, i)[,-c( (ncol(missDat)-5):(ncol(missDat)))] # exclude the derived columns so can continue w older code
    df$A_2 <- as.factor(df$A_2)
    df
  })


  set.seed(seedNum)
  df_list <- ind_imp

  ##############################
  ## Training ##
  ##############################

  train_mod <- lapply(df_list, function(x){
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

  Val_imp_list <- lapply(df_list, function(x){
    ValEst(x, coef_dr=coef_dr, isSplit=FALSE, estVar=FALSE)
  })


  ValRes <- mean (unlist (lapply (Val_imp_list, mean)))

  ##############################
  ## Run the m-out-of-n bootstrap using the minimum 'm' across imputations ##
  ##############################

  check_get_m <- lapply(df_list, function(x){
    get_m_alpha(x, alpha=alpha, coef_dr=coef_dr)})

  ## use the minimum value of m across imputations
  check_get_m <- ceiling( min(unlist(check_get_m)) )

  ## run multiple repetitions of the m-out-of-n bootstrap
  extend_mboot_Val_multi <- lapply(df_list, function(x){
    mboot_Val_multi(m=check_get_m, see_df=x, coef_dr=coef_dr ,reps=reps)}) ## results in reps x numImp df

  extend_mboot_Val_multi <- t(do.call(rbind, extend_mboot_Val_multi))

  ## create (1-eta)% confidence interval results for the m-out-of-n bootstrap
  ## for each imputation
  extend_mboot_CS <- lapply(seq_along(c(1:length(df_list))), function(x){
    get_CS <- mboot_CS(eta=eta, m=check_get_m, mboot_multi_res = extend_mboot_Val_multi[,x], estVal=mean(Val_imp_list[[x]]) )
    return(get_CS)
  })

  impVar_mBoot <- t(as.data.frame(extend_mboot_CS))

  ## lower and upper bounds for the 95% CI with the m-out-of-n bootstrap constructed on each imputed data set
  out1 <- impVar_mBoot
  colnames(out1) <- c("0.025", "0.975")

  ## Final 95% confidence interval with the m-out-of-n bootstrap, averaging the lower and upper bounds across imputations
  mBootCS_avg <- colMeans(out1)

  ## return list of all 95% confidence intervals for the m-out-of-n bootstraps on each imputed data set,
  ## the choice of m, which is the minimum m across imputations
  ## and the final averaged 95% confidence interval from the m-out-of-n bootstrap
  m_opt <- as.data.frame(check_get_m)
  colnames(m_opt) <- "m_opt"
  dfsList <- list(CI_all = out1, m_opt = m_opt, Value_AIPW = ValRes , CI_fin = as.data.frame(t(mBootCS_avg)))
  return(dfsList)

}
