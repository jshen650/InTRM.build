#' Run the double bootstrap procedure for selecting alpha to use with the m-out-of-n bootstrap
#'
#' @param missDat
#' @param numImp
#' @param alpha_see
#' @param seedNum
#' @param iters
#' @param B2
#' @param eta
#' @param runLPC
#'
#' @import mice
#'
#' @return Vector for the recommended alpha and the corresponding coverage
run_double_boot <- function(missDat, numImp, alpha_see, seedNum, iters, B2, eta=0.05, runLPC){

  if(runLPC==TRUE){
    nCores = as.numeric(Sys.getenv('LSB_DJOB_NUMPROC'))
  } else{
    nCores = runLPC
  }


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


  # form list of imputations then form training and seeing lists
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

  ## candidate values for alpha
  alpha_see <- seq(from= 0.025, to=1 ,by=0.025)

  db_Res <- c(0.025, 0)
  for (aVal in alpha_see){
    run_double_boot = mclapply(1:iters, double_boot, eta=eta,estVal=ValRes, B2=B2, see_df_list=df_list,
                               coef_dr=coef_dr, current_alpha=aVal, mc.cores = nCores )
    get_coverage = mean(do.call(rbind, run_double_boot))
    if(get_coverage >= 0.95){
      db_Res <- c(aVal, get_coverage)
      break
    }
  }

  ## output
  db_Res <- t(as.data.frame(db_Res))
  colnames(db_Res) <- c("alpha", "coverage")

  return(db_Res)

}
