#' Double bootstrap for selecting alpha to use with the m-out-of-n bootstrap
#'
#' @param eta Scalar for the (1-eta)% confidence set
#' @param estVal Scalar for the estimate of Value that the confidence set is meant to surround
#' @param i Scalar representing a seed or iteration
#' @param B2 Scalar for number of conditional bootstraps to run
#' @param see_df_list List of imputed data sets
#' @param coef_dr Vector representing decision rule components in the form
#'  of c(intercept, trtmt, main effects, main effects by trtmt interactions)
#' @param current_alpha Scalar for value of alpha used
#'
#' @return Scalar indicator for whether the double centered percentile bootstrap captures the estimate of Value
double_boot <- function(i, eta, estVal, B2, see_df_list, coef_dr, current_alpha){
  set.seed(i)
  n <- nrow(see_df_list[[1]])

  s1_resample_rows = sample(x=nrow(see_df_list[[1]]), size=n, replace=TRUE)

  ## get estimates of m for each of the imputed data sets
  run_get_m <- lapply(see_df_list, function(x){
    s1_samp <- x[s1_resample_rows,]
    s1_m <- get_m_alpha(s1_samp, alpha=current_alpha, coef_dr=coef_dr)
    return(s1_m)
  })

  ## choose the minimum recommended m across imputations
  s1_m_use = floor(min(unlist(run_get_m)))

  ## run bootstrap conditional on initial bootstrap
  run_double = lapply(see_df_list, function(x){
    s1_samp <- x[s1_resample_rows,]
    s1_mVal <- mean( ValEst(s1_samp, coef_dr, isSplit = FALSE, estVar=FALSE) )

    s2_mVal <- mboot_Val_multi(m=s1_m_use, see_df=s1_samp, coef_dr=coef_dr, reps=B2)
    mboot_CS(eta=eta, m=s1_m_use, mboot_multi_res = s2_mVal, estVal=mean( ValEst(x, coef_dr, isSplit = FALSE, estVar=FALSE) ))
  })


  dcpb_CS <- colMeans( t(as.data.frame(run_double) ))

  boot_indicator <- ifelse(dcpb_CS[1] <= estVal & dcpb_CS[2] >= estVal, 1, 0)
  return(boot_indicator)
}
