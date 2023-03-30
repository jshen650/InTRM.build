#' Perform multiple repetitions of estimating Value on different subsamples of size m
#'
#' @param m Scalar for subsample size, where m < n
#' @param see_df Data set used for estimating Value
#' @param coef_dr Vector representing decision rule components in the form
#'  of c(intercept, trtmt, main effects, main effects by trtmt interactions)
#' @param reps Scalar for number of repetitions to perform
#'
#' @return Vector of estimates of Value (AIPW) on different subsamples of size m
mboot_Val_multi <- function(m, see_df, coef_dr, reps){
  store_res <- vector(mode="double", length=reps)
  for(i in 1:reps){
    mboot_Val_res <- mboot_Val(m, see_df, coef_dr)
    store_res[i] <- mboot_Val_res
  }
  return(store_res)
}
