#' Estimate AIPW estimate of Value for a subsample of size m on a given data set
#'
#' @param m Scalar for subsample size, where m < n
#' @param see_df Data set used for estimating Value
#' @param coef_dr Vector representing decision rule components in the form
#'  of c(intercept, trtmt, main effects, main effects by trtmt interactions)
#'
#' @return Scalar AIPW estimate of Value on the subsample of size m
mboot_Val <- function(m, see_df, coef_dr){

  resample_rows <- sample(x=nrow(see_df), size= m, replace = TRUE)
  m_samp <- see_df[resample_rows, ]

  m_Val <- mean( ValEst(m_samp, coef_dr, isSplit = FALSE, estVar=FALSE) )
  return(m_Val)

}
