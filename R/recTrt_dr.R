#' Return vector of recommended treatments under decision rule
#'
#' @param matCov Matrix of covariate information
#' @param coef_dr Vector representing decision rule components in the form
#'  of c(intercept, trtmt, main effects, main effects by trtmt interactions)
#'
#' @return Vector of recommended treatments
#'
recTrt_dr <- function(matCov, coef_dr){
  ## gather components of decision rule

  cov_length <- ncol(matCov)
  lastcov_index <- 2+cov_length
  firstint_index <- 3+cov_length
  lastint_index <- firstint_index+cov_length-1

  intercept_dr <- coef_dr[1]
  trtcoef_dr <- coef_dr[2]
  covcoef_dr <- coef_dr[3:lastcov_index]
  intcoef_dr <- coef_dr[firstint_index:lastint_index]

  assign1 <- 1
  assign2 <- -1

  # consider outcome with trtmt assignment == 1
  Y1 <- intercept_dr + matCov%*%covcoef_dr + assign1*(trtcoef_dr + matCov%*%intcoef_dr)

  # consider outcome with trtmt assignment == -1
  Y2 <- intercept_dr + matCov%*%covcoef_dr

  out <- ifelse(Y1 > Y2, assign1, assign2)

  return(out)

}
