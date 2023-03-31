#' Confidence Set of m-out-of-n bootstrap results
#'
#' @param eta Scalar for the (1-eta)\% confidence set
#' @param m Scalar for subsample size, where m < n
#' @param mboot_multi_res Vector of results from multiple repetitions of the m-out-of-n bootstrap
#' @param estVal Scalar for estimate of Value that the confidence set surrounds
#'
#' @return Vector representing the lower and upper bounds of the (1-eta)% confidence set
mboot_CS <- function(eta, m, mboot_multi_res, estVal){
  lower_p <- eta/2
  upper_p <- 1 - lower_p

  ## get l_hat and u_hat, as defined in the paper from Chakraborty et al. (2013)
  mboot_get <- sqrt(m) * (mboot_multi_res - estVal)
  mboot_quant <- quantile(mboot_get, probs=c(lower_p, upper_p))
  l_hat <- mboot_quant[1]
  u_hat <- mboot_quant[2]

  return( c(estVal - u_hat/sqrt(m), estVal - l_hat/sqrt(m) ) )

}
