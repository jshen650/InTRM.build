#' Estimate of m for the m-out-of-n bootstrap for a given value of alpha
#'
#' @param see_df Data frame of interest
#' @param alpha Scalar for the alpha value used in estimating m for the m-out-of-n bootstrap
#'
#' @return Scalar estimate of m
get_m_alpha <- function(see_df, alpha, coef_dr){
  intercept_dr <- coef_dr[1]
  trtcoef_dr <- coef_dr[2]
  covcoef_dr <- coef_dr[3:7]
  intcoef_dr <- coef_dr[8:12]

  n_see <- nrow(see_df)

  see_X2 <- model.matrix(lm(Y ~ .*A_2, data=see_df))
  see_Sigma2 <- n_see*solve(t(see_X2)%*%see_X2)
  see_Z2 <- diag(as.vector(see_df$Y - see_X2%*%coef_dr)) %*% see_X2 %*% see_Sigma2/sqrt(n_see-dim(see_X2)[2])
  see_Cov2 <- t(see_Z2)%*%see_Z2

  see_Sigma_stage2 <- see_Cov2[c(2:7),c(2:7)] ## the treatment and history variables
  see_H <- as.matrix(cbind(1, see_X2[,-c(1:7)])) ## patient history variables

  see_sigma2 <- diag(see_H%*%see_Sigma_stage2%*%t(see_H))/n_see
  see_Yprime <- trtcoef_dr + see_X2[,-c(1:7)]%*%intcoef_dr

  see_TS <- abs(see_Yprime)/sqrt(see_sigma2)
  presee_lev <- 0.001            # presee level (\nu, in paper)
  see_cutoff <- qnorm(1 - presee_lev/2)
  see_nonregularity <- (see_TS <= see_cutoff )

  see_p <- mean(see_nonregularity)
  see1_alpha <- alpha
  see1_m <- ceiling ( n_see^{ (1 + see1_alpha*(1-see_p))/(1+see1_alpha)} )
  return(see1_m)

}
