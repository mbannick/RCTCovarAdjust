source("~/repos/RCTCovarAdjust/R/pvalues.R")
source("~/repos/RCTCovarAdjust/R/covariance.R")

# Median Unbiased Point Estimation

#' Get a median unbiased point estimate based on the stage-wise ordering
#'
#' @param est Estimate
#' @param sd_K Variance at last stage (or stage of rejection)
#' @param n_K Sample size at the last stage (or stage of rejection)
#' @param u_k Boundaries up until the point of rejection
#'
#' @examples
#' u_k <- c(4.332634, 2.963132, 2.359044, 2.014090)
#' u_k <- cbind(-u_k, u_k)
#' K <- 4
#' n_k <- cumsum(rep(10, K))
#' N <- sum(n_k)
#' corr.1 <- corr.mat(n_k)
#' corr.2 <- corr.mat(n_k, rho=0.9, mis=c(F, F, F, T))
#' corr.3 <- corr.mat(c(n_k, n_k[4]), rho=0.9, mis=c(F, F, F, F, T))
#' alpha <- 0.05
#' get.point.sw(est=0.32, sd_K=1, n_K=sum(n_k), corr=corr.1, alpha=alpha, u_k=u_k[1:(K-1),])
#' get.point.sw(est=0.32, sd_K=1, n_K=sum(n_k), corr=corr.2, alpha=alpha, u_k=u_k[1:(K-1),])
get.point.sw <- function(est, sd_K, n_K, u_k, corr, alpha){
  # Function to translate effect size into z-statistic
  # At the analysis stage K (not at the first stage)
  get.z <- function(eff) sqrt(n_K) * (eff - est) / sd_K

  # Create a function to translate the estimate to a z-statistic
  search.fun <- function(eff){
    z <- get.z(eff)
    p <- get.pvalue.sw(z, u_k=u_k, corr=corr, type="lower")
    return(p - 0.5)
  }
  # Test out a lot of different mean values to see what the p-value is
  est <- uniroot(search.fun, lower=-100, upper=100, trace=1)$root
  return(est)
}
