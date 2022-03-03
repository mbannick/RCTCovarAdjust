source("~/repos/RCTCovarAdjust/R/pvalues.R")

# Confidence intervals

#' Get a confidence interval based on the stage-wise ordering
#'
#' @param est Estimate
#' @param sd_K Variance at last stage (or stage of rejection)
#' @param n_K Sample size at the last stage (or stage of rejection)
#' @param u_k Boundaries up until the point of rejection
#'
#' @examples
#' # ------------------------------------------------------ #
#' u_k <- c(4.332634, 2.963132, 2.359044, 2.014090)
#' u_k1 <- cbind(-u_k, u_k)
#' u_k2 <- cbind(rep(-Inf, 4), u_k)
#'
#' # This gives exactly alpha = 0.05 by putting in the last
#' # boundary.
#' n_k1 <- c(10, 20, 30, 40)
#' n_k2 <- c(10, 20)
#'
#' get.confint.sw(est=1, u_k=u_k1[1:3,], sd_K=1,
#'                n_k=n_k1, alpha=0.05,
#'                ancova_monitor=F, ancova_test=F, last_stage=T)
get.confint.sw <- function(est, sd_K, n_k, u_k, k_r, alpha,
                           rho=1, ancova_monitor, ancova_test,
                           last_stage){

  # Get sample size at the last stage
  n_K <- n_k[length(n_k)]

  # Function to translate effect size into z-statistic
  # At the analysis stage K (not at the first stage)
  get.z <- function(eff) sqrt(n_K) * (eff - est) / sd_K
  # Create a function to translate the estimate to a z-statistic
  search.fun <- function(eff, low=FALSE){
    z <- get.z(eff)
    if(low){
      p <- get.pvalue.sw(z, u_k=u_k, n_k=n_k, k_r=k_r, rho=rho,
                         ancova_monitor=ancova_monitor,
                         ancova_test=ancova_test,
                         last_stage=last_stage, type="lower")
    } else {
      p <- get.pvalue.sw(z, u_k=u_k, n_k=n_k, k_r=k_r, rho=rho,
                         ancova_monitor=ancova_monitor,
                         ancova_test=ancova_test,
                         last_stage=last_stage, type="upper")
    }
    return(p - alpha/2)
  }
  # Test out a lot of different mean values to see what the p-value is
  lower <- uniroot(search.fun, low=FALSE, lower=-100, upper=est, trace=1)$root
  upper <- uniroot(search.fun, low=TRUE, lower=est, upper=100, trace=1)$root
  return(c(lower, upper))
}
