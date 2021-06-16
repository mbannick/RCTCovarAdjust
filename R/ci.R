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
#' u_k <- c(4.332634, 2.963132, 2.359044, 2.014090)
#' K <- length(u_k)
#' n_k <- rep(10, K)
#' N <- sum(n_k)
#' corr.1 <- basic.cov(n_k)
#' corr.2 <- basic.cov(n_k, rho=0.9)
#' corr.3 <- basic.cov(n_k, rho=0.9, extra=TRUE)
#' alpha <- 0.05
#' get.confint.sw(est=0.32, sd_K=1, n_K=sum(n_k), corr=corr.1, alpha=alpha, u_k=u_k[1:(K-1)])
#' get.confint.sw(est=0.32, sd_K=1, n_K=sum(n_k), corr=corr.2, alpha=alpha, u_k=u_k[1:(K-1)])
#'
#' # Monitor with ANOVA and conduct ANCOVA after final ANOVA
#' get.confint.sw(est=0.32, sd_K=1, n_K=sum(n_k), corr=corr.3, alpha=alpha, u_k=u_k)
#'
#' # ------------------------------------------------------ #
#' # FIX IT WITH NEW BOUNDS
#' # These are the alpha spending bounds for alpha = 0.045, so that there's
#' # some type I error left over.
#'
#' u_k2 <- c(4.415989, 3.023536, 2.408191, 2.056245)
#' get.confint.sw(est=0.32, sd_K=1, n_K=sum(n_k), corr=corr.3, alpha=alpha, u_k=u_k2)
get.confint.sw <- function(est, sd_K, n_K, u_k, corr,
                           alpha, algorithm=Miwa(steps=1000)){
  # Function to translate effect size into z-statistic
  # At the analysis stage K (not at the first stage)
  get.z <- function(eff) sqrt(n_K) * (eff - est) / sd_K

  # Create a function to translate the estimate to a z-statistic
  search.fun <- function(eff, low=FALSE){
    z <- get.z(eff)
    if(low){
      p <- get.pvalue.sw(z, u_k=u_k, corr=corr, type="lower")
    } else {
      p <- get.pvalue.sw(z, u_k=u_k, corr=corr, type="upper")
    }
    return(p - alpha/2)
  }
  # Test out a lot of different mean values to see what the p-value is
  lower <- uniroot(search.fun, low=FALSE, lower=-100, upper=est, trace=1)$root
  upper <- uniroot(search.fun, low=TRUE, lower=est, upper=100, trace=1)$root
  return(c(lower, upper))
}
