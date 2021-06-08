source("R/trial-data.R")
source("R/trial-funcs.R")
source("R/pvalues.R")

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
#' corr.2 <- basic.cov(n_k, rho=0.1)
#' alpha <- 0.05
#' get.confint.sw(est=0.32, sd_K=1, n_K=sum(n_k), corr=corr.1, alpha=alpha, u_k=u_k[1:3])
#' get.confint.sw(est=0.32, sd_K=1, n_K=sum(n_k), corr=corr.2, alpha=alpha, u_k=u_k[1:3])
get.confint.sw <- function(est, sd_K, n_K, u_k, corr, alpha, algorithm=Miwa(steps=1000)){
  K <- length(u_k) + 1
  get.z <- function(eff) sqrt(n_K) * (eff - est) / sd_K

  # Create a function to translate the estimate to a z-statistic
  search.fun <- function(eff, low=FALSE){
    z <- get.z(eff)
    if(low){
      p <- get.pvalue.sw(z, u_k=u_k[1:(K-1)], corr=corr, type="lower")
    } else {
      p <- get.pvalue.sw(z, u_k=u_k[1:(K-1)], corr=corr, type="upper")
    }
    return(p - alpha/2)
  }

  # Test out a lot of different mean values to see what the p-value is
  lower <- uniroot(search.fun, low=FALSE, lower=-100, upper=est)$root
  upper <- uniroot(search.fun, low=TRUE, lower=est, upper=100)$root

  return(c(lower, upper))
}

