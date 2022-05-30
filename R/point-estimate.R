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
#' u_k1 <- cbind(-u_k, u_k)
#' u_k2 <- cbind(rep(-Inf, 4), u_k)
#'
#' # This gives exactly alpha = 0.05 by putting in the last
#' # boundary.
#' n_k1 <- c(10, 20, 30, 40)
#' n_k2 <- c(10, 20)
#'
#' get.point.sw(est=1, u_k=u_k1[1:3,], sd_K=0.2,
#'              n_k=n_k1, alpha=0.05,
#'              ancova_monitor=F, ancova_test=F, last_stage=T)
#'
#' get.point.sw(est=1, u_k=matrix(u_k1[1,], nrow=1), sd_K=1,
#'              n_k=c(10), alpha=0.05,
#'              ancova_monitor=F, ancova_test=F, last_stage=F)
get.point.sw <- function(est, sd_K, n_k, ...){
  val <- search.fun.sw(est, sd_K, n_k, alpha=0.5, low=TRUE, ...)
  # fun <- function(est) search.fun.sw(est, sd_K, n_k, alpha=0.5, low=TRUE, ...)
  # effs <- seq(-abs(est)*3, abs(est)*3, by=0.05)
  # lower.vec <- sapply(effs, fun)
  # plot(lower.vec ~ effs, type='l', ylab="p-value", xlab="effect size")
  return(val)
}
