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
#'                ancova_monitor=F, ancova_test=F, last_stage=T,
#'                crossed_lower=FALSE)
get.confint.sw <- function(est, sd_K, n_k, alpha, ...){

  lower <- search.fun.sw(est, sd_K, n_k, alpha=alpha/2, low=FALSE, ...)
  upper <- search.fun.sw(est, sd_K, n_k, alpha=alpha/2, low=TRUE, ...)

  # effs <- seq(-abs(est)*3, abs(est)*3, by=0.02)
  # lower.vec <- sapply(effs, search.fun, low=FALSE) + alpha/2
  # upper.vec <- sapply(effs, search.fun, low=TRUE) + alpha/2
  #
  # plot(lower.vec ~ effs, type='l', ylab="p-value", xlab="effect size")
  # lines(upper.vec ~ effs, col='blue')
  # legend(x=1.5, y=0.8, legend=c("lower", "upper"), col=c("black", "blue"),
  #        lty=c(1,1), cex=0.5)
  # abline(h=0.975, lty='dashed')
  # abline(h=0.025, lty='dashed')

  return(c(lower, upper))
}
