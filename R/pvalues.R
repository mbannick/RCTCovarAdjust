library(magrittr)
library(MASS)
library(mvtnorm)
source("~/repos/RCTCovarAdjust/R/covariance.R")

#' Get a stage-wise p-value
#'
#' @param obs Observed z-statistic
#' @param u_k Vector of boundaries for *previous* stages
#' @param corr Correlation matrix for test statistics
#' @param mean Mean shift vector for the test statistic
#' @param ... Additional arguments for pmvnorm algorithm
#'
#' @examples
#' # OBF boundaries
#' u_k <- c(4.332634, 2.963132, 2.359044, 2.014090)
#'
#' # This gives exactly alpha = 0.05 by putting in the last
#' # boundary.
#' n_k <- rep(100, 4)
#' corr.1 <- basic.cov(n_k)
#' corr.2 <- basic.cov(n_k, rho=0.5)
#' get.pvalue.sw(obs=u_k[4], u_k=u_k[1:3], corr=corr.1)
#' get.pvalue.sw(obs=2.5, u_k=u_k[1:3], corr=corr.1)
#' get.pvalue.sw(obs=1.5, u_k=u_k[1:3], corr=corr.1)
#' get.pvalue.sw(obs=-1.5, u_k=u_k[1:3], corr=corr.1)
#' get.pvalue.sw(obs=-1.5, u_k=u_k[1:3], corr=corr.2)
get.pvalue.sw <- function(obs, u_k, corr, mean=NULL,
                          algorithm=Miwa(steps=1000),
                          type="two-sided"){

  if(!type %in% c("two-sided", "lower", "upper")) stop(
    "Unrecognized p-value type. ",
    "Provide one of two-sided, lower, or upper."
  )

  K <- length(u_k) + 1

  if(is.null(mean)){
    mean <- rep(0, K)
  } else {
    if(length(mean) != K) stop()
  }

  p.upper.tot <- 0
  p.lower.tot <- 0

  for(i in 1:K){
    if(i == 1){

      if(length(u_k) == 0){ # reject at first stage
        tail <- pnorm(abs(obs), mean=mean[i], lower.tail=F)
      } else {
        tail <- pnorm(u_k[i], mean=mean[i], lower.tail=F)
      }

      p.upper <- tail
      p.lower <- tail

    } else {
      corr.i <- corr[1:i, 1:i]

      lower.prev <- -u_k[1:(i-1)]
      upper.prev <- u_k[1:(i-1)]

      if(i == K){
        lower.i <- -abs(obs)
        upper.i <- abs(obs)
      } else {
        lower.i <- -u_k[i]
        upper.i <- u_k[i]
      }

      p.upper <- pmvnorm(
        lower=c(lower.prev, upper.i),
        upper=c(upper.prev, Inf),
        corr=corr.i,
        mean=mean[1:i],
        algorithm=algorithm
      )
      p.lower <- pmvnorm(
        lower=c(lower.prev, -Inf),
        upper=c(upper.prev, lower.i),
        corr=corr.i,
        mean=mean[1:i],
        algorithm=algorithm
      )
    }
    p.upper.tot <- p.upper.tot + p.upper
    p.lower.tot <- p.lower.tot + p.lower
  }

  if(type == "two-sided"){
    return(2*min(p.upper.tot, p.lower.tot))
  } else if(type == "lower"){
    return(p.upper.tot[1])
  } else if(type == "upper"){
    return(p.lower.tot[1])
  }
}
