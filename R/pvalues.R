#' P-values
library(magrittr)
library(MASS)
library(mvtnorm)
source("R/trial-data.R")
source("R/trial-funcs.R")

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
#' corr <- basic.cov(n_k)
#' get.pvalue.sw(obs=u_k[4], u_k=u_k[1:3], corr=corr)
#' get.pvalue.sw(obs=2.5, u_k=u_k[1:3], corr=corr)
#' get.pvalue.sw(obs=1.5, u_k=u_k[1:3], corr=corr)
#' get.pvalue.sw(obs=-1.5, u_k=u_k[1:3], corr=corr)
get.pvalue.sw <- function(obs, u_k, corr, mean=NULL,
                          algorithm=Miwa(steps=1000)){

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

      p.upper <- pnorm(u_k[i], mean=mean[i], lower.tail=F)
      p.lower <- pnorm(u_k[i], mean=mean[i], lower.tail=F)

    } else {
      corr.i <- corr[1:i, 1:i]

      lower.prev <- -u_k[1:(i-1)]
      upper.prev <- u_k[1:(i-1)]

      if(i == K){
        lower.i <- -abs(obs)
        upper.i <- -abs(obs)
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
  return(2*min(p.upper.tot, p.lower.tot))
}
