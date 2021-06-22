library(magrittr)
library(MASS)
library(mvtnorm)
source("~/repos/RCTCovarAdjust/R/covariance.R")

.check.type <- function(type){
  if(!type %in% c("two-sided", "lower", "upper")) stop(
    "Unrecognized p-value type. ",
    "Provide one of two-sided, lower, or upper."
  )
}

#' Get the probability of rejecting at stage k.
#' We need this to return a two-element vector because
#' for the orderings, crossing lower *or* upper boundaries
#' is important.
#'
#' @param u_k Vector of boundaries up through stage k
#' @param corr Correlation matrix for test statistics
#' @param mean Mean shift vector for the test statistic
#' @param algorithm Algorithm argument for pmvnorm
#'
#' @examples
#' u_k <- c(4.332634, 2.963132, 2.359044)
#' u_k1 <- cbind(-u_k, u_k)
#' u_k2 <- cbind(rep(-Inf, 3), u_k)
#' n_k <- c(10, 20, 30)
#' corr.1 <- corr.mat(n_k)
#' corr.2 <- corr.mat(n_k, rho=0.5)
#' .reject.prob.k(u_k=u_k1, mean=c(0, 0, 0), corr=corr.1)
#' .reject.prob.k(u_k=u_k2, mean=c(0, 0, 0), corr=corr.1)
.reject.prob.k <- function(u_k, corr, mean=NULL, algorithm=Miwa(steps=1000)){

  if(!"matrix" %in% class(u_k)) u_k <- matrix(u_k, nrow=1)

  K <- nrow(u_k)

  if(K == 1){
    p.lower <- pnorm(u_k[K, 1], mean=mean[K], lower.tail=T)
    p.upper <- pnorm(u_k[K, 2], mean=mean[K], lower.tail=F)
  } else {
    lower.prev <- u_k[1:(K-1), 1]
    upper.prev <- u_k[1:(K-1), 2]

    lower.K <- u_k[K, 1]
    upper.K <- u_k[K, 2]

    p.lower <- pmvnorm(
      lower=c(lower.prev, -Inf),
      upper=c(upper.prev, lower.K),
      corr=corr,
      mean=mean[1:K],
      algorithm=algorithm
    )
    p.upper <- pmvnorm(
      lower=c(lower.prev, upper.K),
      upper=c(upper.prev, Inf),
      corr=corr,
      mean=mean[1:K],
      algorithm=algorithm
    )
  }
  return(c(p.lower, p.upper))
}

#' Get a stage-wise p-value
#'
#' @param obs Observed z-statistic
#' @param u_k Matrix of boundaries for *previous* stages.
#'            There should be two columns, one for lower bound, one for upper.
#'            This allows asymmetric boundaries.
#' @param corr Correlation matrix for test statistics
#' @param mean Mean shift vector for the test statistic
#' @param algorithm Algorithm argument for pmvnorm
#'
#' @examples
#' # OBF boundaries
#' u_k <- c(4.332634, 2.963132, 2.359044, 2.014090)
#' u_k1 <- cbind(-u_k, u_k)
#' u_k2 <- cbind(rep(-Inf, 4), u_k)
#'
#' # This gives exactly alpha = 0.05 by putting in the last
#' # boundary.
#' n_k <- c(10, 20, 30, 40)
#' corr.1 <- corr.mat(n_k)
#' corr.2 <- corr.mat(n_k, rho=0.5)
#' get.pvalue.sw(obs=u_k1[4,][1], u_k=u_k1[1:3,], corr=corr.1)
#' # 0.05000002
#' get.pvalue.sw(obs=u_k1[4,][1], u_k=u_k2[1:3,], corr=corr.1)
#' get.pvalue.sw(obs=2.5, u_k=u_k1[1:3,], corr=corr.1)
#' # 0.0247869
#' get.pvalue.sw(obs=1.5, u_k=u_k1[1:3,], corr=corr.1)
#' # 0.135164
#' get.pvalue.sw(obs=-1.5, u_k=u_k1[1:3,], corr=corr.1)
#' # 0.135164
#' get.pvalue.sw(obs=-1.5, u_k=u_k1[1:3,], corr=corr.2)
#' # 0.1461026
get.pvalue.sw <- function(obs, u_k, corr, mean=NULL,
                          algorithm=Miwa(steps=1000),
                          type="two-sided"){

  .check.type(type)

  K <- nrow(u_k) + 1

  p.upper.tot <- 0
  p.lower.tot <- 0

  for(i in 1:K){
    if(i == 1){

      if(nrow(u_k) == 0){ # reject at first stage
        p.lower <- pnorm(obs, lower.tail=T)
        p.upper <- pnorm(obs, lower.tail=F)
      } else {
        p.reject <- .reject.prob.k(u_k[i, ], mean=mean, corr=matrix(1))
        p.lower <- p.reject[1]
        p.upper <- p.reject[2]
      }

    } else {
      corr.i <- corr[1:i, 1:i]
      lower.prev <- u_k[1:(i-1), 1]
      upper.prev <- u_k[1:(i-1), 2]

      if(i < K){
        p.reject <- .reject.prob.k(u_k[1:i, ], mean=mean[1:i], corr=corr.i)
        p.lower <- p.reject[1]
        p.upper <- p.reject[2]
      } else {
        p.lower <- pmvnorm(
          lower=c(lower.prev, -Inf),
          upper=c(upper.prev, obs),
          corr=corr.i,
          mean=mean[1:i],
          algorithm=algorithm
        )
        p.upper <- pmvnorm(
          lower=c(lower.prev, obs),
          upper=c(upper.prev, Inf),
          corr=corr.i,
          mean=mean[1:i],
          algorithm=algorithm
        )
      }
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

get.pvalue.sm <- function(obs, u_k, corr, mean=NULL,
                          algorithm=Miwa(steps=1000),
                          type="two-sided"){

  .check.type(type)
  K <- nrow(u_k) + 1

  for(i in 1:K){

  }

  return(p)
}
