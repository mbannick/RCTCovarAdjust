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
#' @param ... Additional arguments for pmvnorm
#'
#' @examples
#' u_k <- c(4.332634, 2.963132, 2.359044)
#' u_k1 <- cbind(-u_k, u_k)
#' u_k2 <- cbind(rep(-Inf, 3), u_k)
#' n_k <- c(10, 20, 30)
#' corr.1 <- corr.mat(n_k)
#' corr.2 <- corr.mat(n_k, rho=0.5)
#' .reject.prob.k(u_k=u_k1, corr=corr.1)
#' .reject.prob.k(u_k=u_k2, corr=corr.1)
.reject.prob.k <- function(u_k, corr, ...){

  if(!"matrix" %in% class(u_k)) u_k <- matrix(u_k, nrow=1)

  K <- nrow(u_k)

  if(K == 1){
    p.lower <- pnorm(u_k[K, 1], lower.tail=T)
    p.upper <- pnorm(u_k[K, 2], lower.tail=F)
  } else {
    lower.prev <- u_k[1:(K-1), 1]
    upper.prev <- u_k[1:(K-1), 2]

    lower.K <- u_k[K, 1]
    upper.K <- u_k[K, 2]

    p.lower <- pmvnorm(
      lower=c(lower.prev, -Inf),
      upper=c(upper.prev, lower.K),
      corr=corr,
      ...
    )
    p.upper <- pmvnorm(
      lower=c(lower.prev, upper.K),
      upper=c(upper.prev, Inf),
      corr=corr,
      ...
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
#' @param ... Additional arguments for pmvnorm
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
get.pvalue.sw <- function(obs, u_k, corr,
                          type="two-sided", ...){

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
        p.reject <- .reject.prob.k(u_k[i, ], corr=matrix(1), ...)
        p.lower <- p.reject[1]
        p.upper <- p.reject[2]
      }

    } else {
      corr.i <- corr[1:i, 1:i]
      lower.prev <- u_k[1:(i-1), 1]
      upper.prev <- u_k[1:(i-1), 2]

      if(i < K){
        p.reject <- .reject.prob.k(u_k[1:i, ], corr=corr.i)
        p.lower <- p.reject[1]
        p.upper <- p.reject[2]
      } else {
        p.lower <- pmvnorm(
          lower=c(lower.prev, -Inf),
          upper=c(upper.prev, obs),
          corr=corr.i,
          ...
        )
        p.upper <- pmvnorm(
          lower=c(lower.prev, obs),
          upper=c(upper.prev, Inf),
          corr=corr.i,
          ...
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

#' Get p-value using the sample mean ordering.
#' Needs to take in the obs as a standardized sample mean.
#'
#' @param obs Standardized sample mean observation
#' @param u_K Boundaries for all stages up through max stage K
#' @param n_K Sample sizes for all stages up through max stage K
#'            This is needed to compute the variance of the sample mean.
#' @param corr Correlation matrix for all stages up through max stage K
#' @param type Type of p-value (upper, lower, two-sided)
#' @param ... Additional arguments for pmvnorm
get.pvalue.sm <- function(obs, u_K, corr, n_K,
                          type="two-sided", ...){

  .check.type(type)
  K <- nrow(n_K)

  p.lower.tot <- 0
  p.upper.tot <- 0

  p.reject.tot <- 0

  for(i in 1:K){

    # Create marginal density function for the sample mean
    # at stage i.
    sm.sd <- 1/n_K[i]**0.5
    sm.F <- function(x) pnorm(x, sd=sm.sd)
    sm.f <- function(x) pnorm(x, sd=sm.sd)

    # Calculate the constant that will normalize the
    # truncated CDF.
    constant <- 1 - sm.F(u_K[i, 2]) + sm.F(u_K[i, 1])

    # Create the truncated CDF.
    # This is the CDF for the sample mean conditional
    # on stopping at stage i.
    sm.F.trunc <- function(x){
      # If the val is in the left tail, just the CDF
      # up until the val
      if(x < u_K[i, 1]){
        val <- sm.F(x)
      # if the val is in the truncated area,
      # just the CDF up until the lower bound
      } else if(x >= u_K[i, 1] & x < u_K[i, 2]){
        val <- sm.F(u_K[i, 1])
      } else {
      # if the val is in the right tail, the CDF
      # for the lower bound and for the upper tail until the val.
        val <- sm.F(u_K[i, 1]) + (1 - sm.F(x) - sm.F(u_K[i, 2]))
      }
      return(val / constant)
    }

    if(i < K){
      # Get correlation matrix up until stage i
      corr.i <- corr[1:i, 1:i]

      # Get probability of rejection at stage i
      p.reject <- .reject.prob.k(u_K[1:i, ], corr=corr.i, ...)

      # The tail probabilities are the tail probabilities for sm.f.trunc
      # times the probability of rejecting in stage i.
      p.lower <- p.reject * sm.F.trunc(obs)
      p.upper <- p.reject * (1 - sm.F.trunc(obs))

      p.reject.tot <- p.reject.tot + p.reject

    } else {

      # (1 - p.reject.tot) is the probability
      # of making it to the last stage K.
      p.lower <- (1 - p.reject.tot) * sm.F(obs)
      p.upper <- (1 - p.reject.tot) * (1 - sm.F(obs))

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

  return(p)
}

#' This is not the function to expose to the user!
#'
#' We will need to calculate the boundaries for u_k
#' from an estimated standard deviation, and depends on
#' the procedure.
#'
#' @param obs Observed statistic
#' @param u_K Critical boundaries through max stage K
#' @param k The stage that it was stopped at
#' @param corr Correlation matrix among the test statistics
#' @param n_K Sample size through max stage K
get.pvalue <- function(obs, u_K, k, corr, n_K=c(),
                       ordering="stage-wise"){

  if(ordering == "sample mean" & length(n_K) == 0){
    stop("If using the stage-wise ordering, need to give the sample size
          for all data looks, including those planned at future stages.")
  }

  if(ordering == "stage-wise"){
    p <- get.pvalue.sw(obs=obs, u_k=u_K[1:k], corr=corr[1:k, 1:k])
  } else if(ordering == "sample mean"){
    p <- get.pvalue.sm(obs=obs, u_K=u_K, corr=corr, n_K=n_K)
  }

  return(p)
}
