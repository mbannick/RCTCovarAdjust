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

.get.pval.from.type <- function(type, lower, upper){
  if(type == "two-sided"){
    return(2*min(upper, lower))
  } else if(type == "lower"){
    return(upper[1])
  } else if(type == "upper"){
    return(lower[1])
  }
}

.check.inputs <- function(u_k, ...){
  dim <- nrow(u_k) + 1
  args <- list(...)
  for(arg in args){
    if(!length(arg) %in% c(1, dim)) stop()
  }
}

.pmvnorm.trycatch <- function(bounds, ...){
  result <- pmvnorm(lower=bounds[, 1], upper=bounds[, 2],
                    ..., algorithm=Miwa(steps=1000))
  if(is.nan(result)){
    result <- pmvnorm(lower=bounds[, 1], upper=bounds[, 2],
                      ..., algorithm=GenzBretz(abseps=1e-8))
  }
  # Some of the values returned are negative, just from precision
  # issues. Need to take max of 0 and p-value.
  return(max(result, 0))
}

.pmvnorm.list <- function(blocks, ...){

  if(!is.list(blocks)) blocks <- list(blocks)
  val <- 0

  for(bound in blocks){
    result <- .pmvnorm.trycatch(bounds=bound, ...)
    val <- val + result
  }
  return(val)
}

#' Get the probability of rejecting at stage k.
#' We need this to return a two-element vector because
#' for the orderings, crossing lower *or* upper boundaries
#' is important.
#'
#' @param u_k Vector of boundaries up through stage k
#' @param corr Correlation matrix for test statistics
#'
#' @examples
#' u_k <- c(4.332634, 2.963132, 2.359044)
#' u_k1 <- cbind(-u_k, u_k)
#' u_k2 <- cbind(rep(-Inf, 3), u_k)
#' n_k <- c(10, 20, 30)
#' corr.1 <- corr.mat(n_k)
#' corr.2 <- corr.mat(n_k, rho=0.5, mis=c(F, F, T))
#' .reject.prob.k(u_k=u_k1, corr=corr.1)
#' .reject.prob.k(u_k=u_k2, corr=corr.1)
.reject.prob.k <- function(u_k, corr){

  if(!"matrix" %in% class(u_k)) u_k <- matrix(u_k, nrow=1)

  K <- nrow(u_k)

  if(K == 1){
    p.lower <- pnorm(u_k[K, 1], lower.tail=T)
    p.upper <- pnorm(u_k[K, 2], lower.tail=F)
  } else {
    prev <- u_k[1:(K-1), ]

    lower.K <- u_k[K, 1]
    upper.K <- u_k[K, 2]

    lower.block <- rbind(prev, c(-Inf, lower.K))
    upper.block <- rbind(prev, c(upper.K, Inf))

    p.lower <- .pmvnorm.list(blocks=lower.block, corr=corr)
    p.upper <- .pmvnorm.list(blocks=upper.block, corr=corr)
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
#' n_k1 <- c(10, 20, 30, 40)
#' n_k2 <- c(10, 20)
#'
#' get.pvalue.sw(obs=u_k1[4,][1], u_k=u_k1[1:3,], n_k=n_k1, rho=0.1,
#'               ancova_monitor=F, ancova_test=F, last_stage=T)
#' # 0.05000002
#' get.pvalue.sw(obs=u_k1[4,][1], u_k=u_k1[1:3,], n_k=n_k1, rho=0.1,
#'               ancova_monitor=F, ancova_test=T, last_stage=T)
#' # 0.06235481
#' get.pvalue.sw(obs=u_k1[4,][1], u_k=u_k1[1:3,], n_k=n_k1, rho=0.9,
#'               ancova_monitor=F, ancova_test=T, last_stage=T)
#' # 0.05304498
#' get.pvalue.sw(obs=1.5, u_k=u_k1[1:3,], n_k=n_k1, ancova_monitor=F,
#'               last_stage=T)
#' # 0.135164
#' get.pvalue.sw(obs=-1.5, u_k=u_k1[1:3,], n_k=n_k1, ancova_monitor=F,
#'               last_stage=T)
#' # 0.135164
#' get.pvalue.sw(obs=-1.5, u_k=u_k1[1:3,], n_k=n_k1, ancova_monitor=F,
#'               rho=0.8, last_stage=T)
#' # 0.1399254
#' get.pvalue.sw(obs=2.35, u_k=matrix(u_k1[1:1,], nrow=1),
#'               n_k=n_k2, ancova_monitor=F,
#'               rho=0.8, last_stage=T)
#' 0.01877924
#' get.pvalue.sw(obs=qnorm(0.975), u_k=NULL,
#'               n_k=c(10), ancova_monitor=F,
#'               rho=0.9, last_stage=T)
#' # 0.05
#' get.pvalue.sw(obs=qnorm(0.975), u_k=matrix(u_k1[1:1,], nrow=1),
#'               n_k=n_k2, ancova_monitor=T,
#'               rho=0.8, last_stage=F)
get.pvalue.sw <- function(obs, u_k, n_k, rho=1,
                          ancova_monitor=F,
                          ancova_test=T,
                          last_stage=F,
                          type="two-sided"){
  .check.type(type)

  # K is not necessarily the number of stages, it is
  # the number of hypothesis tests being conducted.
  K <- length(n_k)
  switch <- ancova_monitor != ancova_test

  if(!last_stage & switch){
    if(nrow(u_k) != length(n_k)) stop()
  }
  if(last_stage & !is.null(u_k)){
    if(nrow(u_k) != (length(n_k)-1)) stop()
  }

  p.upper.tot <- 0
  p.lower.tot <- 0

  for(i in 1:K){
    if(i == 1){

      if(is.null(u_k)){ # reject at first stage
        p.lower <- pnorm(obs, lower.tail=T)
        p.upper <- pnorm(obs, lower.tail=F)
      } else {
        p.reject <- .reject.prob.k(u_k[i, ], corr=matrix(1))
        p.lower <- p.reject[1]
        p.upper <- p.reject[2]
      }

    } else {
      if(i < K){

        corr.i <- corr.mat(n_k[1:i], rho=rho)
        p.reject <- .reject.prob.k(u_k[1:i, ], corr=corr.i)
        p.lower <- p.reject[1]
        p.upper <- p.reject[2]

      } else {

        prev <- u_k[1:(i-1), ]

        if(last_stage | !switch){

          if(!switch){
            mis <- rep(F, K)
          } else {
            # If last stage and switch
            mis <- c(rep(F, K-1), T)
          }

          corr.i <- corr.mat(n_k, rho=rho, mis=mis)

          lower.block <- rbind(prev, c(-Inf, obs))
          upper.block <- rbind(prev, c(obs, Inf))

        } else {
          corr.i <- corr.mat(n_k=c(n_k, n_k[K]), rho=rho,
                             mis=c(rep(F, i), switch))

          l.here <- u_k[K, 1]
          u.here <- u_k[K, 2]

          lower.reject <- c(-Inf, l.here)
          upper.reject <- c(u.here, Inf)

          lower.tail <- c(-Inf, obs)
          upper.tail <- c(obs, Inf)

          lower.block.1 <- rbind(prev, lower.reject, lower.tail)
          upper.block.1 <- rbind(prev, upper.reject, upper.tail)

          lower.block.2 <- rbind(prev, upper.reject, lower.tail)
          upper.block.2 <- rbind(prev, lower.reject, upper.tail)

          lower.block <- list(lower.block.1, lower.block.2)
          upper.block <- list(upper.block.1, upper.block.2)

        }

        p.lower <- .pmvnorm.list(blocks=lower.block, corr=corr.i)
        p.upper <- .pmvnorm.list(blocks=upper.block, corr=corr.i)

      }
    }
    p.upper.tot <- p.upper.tot + p.upper
    p.lower.tot <- p.lower.tot + p.lower
  }

  val <- .get.pval.from.type(type, lower=p.lower.tot, upper=p.upper.tot)
  return(val)
}

#' Get p-value using the sample mean ordering.
#' Needs to take in the obs as a standardized sample mean.
#'
#' @param obs Standardized sample mean observation
#' @param u_K Boundaries for all stages except the last stage K
#' @param n_K Sample sizes for all stages up through max stage K
#'            This is needed to compute the variance of the sample mean.
#' @param type Type of p-value (upper, lower, two-sided)
#'
#' @examples
#' # OBF boundaries
#' u_k <- c(4.332634, 2.963132, 2.359044, 2.014090)
#' u_k1 <- cbind(-u_k, u_k)
#' u_k2 <- cbind(rep(-Inf, 4), u_k)
#'
#' # This gives exactly alpha = 0.05 by putting in the last
#' # boundary.
#' n_K <- c(10, 20, 30, 40)
#'
#' get.pvalue.sm(obs=-2.01409/sqrt(40), u_K=u_k1[1:3,], n_K=n_K,
#'               ancova_monitor=F, ancova_test=F, rho=1.0)
#' get.pvalue.sm(obs=-2.01409/sqrt(40), u_K=u_k1[1:3,], n_K=n_K,
#'               ancova_monitor=F, ancova_test=T, rho=0)
#' # 0.05000002
#' get.pvalue.sm(obs=-2.01409/sqrt(40), u_K=u_k1[1:3,], n_K=n_K,
#'               ancova_monitor=F, ancova_test=T, rho=0.4)
#' get.pvalue.sm(obs=-2.01409/sqrt(40), u_K=u_k1[1:3,], n_K=n_K,
#'               ancova_monitor=F, ancova_test=T, rho=0.8)
#' for(i in seq(0.01, 0.9, 0.1)){
#'     print(get.pvalue.sm(obs=-2.01409/sqrt(40), u_K=u_k1[1:3,], n_K=n_K,
#'               ancova_monitor=F, ancova_test=T, rho=i))
#' }
#'
#'
#' get.pvalue.sm(obs=-5/sqrt(40), u_K=u_k1[1:3,], n_K=n_K,
#'               ancova_monitor=F, ancova_test=F, rho=1.0)
get.pvalue.sm <- function(obs, u_K, n_K, rho=1,
                          ancova_monitor=F,
                          ancova_test=F,
                          type="two-sided"){
  .check.type(type)
  .check.inputs(u_K, n_K)

  # Get max stages
  K <- length(n_K)

  # Convert sample mean to likelihood ratio
  lr_obs <- obs * sqrt(n_K[K])

  # Indicator for whether we're switching between ANOVA and ANCOVA
  # at the testing stage.
  switch <- ancova_monitor != ancova_test

  p.upper.tot <- 0
  p.lower.tot <- 0

  for(i in 1:K){

    if(i == 1){
      prev <- c()
    } else {
      prev <- u_K[1:(i-1), ]
    }

    if(i < K){

      if(switch){

        corr.i <- corr.mat(n_k=c(n_K[1:i], n_K[i]), rho=rho,
                           mis=c(rep(F, i), T),
                           sme=c(rep(F, i), T))

        # Get the critical value for rejection
        # at this stage.
        l.here <- u_K[i, 1]
        u.here <- u_K[i, 2]

        lower.reject <- c(-Inf, l.here)
        upper.reject <- c(u.here, Inf)

        lower.tail <- c(-Inf, obs)
        upper.tail <- c(obs, Inf)

        lower.block.1 <- rbind(prev, lower.reject, lower.tail)
        upper.block.1 <- rbind(prev, upper.reject, upper.tail)

        # When test rejection and observed don't match up
        # in terms of direction -- should happen rarely.
        # These probabilities should be really small.
        lower.block.2 <- rbind(prev, upper.reject, lower.tail)
        upper.block.2 <- rbind(prev, lower.reject, upper.tail)

        lower.block <- list(lower.block.1, lower.block.2)
        upper.block <- list(upper.block.1, upper.block.2)

      } else {

        # When performing a test and wanting a SM p-value, (before the last stage)
        # need to only get the correlation matrix for the Z-statistics
        # rather than the sample mean because otherwise, it's degenerate.

        corr.i <- corr.mat(n_k=c(n_K[1:i]),
                           mis=rep(F, i), sme=F)

        l.here <- min(u_K[i, 1], lr_obs)
        u.here <- max(u_K[i, 2], lr_obs)

        lower.tail <- c(-Inf, l.here)
        upper.tail <- c(u.here, Inf)

        lower.block <- rbind(prev, lower.tail)
        upper.block <- rbind(prev, upper.tail)

      }

    } else {
      corr.i <- corr.mat(n_k=n_K, rho=rho,
                         mis=c(rep(F, K-1), switch),
                         sme=c(rep(F, K-1), T))

      lower.block <- rbind(prev, c(-Inf, obs))
      upper.block <- rbind(prev, c(obs, Inf))
    }

    p.lower <- .pmvnorm.list(blocks=lower.block, sigma=corr.i)
    p.upper <- .pmvnorm.list(blocks=upper.block, sigma=corr.i)

    p.upper.tot <- p.upper.tot + p.upper
    p.lower.tot <- p.lower.tot + p.lower
  }

  val <- .get.pval.from.type(type, lower=p.lower.tot, upper=p.upper.tot)
  return(val)
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
