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

.check.inputs <- function(u_k, ...){
  dim <- nrow(u_k) + 1
  args <- list(...)
  for(arg in args){
    if(!length(arg) %in% c(1, dim)) stop()
  }
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
    lower.prev <- u_k[1:(K-1), 1]
    upper.prev <- u_k[1:(K-1), 2]

    lower.K <- u_k[K, 1]
    upper.K <- u_k[K, 2]

    p.lower <- pmvnorm(
      lower=c(lower.prev, -Inf),
      upper=c(upper.prev, lower.K),
      corr=corr,
      algorithm=Miwa(steps=1000)
    )
    p.upper <- pmvnorm(
      lower=c(lower.prev, upper.K),
      upper=c(upper.prev, Inf),
      corr=corr,
      algorithm=Miwa(steps=1000)
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
#' n_k1 <- c(10, 20, 30, 40)
#' n_k2 <- c(10, 20)
#'
#' ancova1 <- c(F, F, F, F)
#' ancova2 <- c(F, T)
#'
#' get.pvalue.sw(obs=u_k1[4,][1], u_k=u_k1[1:3,], n_k=n_k1,
#'               ancova_monitor=ancova1, ancova_test=F, last_stage=T)
#' # 0.05000002
#' get.pvalue.sw(obs=2.5, u_k=u_k1[1:3,], n_k=n_k1, ancova_monitor=ancova1,
#'               last_stage=T)
#' # 0.0247869
#' get.pvalue.sw(obs=1.5, u_k=u_k1[1:3,], n_k=n_k1, ancova_monitor=ancova1,
#'               last_stage=T)
#' # 0.135164
#' get.pvalue.sw(obs=-1.5, u_k=u_k1[1:3,], n_k=n_k1, ancova_monitor=ancova1,
#'               last_stage=T)
#' # 0.135164
#' get.pvalue.sw(obs=-1.5, u_k=u_k1[1:4,], n_k=n_k1, ancova_monitor=ancova1,
#'               rho=0.8, last_stage=T)
#' # 0.1399254
#' get.pvalue.sw(obs=2.963132, u_k=u_k1[1:2,], n_k=n_k2, ancova_monitor=ancova2,
#'               rho=0.8, last_stage=F)
#' # 0.005240675
get.pvalue.sw <- function(obs, u_k, n_k, rho=1,
                          ancova_monitor=F,
                          ancova_test=T,
                          last_stage=F,
                          type="two-sided"){
  .check.type(type)

  # K is not necessarily the number of stages, it is
  # the number of hypothesis tests being conducted.
  K <- length(n_k)
  if(length(ancova_monitor) == 1) ancova_monitor <- rep(ancova_monitor, K)

  # If in the last stage, replace whatever was going to happen with monitoring
  # because we ignore the last test in the ordering.
  if(last_stage){
    ancova_monitor[K] <- ancova_test
  }
  if(!last_stage){
    if(nrow(u_k) != length(n_k)) stop()
  }

  p.upper.tot <- 0
  p.lower.tot <- 0

  for(i in 1:K){
    if(i == 1){

      if(nrow(u_k) == 0){ # reject at first stage
        p.lower <- pnorm(obs, lower.tail=T)
        p.upper <- pnorm(obs, lower.tail=F)
      } else {
        p.reject <- .reject.prob.k(u_k[i, ], corr=matrix(1))
        p.lower <- p.reject[1]
        p.upper <- p.reject[2]
      }

    } else {

      lower.prev <- u_k[1:(i-1), 1]
      upper.prev <- u_k[1:(i-1), 2]

      if(i < K){
        corr.i <- corr.mat(n_k[1:i],
                           rho=rho,
                           mis=ancova_monitor[1:i])
        p.reject <- .reject.prob.k(u_k[1:i, ], corr=corr.i)
        p.lower <- p.reject[1]
        p.upper <- p.reject[2]
      } else {
        if(last_stage){

          corr.i <- corr.mat(n_k,
                             rho=rho,
                             mis=ancova_monitor)

          p.lower <- pmvnorm(
            lower=c(lower.prev, -Inf),
            upper=c(upper.prev, obs),
            corr=corr.i,
            algorithm=Miwa(steps=1000)
          )
          p.upper <- pmvnorm(
            lower=c(lower.prev, obs),
            upper=c(upper.prev, Inf),
            corr=corr.i,
            algorithm=Miwa(steps=1000)
          )
        } else {
          browser()

          # THIS IS WHAT I HAD OVERLOOKED BEFORE
          corr.i <- corr.mat(c(n_k, n_k[K]),
                             rho=rho,
                             mis=c(ancova_monitor, ancova_test))

          l.here <- u_k[K, 1]
          u.here <- u_k[K, 2]

          p.lower.1 <- pmvnorm(
            lower=c(lower.prev, -Inf, -Inf),
            upper=c(upper.prev, l.here, obs),
            corr=corr.i,
            algorithm=Miwa(steps=1000)
          )
          p.lower.2 <- pmvnorm(
            lower=c(lower.prev, u.here, -Inf),
            upper=c(upper.prev, Inf, obs),
            corr=corr.i,
            algorithm=Miwa(steps=1000)
          )
          p.upper.1 <- pmvnorm(
            lower=c(lower.prev, -Inf, obs),
            upper=c(upper.prev, l.here, Inf),
            corr=corr.i,
            algorithm=Miwa(steps=1000)
          )
          p.upper.2 <- pmvnorm(
            lower=c(lower.prev, u.here, obs),
            upper=c(upper.prev, Inf, Inf),
            corr=corr.i,
            algorithm=Miwa(steps=1000)
          )
          p.lower <- p.lower.1 + p.lower.2
          p.upper <- p.upper.1 + p.upper.2
        }
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
#' n_k1 <- c(10, 20, 30, 40)
#' n_k2 <- c(10, 20, 30, 40, 40)
#' n_k3 <- c(10, 20, 20)
#'
#' mis1 <- c(F, F, F, F)
#' mis2 <- c(F, F, F, F, T)
#' mis3 <- c(F, F, T)
#'
#' get.pvalue.sm(obs=u_k1[4,][1], u_k=u_k1[1:3,], n_k=n_k1, mis=mis1)
#' # 0.05000002
#' get.pvalue.sm(obs=2.5, u_k=u_k1[1:3,], n_k=n_k1, mis=mis1)
#' # 0.0247869
#' get.pvalue.sm(obs=1.5, u_k=u_k1[1:3,], n_k=n_k1, mis=mis1)
#' # 0.135164
#' get.pvalue.sm(obs=-1.5, u_k=u_k1[1:3,], n_k=n_k1, mis=mis1)
#' # 0.135164
#' get.pvalue.sm(obs=-1.5, u_k=u_k1[1:4,], n_k=n_k2, rho=0.8, mis=mis2)
#' # 0.1491
#' get.pvalue.sm(obs=2.963132, u_k=u_k1[1:2,], n_k=n_k3, rho=0.8, mis=mis3)
#' # 0.005240675
get.pvalue.sm <- function(obs, u_K, n_K, rho=1, mis=F,
                          type="two-sided"){
  .check.type(type)
  .check.inputs(u_k, n_k, mis)

  # K is not necessarily the number of stages, it is
  # the number of hypothesis tests being conducted.
  K <- length(n_k)
  if(length(mis) == 1) mis <- rep(mis, K)
  mis.end <- all(mis) | all(!mis)

  p.upper.tot <- 0
  p.lower.tot <- 0

  for(i in 1:(K-1)){

    corr.i <- corr.mat(n_k=c(n_k[1:i], n_k[i]),
                       rho=rho,
                       mis=c(mis[1:i], mis.end),
                       sme=c(rep(F, i), T))

    lower.prev <- u_k[1:(i-1), 1]
    upper.prev <- u_k[1:(i-1), 2]

    lower.this <- u_k[i, 2]
    upper.this <-

    p.lower <- pmvnorm(
      lower=c(lower.prev, -Inf),
      upper=c(upper.prev, obs),
      Sigma=corr.i,
      algorithm=Miwa(steps=1000)
    )
    p.upper <- pmvnorm(
      lower=c(lower.prev, obs),
      upper=c(upper.prev, Inf),
      Sigma=corr.i,
      algorithm=Miwa(steps=1000)
    )

    p.upper.tot <- p.upper.tot + p.upper
    p.lower.tot <- p.lower.tot + p.lower

  }

  corr.i <- corr.mat(n_k=c(n_k[1:i], n_k[i]),
                     rho=rho,
                     mis=c(mis[1:i], mis.end),
                     sme=c(rep(F, i), T))

  lower.prev <- u_k[1:(K-1), 1]
  upper.prev <- u_k[1:(K-1), 2]

  p.lower <- pmvnorm(
    lower=c(lower.prev, -Inf),
    upper=c(upper.prev, obs),
    Sigma=corr.i,
    algorithm=Miwa(steps=1000)
  )
  p.upper <- pmvnorm(
    lower=c(lower.prev, obs),
    upper=c(upper.prev, Inf),
    Sigma=corr.i,
    algorithm=Miwa(steps=1000)
  )

  p.upper.tot <- p.upper.tot + p.upper
  p.lower.tot <- p.lower.tot + p.lower

  if(type == "two-sided"){
    return(2*min(p.upper.tot, p.lower.tot))
  } else if(type == "lower"){
    return(p.upper.tot[1])
  } else if(type == "upper"){
    return(p.lower.tot[1])
  }
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
