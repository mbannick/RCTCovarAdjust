library(magrittr)
library(MASS)
library(mvtnorm)
source("R/trial-data.R")
source("R/trial-funcs.R")

#' Root solve for a particular alpha level, two-sided.
#'
#' @param power The cumulative type I error (or power generally) to solve for.
#' @param corr Correlation matrix for the test statistics.
#' @param mean Mean vector for the multivariate normal. Defaults to 0.
#' @param u_k Previous boundaries.
#' @param tol Tolerance parameter for binary search.
#' @param ... Additional arguments for pmvnorm algorithm.
solve.boundary <- function(power, mean=NULL, corr=NULL, u_k=NULL,
                           tol=.Machine$double.eps, ...){

  if(is.null(mean)) mean <- rep(0, length(u_k) + 1)

  if(is.null(u_k)){
    bin <- function(x) 2 * pnorm(x, mean=mean, sd=1, lower.tail=F) - power
  } else {
    bin <- function(x) 1 - power - pmvnorm(lower=c(-u_k, -x),
                                           upper=c(u_k, x),
                                           mean=mean, corr=corr, ...)
  }
  u <- uniroot(bin, c(0, 100), tol=tol)
  return(u$root)
}

#' Derive boundaries with alpha-spending.
#'
#' Get the boundaries required for a particular alpha-spending function
#' and an observed information rate out of a total sample size.
#'
#' Your data generator in this case needs to be the data generation
#' under the null hypothesis in order for this function to work correctly
#' with alpha-spending.
#'
#' @export
#' @param a.func A continuous, monotonic increasing function of t
#'   where a.func(a, t=0) = 0 and a.func(a, t=1) = a
#'   where a is the type I error desired
#' @param a The type I error desired
#' @param rates A vector of information rates (between 0 and 1)
#' @param N maximum total sample size
#' @param n_sims Number of Monte-Carlo simulations
#' @param rho Fraction of variance explained by fitting ANCOVA.
#'
#' @examples
#' set.seed(101)
#' # information rates
#' t <- 1:4/4
#'
#' # approximate Pocock boundaries
#' a.func.pocock <- function(a, t) a * log(1 + (exp(1) - 1) * t)
#' get.boundaries(a.func=a.func.pocock, a=0.05,
#'                rates=t, N=1000)
#' get.boundaries(a.func=a.func.pocock, a=0.05,
#'                rates=t, N=1000, algorithm=Miwa(steps=1000))
#'
#' # approximate O'Brien-Fleming boundaries
#' a.func.obf <- function(a, t) 4 * (1 - pnorm(qnorm(1-a/4)/sqrt(t)))
#' get.boundaries(a.func=a.func.obf, a=0.05,
#'                rates=t, N=1000)
#' get.boundaries(a.func=a.func.obf, a=0.05,
#'                rates=t, N=1000, algorithm=Miwa(steps=1000))
get.boundaries <- function(a.func, a, rates, N,
                           u_k=c(), rho=1, algorithm=Miwa(steps=1000)){

  # Get sample size increments
  K <- length(rates)
  n_k <- round(rates*N)
  n_k <- c(n_k[1], diff(n_k))

  # Number of *fixed* previous bounds
  K_prev <- length(u_k)

  # Create alpha spending function
  a.spend <- function(t) a.func(a=a, t)

  # Append 0 onto the rates
  a.cuml <- a.spend(rates)

  bounds <- c()
  for(i in 1:K){
    if(i > K_prev){

      # Create covariance matrix
      Sigma <- basic.cov(n_k[1:i])
      if(i == K){
        Sigma[K, 1:(K-1)] <- Sigma[K, 1:(K-1)] * rho
        Sigma[1:(K-1), K] <- Sigma[1:(K-1), K] * rho
      }
      bound <- solve.boundary(power=a.cuml[i], corr=Sigma,
                              u_k=bounds, algorithm=algorithm)
    } else {
      bound <- u_k[i]
    }
    bounds <- c(bounds, bound)
  }
  return(bounds)
}
