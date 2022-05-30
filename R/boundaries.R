library(magrittr)
library(MASS)
library(mvtnorm)

source("~/repos/RCTCovarAdjust/R/constants.R")
source("~/repos/RCTCovarAdjust/R/covariance.R")

.get.power <- function(c, rates, obf=FALSE, rho=1){

  K <- length(rates)
  u_k <- rep(c, K)

  if(obf) u_k <- u_k / sqrt(1:K)

  Sigma <- corr.mat(rates, rho=rho, mis=c(rep(F, K-1), T))
  val <- 1 - pmvnorm(lower=-u_k,
                     upper=u_k,
                     mean=rep(0, K), corr=Sigma,
                     algorithm=Miwa(steps=1000))

  return(val)
}

#' Get boundaries
#'
#' @param K number of stages
#' @param obf O'Brien-Fleming bounds (TRUE) or Pocock (FALSE)
#' @param rho Reduction in variance due to ANCOVA
#' @param power Alpha-level
#'
#' @examples
#' get.bound(3, obf=TRUE)
get.bound <- function(rates, obf=FALSE, rho=1, power=0.05){

  K <- length(rates)
  f <- function(x) .get.power(x, rates=rates, obf=obf, rho=rho) - power
  s <- uniroot(f, interval=c(0, 100))

  if(obf){
    return(s$root / sqrt(1:K))
  } else {
    return(rep(s$root, K))
  }
}

get.bound.by.corr <- function(corr, obf=FALSE, power=0.05){
  K <- nrow(corr)
  get.power <- function(c){

    u_k <- rep(c, K)
    if(obf) u_k <- u_k / sqrt(1:K)

    val <- 1 - pmvnorm(lower=-u_k,
                       upper=u_k,
                       mean=rep(0, K), corr=corr,
                       algorithm=Miwa(steps=1000))
  }
  f <- function(x) get.power(x) - power
  s <- uniroot(f, interval=c(0, 100))

  if(obf){
    return(s$root / sqrt(1:K))
  } else {
    return(rep(s$root, K))
  }

}

#' Root solve for a particular alpha level, two-sided.
#'
#' @param power The cumulative type I error (or power generally) to solve for.
#' @param corr Correlation matrix for the test statistics.
#' @param mean Mean vector for the multivariate normal. Defaults to 0.
#' @param u_k Previous boundaries.
#' @param tol Tolerance parameter for binary search.
#' @param ... Additional arguments for pmvnorm algorithm.
solve.boundary <- function(power, corr=NULL, u_k=NULL,
                           tol=.Machine$double.eps, algorithm=Miwa(steps=1000)){

  if(is.null(u_k)){
    bin <- function(x) 2 * pnorm(x, mean=0, sd=1, lower.tail=F) - power
  } else {
    bin <- function(x) 1 - power - pmvnorm(lower=c(u_k[, 1], -x),
                                           upper=c(u_k[, 2], x),
                                           mean=rep(0, nrow(u_k) + 1),
                                           corr=corr,
                                           algorithm=algorithm)
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
#'   where a.func(0) = 0 and a.func(1) = a
#'   where a is the type I error desired
#' @param a The type I error desired
#' @param rates A vector of information rates (between 0 and 1)
#' @param N maximum total sample size
#' @param n_sims Number of Monte-Carlo simulations
#' @param rho Fraction of variance explained by fitting ANCOVA.
#'
#' @examples
#' # information rates
#' t <- 1:4/4
#'
#' # approximate Pocock boundaries
#' get.boundaries(a.func=POCOCK.SPEND(0.05), rates=t)
#' get.boundaries(a.func=POCOCK.SPEND(0.05), rates=t, rho=0.5)
#'
#' # approximate O'Brien-Fleming boundaries
#' get.boundaries(a.func=OBF.SPEND(0.05), rates=t)
#' get.boundaries(a.func=OBF.SPEND(0.05), rates=t,
#'                u_k=c(4.332634, 2.963132, 2.359044))
#' get.boundaries(a.func=OBF.SPEND(0.05), rates=t, rho=0.5)
get.boundaries <- function(a.func, rates,
                           u_k=c(), rho=1, algorithm=Miwa(steps=1000),
                           extra=FALSE){

  # Get sample size increments
  K <- length(rates)

  # Number of *fixed* previous bounds
  K_prev <- length(u_k)

  # Append 0 onto the rates
  a.cuml <- a.func(rates)

  bounds <- c()
  for(i in 1:K){
    if(i > K_prev){

      # Create covariance matrix
      Sigma <- corr.mat(rates[1:i])
      if(i == K){
        Sigma <- corr.mat(rates[1:i], rho=rho)
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
