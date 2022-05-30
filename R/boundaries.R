library(magrittr)
library(MASS)
library(mvtnorm)

source("~/repos/RCTCovarAdjust/R/constants.R")
source("~/repos/RCTCovarAdjust/R/covariance.R")

#' Get Pocock or OBF-type boundaries by correlation matrix.
#' Compute boundaries based on a correlation matrix between
#' the test statistics.
get.bound.by.corr <- function(corr, obf=FALSE,
                              power=0.05, algorithm=Miwa(steps=1000)){
  K <- nrow(corr)
  get.power <- function(c){

    u_k <- rep(c, K)
    if(obf) u_k <- u_k / sqrt(1:K)

    val <- 1 - pmvnorm(lower=-u_k,
                       upper=u_k,
                       mean=rep(0, K), corr=corr,
                       algorithm=algorithm)
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
#' @param u_k Previous boundaries, in a matrix.
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

#' Get Pocock or OBF-style boundaries at the design stage
#'
#' @export
#' @param obf
#' @param rates A vector of information rates (between 0 and 1)
#' @param u_k An optional matrix of previous boundaries
#' @param rho Fraction of variance explained by fitting ANCOVA.
#' @param change A vector indicating which stages use ANOVA v. ANCOVA.
#' @examples
#' t <- 1:4/4
#'
#' # approximate Pocock boundaries
#' get.boundaries.design(rates=t, obf=FALSE)
#' get.boundaries.design(rates=t, obf=FALSE, rho=0.1, change=c(0, 1, 0, 0))
#' get.boundaries.design(rates=t, obf=TRUE, rho=0.1, change=c(1, 0, 0, 0))
get.boundaries.design <- function(rates, obf,
                                  rho=1, change=0,
                                  algorithm=Miwa(steps=1000)){

  if(length(change) == 1) change <- rep(change, length(rates))
  if(length(change) != length(rates)) stop("Change vector needs to be
                                           the same length as the number of stages.")

  corr <- corr.mat(rates, rho=rho, mis=as.logical(change))
  bounds <- get.bound.by.corr(corr, obf=obf)
  return(bounds)
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
#' @param rates A vector of information rates (between 0 and 1)
#' @param u_k An optional matrix of previous boundaries
#' @param rho Fraction of variance explained by fitting ANCOVA.
#' @param change A vector indicating which stages use ANOVA v. ANCOVA.
#' @examples
#' # information rates
#' t <- 1:4/4
#'
#' # approximate Pocock boundaries
#' get.boundaries.aspend(a.func=pocock.spend(0.05), rates=t)
#' get.boundaries.aspend(a.func=pocock.spend(0.05), rates=t,
#'                       rho=0.5, change=c(0, 0, 0, 1))
get.boundaries.aspend <- function(a.func, rates,
                                  u_k=NULL, rho=1, change=0,
                                  algorithm=Miwa(steps=1000)){

  # Get sample size increments
  K <- length(rates)

  # Number of *fixed* previous bounds
  K_prev <- length(u_k)

  # Get cumulative alpha
  a.cuml <- a.func(rates)

  if(length(change) == 1) change <- rep(change, K)
  if(length(change) != K) stop("Change vector needs to be the same length
                               as the number of stages.")

  bounds <- c()
  for(i in 1:K){
    if(i > K_prev){
      # Create covariance matrix
      corr <- corr.mat(n_k=rates[1:i], rho=rho, mis=as.logical(change[1:i]))
      bound <- solve.boundary(power=a.cuml[i], corr=corr,
                              u_k=bounds, algorithm=algorithm)
    } else {
      bound <- u_k[i]
    }
    bounds <- rbind(bounds, c(-bound, bound))
  }
  return(bounds)
}
