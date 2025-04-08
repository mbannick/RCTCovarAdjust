source("~/repos/RCTCovarAdjust/R/constants.R")
source("~/repos/RCTCovarAdjust/R/covariance.R")

# #' Get Pocock or OBF-type boundaries by correlation matrix.
# #' Compute boundaries based on a correlation matrix between
# #' the test statistics.
# #'
# #' @import MASS
# #' @import mvtnorm
get.bound.by.corr <- function(corr, obf=FALSE, unequal_type=FALSE,
                              t_k=NULL, u_k=NULL,
                              power=0.05, algorithm=Miwa(steps=1000)){
  K <- nrow(corr)
  if(length(t_k) != K) stop("Information fractions vector must
                            be same size as correlation matrix.")

  bound.func <- function(c){
    if(obf){
      if(!unequal_type){
        u <- c / sqrt(1:K)
      } else {
        if(is.null(t_k)){
          stop("Must provide information fractions if desire OBF bounds
             with unequal sample size constant correction (eq. 3.7 in Wassmer and Brannath).")
        } else {
          u <- c / sqrt(t_k/t_k[1])
        }
      }
    } else {
      u <- rep(c, K)
    }

    if(!is.null(u_k)){
      u[1:length(u_k)] <- u_k
    }
    return(u)
  }
  get.power <- function(c){
    us <- bound.func(c)

    val <- 1 - pmvnorm(lower=-us,
                       upper=us,
                       mean=rep(0, K), corr=corr,
                       algorithm=algorithm)
    return(val)
  }
  f <- function(x) get.power(x) - power
  s <- uniroot(f, interval=c(0, 100))

  return(bound.func(s$root))
}

# #' Root solve for a particular alpha level, two-sided.
# #'
# #' @param power The cumulative type I error (or power generally) to solve for.
# #' @param corr Correlation matrix for the test statistics.
# #' @param mean Mean vector for the multivariate normal. Defaults to 0.
# #' @param u_k Previous boundaries, in a matrix.
# #' @param tol Tolerance parameter for binary search.
# #' @param ... Additional arguments for pmvnorm algorithm.
# #'
# #' @import MASS
# #' @import mvtnorm
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
#' @param obf Whether to use OBF (TRUE) or Pocock (FALSE)
#' @param rates A vector of information rates (between 0 and 1)
#' @param unequal_type Correction for unequal sample sizes across stage
#' @param rho Fraction of variance explained by fitting ANCOVA.
#' @param change A vector indicating which stages use ANOVA v. ANCOVA.
#'
#' @import MASS
#' @import mvtnorm
#' @export
#' @examples
#' # Information fractions
#' t <- 1:3/3
#'
#' # OBF-type boundaries
#' get.boundaries.design(rates=t, obf=TRUE)
#' get.boundaries.design(rates=c(0.3, 0.9, 1.0), obf=TRUE)
#' get.boundaries.design(rates=c(0.3, 0.9, 1.0), obf=TRUE, unequal_type=TRUE)
#'
#' # ANCOVA at last stage, R^2 = 0.5
#' get.boundaries.design(rates=c(0.3, 0.9, 1.0), obf=TRUE, unequal_type=TRUE, rho=sqrt(0.5), change=c(0, 0, 1))
#'
#' # ANCOVA at last two stages, R^2 = 0.5
#' get.boundaries.design(rates=t, obf=TRUE, rho=sqrt(0.5), change=c(1, 0, 0))
get.boundaries.design <- function(rates, obf, unequal_type=FALSE,
                                  rho=1, change=0, u_k=NULL,
                                  algorithm=Miwa(steps=1000)){

  if(length(change) == 1) change <- rep(change, length(rates))
  if(length(change) != length(rates)) stop("Change vector needs to be
                                           the same length as the number of stages.")

  corr <- corr.mat(rates, rho=rho, mis=as.logical(change))
  bounds <- get.bound.by.corr(corr, t_k=rates, obf=obf, unequal_type=unequal_type, u_k=u_k)
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
#' @import MASS
#' @import mvtnorm
#' @examples
#' # Information rates
#' t <- 1:4/4
#'
#' # Approximate Pocock boundaries, with R^2 = 0.5, ANCOVA at last stage
#' get.boundaries.aspend(a.func=pocock.spend(0.05), rates=t)
#' get.boundaries.aspend(a.func=pocock.spend(0.05), rates=t,
#'                       rho=sqrt(0.5), change=c(0, 0, 0, 1))
get.boundaries.aspend <- function(a.func, rates,
                                  u_k=NULL, rho=1, change=0,
                                  algorithm=Miwa(steps=1000)){

  browser()

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
