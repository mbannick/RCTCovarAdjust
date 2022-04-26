source("~/repos/RCTCovarAdjust/R/covariance.R")
source("~/repos/RCTCovarAdjust/R/boundaries.R")
source("~/repos/RCTCovarAdjust/R/constants.R")

#' Generating function to get boundaries
#' based on previous bounds and a correlation matrix
#' between the test statistics.
#'
#' Takes in an alpha spending function, information
#' fractions, and total sample size.
#'
#' @examples
#' func <- obf.spend(0.05)
#' rates <- 1:4/4
#'
#' bound.func <- get.boundary.closure(func, rates)
#' prev_bounds <- c(4.332634, 2.963132, 2.359044)
#' prev_bounds <- cbind(-prev_bounds, prev_bounds)
#' corr <- corr.mat(rates, rho=0.8, mis=c(F, F, F, T))
#' bound.func(prev_bounds, corr)
get.boundary.closure <- function(a.func, rates, est.bounds, a.type){

  if(est.bounds){
    get.boundary <- function(prev_bounds, corr){

      if(is.null(prev_bounds)){
        i <- 1
      } else {
        i <- nrow(prev_bounds) + 1
      }

      if(i > 0){
        if(any(dim(corr) != c(i, i))) stop(
          "Mismatch between previous bounds
       and correlation dimension."
        )
      }

      a_spend <- a.func(rates)
      bound <- solve.boundary(power=a_spend[i],
                              corr=corr, u_k=prev_bounds)

      return(bound)
    }
  } else {
    get.boundary <- function(corr){
      bounds <- get.bound.by.corr(corr, obf=(a.type=="obf"))
      return(bounds)
    }
  }

}
