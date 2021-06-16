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
#' func <- OBF.SPEND(0.05)
#' rates <- 1:4/4
#'
#' bound.func <- get.boundary.closure(func, rates)
#' prev_bounds <- c(4.332634, 2.963132, 2.359044)
#' corr <- corr.mat(rates)
#' bound.func(prev_bounds, corr)
get.boundary.closure <- function(a.func, rates){

  get.boundary <- function(prev_bounds, corr){

    i <- length(prev_bounds) + 1

    if(any(dim(corr) != c(i, i))) stop(
      "Mismatch between previous bounds
       and correlation dimension."
    )

    a_spend <- a.func(rates)
    bound <- solve.boundary(power=a_spend[i],
                            corr=corr, u_k=prev_bounds)

    return(bound)
  }

}
