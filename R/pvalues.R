library(magrittr)
library(MASS)
library(mvtnorm)
# source("~/repos/RCTCovarAdjust/R/covariance.R")
# source("~/repos/RCTCovarAdjust/R/boundaries.R")

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
    return(lower[1])
  } else if(type == "upper"){
    return(upper[1])
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

  suppressWarnings(
    result <- pmvnorm(lower=bounds[, 1], upper=bounds[, 2],
                      ..., algorithm=Miwa(steps=1000))
  )

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

# #' Get the probability of rejecting at stage k.
# #' We need this to return a two-element vector because
# #' for the orderings, crossing lower *or* upper boundaries
# #' is important.
# #'
# #' @param u_k Vector of boundaries up through stage k
# #' @param corr Correlation matrix for test statistics
# #'
# #' @examples
# #' u_k <- c(4.332634, 2.963132, 2.359044)
# #' u_k1 <- cbind(-u_k, u_k)
# #' u_k2 <- cbind(rep(-Inf, 3), u_k)
# #' n_k <- c(10, 20, 30)
# #' corr.1 <- corr.mat(n_k)
# #' corr.2 <- corr.mat(n_k, rho=0.5, mis=c(F, F, T))
# #' .reject.prob.k(u_k=u_k1, corr=corr.1)
# #' .reject.prob.k(u_k=u_k2, corr=corr.1)
.reject.prob.k <- function(u_k, corr, alt){

  if(!"matrix" %in% class(u_k)) u_k <- matrix(u_k, nrow=1)

  K <- nrow(u_k)

  if(K == 1){
    p.lower <- pnorm(u_k[K, 1], mean=alt, lower.tail=T)
    p.upper <- pnorm(u_k[K, 2], mean=alt, lower.tail=F)
  } else {
    prev <- u_k[1:(K-1), ]

    lower.K <- u_k[K, 1]
    upper.K <- u_k[K, 2]

    lower.block <- rbind(prev, c(-Inf, lower.K))
    upper.block <- rbind(prev, c(upper.K, Inf))

    p.lower <- .pmvnorm.list(blocks=lower.block, corr=corr, mean=alt)
    p.upper <- .pmvnorm.list(blocks=upper.block, corr=corr, mean=alt)
  }
  return(c(p.lower, p.upper))
}

.one.test.k <- function(obs, u_k, alt, corr=matrix(1), last_stage=TRUE){
  if(!"matrix" %in% class(u_k)) u_k <- matrix(u_k, nrow=1)

  K <- nrow(u_k)

  if(nrow(corr) != K) stop("Correlation matrix needs to be the same
                           dimension as the boundaries.")
  if(K == 1){
    if(last_stage){
      p.lower <- pnorm(obs, mean=alt, lower.tail=T)
      p.upper <- pnorm(obs, mean=alt, lower.tail=F)
    } else {
      # Stop in this stage and more negative, or do not stop

      # Probability of not stopping
      p.notstop <- pnorm(u_k[1, 2], mean=alt, lower.tail=T) -
                   pnorm(u_k[1, 1], mean=alt, lower.tail=T)

      # Stopping in this stage and more negative
      p.L.neg <- pnorm(min(obs, u_k[1, 1]), mean=alt, lower.tail=T)
      p.L.pos <- pnorm(max(obs, u_k[1, 2]), mean=alt, lower.tail=T) -
                 pnorm(u_k[1, 2], mean=alt, lower.tail=T)


      # Stop in this stage, and more positive
      p.U.neg <- pnorm(min(obs, u_k[1, 1]), mean=alt, lower.tail=F) -
                 pnorm(u_k[1, 1], mean=alt, lower.tail=F)
      p.U.pos <- pnorm(max(obs, u_k[1, 2]), mean=alt, lower.tail=F)

      p.lower <- p.L.neg + p.L.pos
      p.upper <- p.U.neg + p.U.pos
      if(obs <= u_k[1, 1]) p.upper <- p.upper + p.notstop
      if(obs >= u_k[1, 2]) p.lower <- p.lower + p.notstop
    }
  } else {

    # Previous boundaries
    prev <- u_k[1:(K-1), ]

    if(last_stage){

      lower.blocks <- rbind(prev, c(-Inf, obs))
      upper.blocks <- rbind(prev, c(obs, Inf))

    } else {

      # Stop in this stage and more negative
      p.L.neg <- rbind(prev, c(-Inf, min(obs, u_k[K, 1])))
      p.L.pos <- rbind(prev, c(u_k[K, 2], max(obs, u_k[K, 2])))

      # Stop in this stage, and more positive
      p.U.neg <- rbind(prev, c(min(u_k[K, 1], obs), u_k[K, 1]))
      p.U.pos <- rbind(prev, c(max(obs, u_k[K, 2]), Inf))

      # Do not stop
      p.notstop <- rbind(prev, u_k[K,])

      if(obs <= u_k[K, 1]){
        lower.blocks <- list(p.L.neg, p.L.pos)
        upper.blocks <- list(p.U.neg, p.U.pos, p.notstop)
      } else if(obs >= u_k[K, 2]) {
        lower.blocks <- list(p.L.neg, p.L.pos, p.notstop)
        upper.blocks <- list(p.U.neg, p.U.pos)
      } else {
        lower.blocks <- list(p.L.neg, p.L.pos)
        upper.blocks <- list(p.U.neg, p.U.pos)
      }
    }
    p.lower <- .pmvnorm.list(blocks=lower.blocks, corr=corr, mean=alt)
    p.upper <- .pmvnorm.list(blocks=upper.blocks, corr=corr, mean=alt)
  }

  return(c(p.lower, p.upper))
}

.two.tests.k <- function(obs, u_k, alt, corr, crossed_lower=FALSE){

  if(!"matrix" %in% class(u_k)) u_k <- matrix(u_k, nrow=1)
  K <- nrow(u_k)

  if(nrow(corr) != K+1) stop("Correlation matrix needs to be
                             higher dimension than boundaries.")

  if(K > 1){
    prev <- u_k[1:(K-1), ]
  } else {
    prev <- NULL
  }

  l.K <- u_k[K, 1] # lower bound at this stage
  u.K <- u_k[K, 2] # upper bound at this stage

  # Stop in this stage with MONITOR and more negative TEST, or do not stop
  lower.cross <- c(-Inf, l.K)
  upper.cross <- c(u.K, Inf)
  no.cross <- c(l.K, u.K)

  # More negative TEST
  lower.tail <- c(-Inf, obs)
  upper.tail <- c(obs, Inf)

  lower.block.1 <- rbind(prev, lower.cross, lower.tail)
  upper.block.1 <- rbind(prev, upper.cross, upper.tail)

  lower.block.2 <- rbind(prev, upper.cross, lower.tail)
  upper.block.2 <- rbind(prev, lower.cross, upper.tail)

  nocross.block <- rbind(prev, no.cross, c(-Inf, Inf))

  if(crossed_lower){
    lower.blocks <- list(lower.block.1, lower.block.2)
    upper.blocks <- list(upper.block.1, upper.block.2, nocross.block)
  } else {
    lower.blocks <- list(lower.block.1, lower.block.2, nocross.block)
    upper.blocks <- list(upper.block.1, upper.block.2)
  }

  p.lower <- .pmvnorm.list(blocks=lower.blocks, corr=corr, mean=alt)
  p.upper <- .pmvnorm.list(blocks=upper.blocks, corr=corr, mean=alt)

  return(c(p.lower, p.upper))
}

#' Get a stage-wise p-value after the trial is complete.
#'
#' @param obs Observed z-statistic
#' @param u_k Matrix of boundaries for *previous* stages (if last_stage = F)
#'            There should be two columns, one for lower bound, one for upper.
#'            This allows asymmetric boundaries.
#' @param k_r Stage trial ended
#' @param rho sqrt(R^2)
#' @param ancova_monitor Whether to monitor with ANCOVA
#' @param ancova_test Whether to get final inferene with ANCOVA
#' @param crossed_lower Was it the lower or upper boundary that was crossed
#'                      if there was early stopping (if not, arg ignored).
#' @param type Type of p-value ("two-sided", "lower", or "upper")
#' @export
#'
#' @examples
#' # OBF boundaries in 3-stage trial
#' bounds <- get.boundaries.design(rates=1:3/3, obf=FALSE)
#' u_k <- cbind(-bounds, bounds)
#' n_k <- c(10, 20, 30)
#'
#' # Hit the last boundary at the last stage, should be p-value = 0.05
#' get.pvalue.sw(obs=bounds[3],
#'               u_k=u_k, k_r=3, n_k=n_k, last_stage=T,
#'               ancova_monitor=F, ancova_test=F)
#'
#' # Switch to using ANCOVA at last stage
#' get.pvalue.sw(obs=bounds[3],
#'               u_k=u_k, k_r=3, n_k=n_k, last_stage=T,
#'               ancova_monitor=F, ancova_test=T, rho=sqrt(0.5))
#'
#' # What if we had been using corrected boundaries
#' bounds.c <- get.boundaries.design(rates=1:3/3, obf=FALSE,
#'                                   rho=sqrt(0.5), change=c(0, 0, 1))
#' u_k <- cbind(-bounds.c, bounds.c)
#' get.pvalue.sw(obs=bounds.c[3],
#'               u_k=u_k, k_r=3, n_k=n_k, last_stage=T,
#'               ancova_monitor=F, ancova_test=T, rho=sqrt(0.5))
get.pvalue.sw <- function(obs, u_k, n_k, k_r,
                          alt=0,
                          rho=1,
                          ancova_monitor=F,
                          ancova_test=T,
                          last_stage=F,
                          crossed_lower=FALSE,
                          type="two-sided"){
  .check.type(type)

  # If rho is 1, that means that ANCOVA == ANOVA.
  if(rho == 1.){
    ancova_test <- F
    ancova_monitor <- F
  }
  # Indicator for switching between ANOVA and ANCOVA
  # at the last stage.
  switch <- ancova_monitor != ancova_test

  if(length(alt) == 1){
    if(switch & !last_stage){
      alt <- rep(alt, k_r + 1)
    } else {
      alt <- rep(alt, k_r)
    }
  }

  if(k_r < 1) stop("Need to have rejected at stage 1 or higher.")

  if(nrow(u_k) != k_r) stop("There need to be as many bounds
                            as the stage at which you rejected, k_r.")

  if(length(n_k) != k_r) stop("There need to be as many sample sizes
                            as the stage at which you rejected, k_r.")
  if(k_r == 1){
    # Rejected at the first stage
    if((!switch) | last_stage){
      final.p <- .one.test.k(obs=obs, u_k=u_k[1,], alt=alt)
    } else {
      # Switch between methods
      corr <- corr.mat(rep(n_k[1], 2), rho=rho, mis=c(F, T))
      final.p <- .two.tests.k(obs=obs, u_k=u_k[1,], corr=corr, alt=alt,
                              crossed_lower=crossed_lower)
    }
  } else {
    # Rejected at some time beyond the first stage
    final.p <- c(0, 0)

    for(i in 1:k_r){
      if(i == 1){
        p.i <- .reject.prob.k(u_k[1, ], corr=matrix(1), alt=alt[1])
      } else if(i < k_r){
        corr <- corr.mat(n_k[1:i])
        p.i <- .reject.prob.k(u_k[1:i,], corr=corr, alt=alt[1:i])
      } else {
        if(last_stage){
          if(switch){
            corr <- corr.mat(n_k[1:k_r], rho=rho, mis=c(rep(F, k_r-1), T))
          } else {
            corr <- corr.mat(n_k[1:k_r], rho=rho, mis=rep(F, k_r))
          }
          p.i <- .one.test.k(obs=obs, u_k=u_k, corr=corr, alt=alt[1:i], last_stage=TRUE)
        } else {
          if(switch){
            corr <- corr.mat(c(n_k[1:k_r], n_k[k_r]), rho=rho, mis=c(rep(F, k_r), T))
            p.i <- .two.tests.k(obs=obs, u_k=u_k, corr=corr, alt=alt[1:(i+1)],
                                crossed_lower=crossed_lower)
          } else {
            corr <- corr.mat(c(n_k[1:k_r]))
            p.i <- .one.test.k(obs=obs, u_k=u_k, corr=corr, alt=alt[1:i], last_stage=FALSE)
          }
        }
      }
      final.p <- final.p + p.i
    }
  }
  val <- .get.pval.from.type(type, lower=final.p[1], upper=final.p[2])
  return(val)
}

search.fun.sw <- function(est, sd_K, u_k, n_k, k_r, alpha,
                          rho, ancova_monitor, ancova_test,
                          last_stage, crossed_lower, low){

  # Get sample size at the last stage
  n_K <- n_k[length(n_k)]

  # Function to translate effect size into z-statistic
  # At the analysis stage K (not at the first stage)
  switch <- ancova_monitor != ancova_test
  if((rho < 1.0) & switch){

    sd_test <- sd_K
    sd_monitor <- ifelse(ancova_test, sd_K/rho, sd_K*rho)
    if(!last_stage){
      # Add another of the sample sizes on to the end
      n_k_forz <- c(n_k, n_K)
      sd_k_forz <- c(rep(sd_monitor, k_r), sd_test)
    } else {
      n_k_forz <- n_k
      sd_k_forz <- c(rep(sd_monitor, k_r-1), sd_test)
    }
  } else {
    n_k_forz <- n_k
    sd_k_forz <- rep(sd_K, k_r)
  }
  get.z <- function(eff) sqrt(n_k_forz) * (eff) / sd_k_forz

  # Observed test statistic is
  obs.z <- sqrt(n_K) * (est) / sd_K

  # Create a function to translate the estimate to a z-statistic
  search.fun <- function(eff){

    test.z <- get.z(eff)

    if(low){
      p <- get.pvalue.sw(obs.z, alt=test.z, u_k=u_k, n_k=n_k, k_r=k_r, rho=rho,
                         ancova_monitor=ancova_monitor,
                         ancova_test=ancova_test,
                         last_stage=last_stage,
                         crossed_lower=crossed_lower,
                         type="lower")
    } else {
      p <- get.pvalue.sw(obs.z, alt=test.z, u_k=u_k, n_k=n_k, k_r=k_r, rho=rho,
                         ancova_monitor=ancova_monitor,
                         ancova_test=ancova_test,
                         last_stage=last_stage,
                         crossed_lower=crossed_lower,
                         type="upper")
    }
    return(p - alpha)
  }
  val <- uniroot(search.fun, lower=-500, upper=500, trace=1)$root
  # out.dir <- "~/repos/RCTCovarAdjust/debug2/"
  # num <- length(list.files(out.dir))
  # pdf(paste0(out.dir, "sim-", num+1, ".pdf"), height=6, width=6)
  # xs <- seq(-0.7, 0.7, by=0.01)
  # ps.L <- sapply(xs, function(x) get.pvalue.sw(obs.z, alt=get.z(x),
  #                                              u_k=u_k, n_k=n_k, k_r=k_r, rho=rho,
  #                                              ancova_monitor=ancova_monitor,
  #                                              ancova_test=ancova_test,
  #                                              last_stage=last_stage,
  #                                              crossed_lower=crossed_lower,
  #                                              type="lower"))
  # ps.U <- sapply(xs, function(x) get.pvalue.sw(obs.z, alt=get.z(x),
  #                                              u_k=u_k, n_k=n_k, k_r=k_r, rho=rho,
  #                                              ancova_monitor=ancova_monitor,
  #                                              ancova_test=ancova_test,
  #                                              crossed_lower=crossed_lower,
  #                                              last_stage=last_stage, type="upper"))
  # plot(ps.L ~ xs, type='l')
  # lines(ps.U ~ xs, col='blue')
  # abline(v=0.1)
  # dev.off()
  return(val)
}

# #' Get p-value using the sample mean ordering.
# #' Needs to take in the obs as a standardized sample mean.
# #'
# #' @param obs Standardized sample mean observation
# #' @param u_K Boundaries for all stages except the last stage K
# #' @param n_K Sample sizes for all stages up through max stage K
# #'            This is needed to compute the variance of the sample mean.
# #' @param type Type of p-value (upper, lower, two-sided)
# #'
# #' @examples
# #' # OBF boundaries
# #' u_k <- c(4.332634, 2.963132, 2.359044, 2.014090)
# #' u_k1 <- cbind(-u_k, u_k)
# #' u_k2 <- cbind(rep(-Inf, 4), u_k)
# #'
# #' # This gives exactly alpha = 0.05 by putting in the last
# #' # boundary.
# #' n_K <- c(10, 20, 30, 40)
# #'
# #' get.pvalue.sm(obs=-2.01409/sqrt(40), u_K=u_k1[1:3,], n_K=n_K,
# #'               ancova_monitor=F, ancova_test=F, rho=1.0)
# #' get.pvalue.sm(obs=-2.01409/sqrt(40), u_K=u_k1[1:3,], n_K=n_K,
# #'               ancova_monitor=F, ancova_test=T, rho=0)
# #' # 0.05000002
# #' get.pvalue.sm(obs=-2.01409/sqrt(40), u_K=u_k1[1:3,], n_K=n_K,
# #'               ancova_monitor=F, ancova_test=T, rho=0.4)
# #' get.pvalue.sm(obs=-2.01409/sqrt(40), u_K=u_k1[1:3,], n_K=n_K,
# #'               ancova_monitor=F, ancova_test=T, rho=0.8)
# #' for(i in seq(0.01, 0.9, 0.1)){
# #'     print(get.pvalue.sm(obs=-2.01409/sqrt(40), u_K=u_k1[1:3,], n_K=n_K,
# #'               ancova_monitor=F, ancova_test=T, rho=i))
# #' }
# #'
# #'
# #' get.pvalue.sm(obs=-5/sqrt(40), u_K=u_k1[1:3,], n_K=n_K,
# #'               ancova_monitor=F, ancova_test=F, rho=1.0)
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
# get.pvalue <- function(obs, u_K, k, corr, n_K=c(),
#                        ordering="stage-wise"){
#
#   if(ordering == "sample mean" & length(n_K) == 0){
#     stop("If using the stage-wise ordering, need to give the sample size
#           for all data looks, including those planned at future stages.")
#   }
#
#   if(ordering == "stage-wise"){
#     p <- get.pvalue.sw(obs=obs, u_k=u_K[1:k], corr=corr[1:k, 1:k])
#   } else if(ordering == "sample mean"){
#     p <- get.pvalue.sm(obs=obs, u_K=u_K, corr=corr, n_K=n_K)
#   }
#
#   return(p)
# }
