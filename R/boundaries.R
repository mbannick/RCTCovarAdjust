library(magrittr)
library(MASS)

#' Calculate the type I error based on a sequence of boundaries.
#' Inputs include (1) trial parameters, (2) data generating procedure parameters,
#' and (3) a function to compute the test statistic based on observations and stage.
#'
#' @param u Constant boundary
#' @param simulations Output from create.simulations function.
#'   Should be long by simulation and wide by stage.
#' @param previous A vector of previous conditions.
#'
#' @examples
#' simulations <- matrix(rnorm(100), nrow=100)
#' get.power(2.7, simulations=simulations)
#'
#' previous <- rep(TRUE, 10)
get.power <- function(u, simulations, previous=NULL){

  # Check to see if statistics exceeded the boundaries,
  # for each simulation and stage
  reject <- abs(simulations) >= u

  # If there are previous conditions to add to this
  # probability, then do so here
  if(!is.null(previous)) reject <- previous & reject

  # For each simulation, check if it exceeded at any stage and
  # take the average, that's the power
  reject <- apply(reject, MARGIN=1, FUN=any) %>% mean
  return(reject)
}

#' Root solve for a particular alpha level.
#'
#' @param power The type I error (or power generally) to solve for
#' @param simulations Output from create.simulations function.
solve.boundary <- function(power, simulations, previous=NULL){

  # Solve for the u that results in desired power
  u.search <- function(u) get.power(u, simulations, previous=previous) - power
  u <- uniroot(u.search, lower=0, upper=50)$root

  return(u)
}

#' Derive constant boundaries.
#'
#' Get the boundaries required for a particular power.
#' The boundaries are constant across stage, *for a particular stat.func*.
#'
#' E.g. if you want Pocock boundaries, you would have a likelihood ratio stat.func
#' E.g. if you want OBF boundaries, you would have a score statistic stat.func.
#'
#' @export
#' @param power The type I error desired (or power)
#' @param n Sample size per stage
#' @param n_sims Number of Monte-Carlo simulations
#' @param K Total number of stages
#' @param sim.generator A generator for simulated test statistics
#'
#' @examples
#' # replicate Pocock boundaries
#' sim.generator <- function(n_sims, n_k) mvrnorm(n=n_sims,
#'                                                mu=rep(0, length(n_k)),
#'                                                Sigma=basic.cov(n_k))
#' # this should be close to 2.1783
#' get.boundaries(power=0.05, n_sims=100000, K=2, n=1000,
#'                sim.generator=sim.generator)
get.boundaries <- function(power, n_sims, K, n, sim.generator, ...){

  # Create sample size vector
  n_k <- rep(n, K)

  simulations <- sim.generator(n_sims=n_sims, n_k=n_k, ...)
  u <- solve.boundary(power, simulations)

  return(u)
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
#' @param sim.generator A generator for simulated test statistics
#'
#' @examples
#' set.seed(101)
#' # information rates
#' t <- 1:3/3
#'
#' sim.generator <- function(n_sims, n_k) mvrnorm(n=n_sims,
#'                                                mu=rep(0, length(n_k)),
#'                                                Sigma=basic.cov(n_k))
#'
#' # approximate Pocock boundaries
#' a.func.pocock <- function(a, t) a * log(1 + (exp(1) - 1) * t)
#' get.boundaries.aspend(a.func=a.func.pocock, a=0.05,
#'                       rates=t, N=1000, n_sims=100000,
#'                       sim.generator=sim.generator)
#'
#' # approximate O'Brien-Fleming boundaries
#' a.func.obf <- function(a, t) 4 * (1 - pnorm(qnorm(1-a/4)/sqrt(t)))
#' get.boundaries.aspend(a.func=a.func.obf, a=0.05,
#'                       rates=t, N=100000, n_sims=10000,
#'                       sim.generator=sim.generator,
#'                       u_k=c(3.471, 2.454))
get.boundaries.aspend <- function(a.func, a, rates, N, n_sims,
                                  u_k=c(), sim.generator, ...){

  # Get sample size increments
  K <- length(rates)
  n_k <- round(rates*N)
  n_k <- c(n_k[1], diff(n_k))

  # Number of *fixed* previous bounds
  K_prev <- length(u_k)

  # Create the simulations
  simulations <- sim.generator(n_sims=n_sims, n_k=n_k, ...)

  # Create alpha spending function
  a.spend <- function(t) a.func(a=a, t)

  # Append 0 onto the rates
  if(!0 %in% rates) rates <- c(0, rates)
  a.diff <- a.spend(rates) %>% diff

  # Create a "previous" vector that we will use to
  # indicate having not rejected at previous stages
  bounds <- c()
  previous <- rep(TRUE, n_sims)

  # For each stage, solve for the boundary, unless
  # you passed in fixed boundaries
  for(i in 1:K){
    sim <- simulations[, i] %>% as.matrix
    if(i > K_prev){
      bound <- solve.boundary(power=a.diff[i],
                              simulations=sim,
                              previous=previous)
    } else {
      bound <- u_k[i]
    }

    # Update the previously rejected vector
    # with information from this stage
    previous <- previous & (abs(sim) < bound)
    bounds <- c(bounds, bound)
  }

  return(bounds)
}
