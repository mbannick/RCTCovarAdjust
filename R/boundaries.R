library(magrittr)

#' Validate parameters to use in the type I error functions.
#' Don't want to include this in the function itself because
#' it would be run multiple times when solving for the boundaries.
#'
#' @param n_k Sample size or sequence of sample sizes for each stage.
#' @param stat.func Function (or list of functions)
#'   to calculate the test statistic based on the generated data. If this is a list
#'   of functions, it needs to be the same length as the sequence of boundaries.
#' @param data.generator Function to generate data. The output of this function
#'   should be a matrix, which is the input for stat.func.
#'   First argument of this function needs to be n, the sample size, and the function
#'   needs to be vectorized.
#' @param ... Additional arguments to the \code{data.generator} function
#'
#' @examples
#' validate.trial.params(n_k=rep(10, 5), stat.func=mean, data.generator=rnorm)
#'
validate.trial.params <- function(n_k, stat.func, data.generator, ...) {

  # CHECK FOR STAGES
  K <- length(n_k)
  if(K == 1) stop("Must have at least two stages.")

  # CHECK TEST STATISTIC FUNCTIONS
  if(length(stat.func) > 1){
    if(length(stat.func) != K) stop("Sequence of test statistic functions
                                  must be of the same length as the sequence of boundaries.")
  } else {
    stat.func <- list(stat.func)
  }

  # CHECK IF DATA GENERATOR AND TEST STATISTIC ARE COMPATIBLE
  x <- data.generator(n=10, ...)
  msg <- "Data generator and test statistic are not compatible."
  for(stat in stat.func){
    tryCatch(stat(x),
             error=function(c) stop(msg),
             warning=function(c) stop(msg))
  }
}

#' Creates simulations of group sequential trial data and
#' calculates test statistic based on user-defined stat function.
#'
#' @export
#' @param n_sims Number of simulations
#' @param n_k Sample size per stage
#' @param stat.func Function (or list of functions)
#'   to calculate the test statistic based on the generated data, at each
#'   observed information rates. If this is a list
#'   of functions, it needs to be the same length as the sequence of boundaries.
#' @param data.generator Function to generate data. The output of this function
#'   should be a matrix, which is the input for stat.func.
#'   First argument of this function needs to be n, the sample size, and the function
#'   needs to be vectorized.
#' @param ... Additional arguments to be passed to \code{data.generator}
#' @examples
#' create.simulations(10, n_k=c(10, 10, 10),
#'                    stat.func=list(mean, sum, mean),
#'                    data.generator=function(n) as.matrix(rnorm(n)))
create.simulations <- function(n_sims, n_k, stat.func, data.generator, ...){

  # Set up the sample sizes
  K <- length(n_k)
  n <- sum(n_k)
  n_cuml <- cumsum(n_k)

  if(class(stat.func) == "function") stat.func <- list(stat.func)

  test.stat <- function() {
    # Generate data and calculate test statistics
    data <- data.generator(n=n, ...)

    # Calculate test statistics with a list of test statistic functions
    # and a list of cumulative sample sizes to subset
    apply.it <- function(f, x) f(data[1:x, ])
    t.stats <- mapply(apply.it, f=stat.func, x=n_cuml)

    return(t.stats)
  }

  results <- replicate(n_sims, test.stat()) %>% t
  return(results)
}

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
#' @param stat.func Function (or list of functions)
#'   to calculate the test statistic based on the generated data, at each
#'   observed information rates. If this is a list
#'   of functions, it needs to be the same length as the sequence of boundaries.
#' @param data.generator Function to generate data. The output of this function
#'   should be a matrix, which is the input for stat.func.
#'   First argument of this function needs to be n, the sample size, and the function
#'   needs to be vectorized.
#'
#' @examples
#' set.seed(101)
#' # replicate Pocock boundaries
#' get.boundaries(power=0.05, n_sims=1000, K=2, n=100,
#'                stat.func=function(x) mean(x) * sqrt(length(x)),
#'                data.generator=function(n) matrix(rnorm(n)))
#'
#' # replicate OBF boundaries
#' # notice the sqrt(2) in the stat.list
#' z.stat <- function(x) mean(x) * sqrt(length(x))
#' stat.list <- list(z.stat, function(...) z.stat(...) * sqrt(2))
#' get.boundaries(power=0.05, n_sims=1000, K=2, n=100,
#'                stat.func=stat.list,
#'                data.generator=function(n) matrix(rnorm(n)))
get.boundaries <- function(power, n_sims, K, n, stat.func, data.generator, ...){

  # Create sample size vector
  n_k <- rep(n, K)
  # Validate parameters
  validate.trial.params(n_k, stat.func, data.generator, ...)

  simulations <- create.simulations(n_sims=n_sims, n_k=n_k,
                                    stat.func=stat.func,
                                    data.generator=data.generator, ...)

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
#' @param stat.func Function (or list of functions)
#'   to calculate the test statistic based on the generated data, at each
#'   observed information rates. If this is a list
#'   of functions, it needs to be the same length as the sequence of boundaries.
#' @param data.generator Function to generate data. The output of this function
#'   should be a matrix, which is the input for stat.func.
#'   First argument of this function needs to be n, the sample size, and the function
#'   needs to be vectorized.
#'
#' @examples
#' set.seed(101)
#' # information rates
#' t <- 1:3/3
#'
#' # approximate Pocock boundaries
#' a.func.pocock <- function(a, t) a * log(1 + (exp(1) - 1) * t)
#' get.boundaries.aspend(a.func=a.func.pocock, a=0.05,
#'                       rates=t, N=1000, n_sims=1000,
#'                       stat.func=function(x) mean(x) * sqrt(length(x)),
#'                       data.generator=function(n) matrix(rnorm(n)))
#'
#' # approximate O'Brien-Fleming boundaries
#' a.func.obf <- function(a, t) 4 * (1 - pnorm(qnorm(1-a/4)/sqrt(t)))
#' get.boundaries.aspend(a.func=a.func.obf, a=0.05,
#'                       rates=t, N=1000, n_sims=1000,
#'                       stat.func=function(x) mean(x) * sqrt(length(x)),
#'                       u_k=c(3.35, 2.47),
#'                       data.generator=function(n) matrix(rnorm(n)))
get.boundaries.aspend <- function(a.func, a, rates, N, n_sims,
                                  stat.func, data.generator, u_k=c(), ...){

  # Get sample size increments
  K <- length(rates)
  n_k <- round(rates*N)
  n_k <- c(n_k[1], diff(n_k))

  # Number of *fixed* previous bounds
  K_prev <- length(u_k)

  # Validate parameters
  validate.trial.params(n_k, stat.func, data.generator, ...)

  # Create the simulations
  simulations <- create.simulations(n_sims=n_sims, n_k=n_k,
                                    stat.func=stat.func,
                                    data.generator=data.generator, ...)

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
