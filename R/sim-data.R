library(magrittr)

#' Generate a function that will simulate data with fixed
#' simulation parameters with sample size as the argument.
#'
#' @example
#' sim <- sim.data.closure(1, 2, 0.4)
#' sim(1000)
sim.data.closure <- function(delta, n_cov, rho){

  # Solve for the (single) coefficient based on rho and n_cov
  beta <- rep(sqrt((1-rho**2)/n_cov), n_cov)

  sim.data <- function(n){
    # intercept
    int <- rep(1, n)

    # n covariates
    Z <- matrix(rnorm(n*n_cov, mean=0, sd=1), ncol=n_cov)

    # treatment indicator
    t <- rbinom(n, size=1, prob=0.5) %>% ifelse(1, -1)

    # noise
    e <- rnorm(n, 0, sd=rho)

    # outcome -- no intercept
    y <- int + t * delta + Z %*% beta + e %>% matrix
    colnames(y) <- "y"
    colnames(Z) <- paste0("cov_", 1:n_cov)

    # design matrix
    X <- cbind(t, int, Z)

    return(list(X=X, y=y))
  }
  return(sim.data)
}

sim.trial.closure <- function(N, rates){
  sim.trial <- function(sim.data){
    n_cuml <- N * rates
    n_k <- c(n_cuml[1], diff(n_cuml))

    dat <- lapply(n_k, sim.data)

    datX <- lapply(dat, function(x) x$X)
    daty <- lapply(dat, function(x) x$y)

    cumlX <- lapply(1:length(dat), function(x) do.call(rbind, datX[1:x]))
    cumly <- lapply(1:length(dat), function(x) do.call(rbind, daty[1:x]))

    dat_list <- mapply(function(X, y) list(X=X, y=y),
                       cumlX, cumly, SIMPLIFY=FALSE)

    return(dat_list)
  }
  return(sim.trial)
}
