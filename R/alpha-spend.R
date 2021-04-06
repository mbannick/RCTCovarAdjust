library(magrittr)
library(data.table)

#' THESE ARE FUNCTIONS TO FACILITATE
#' THE HYPOTHESIS TESTING IN AN ONGOING
#' TRIAL WITH ALPHA SPENDING

#' Function to generate example ANCOVA data for an alpha
#' spending hypothesis testing trial.
#'
#' @example
#' sim.ancova(10, 1, c(1, 2), 1, 1, 1)
sim.ancova <- function(n, delta, gamma, theta, sigma2, b0){

  # intercept
  int <- rep(1, n)

  # n covariates
  n_cov <- length(gamma)
  Z <- matrix(rnorm(n*n_cov, mean=0, sd=theta), ncol=n_cov)

  # treatment indicator
  t <- rbinom(n, size=1, prob=0.5) %>% ifelse(1, -1)

  # noise
  e <- rnorm(n, 0, sigma2**0.5)

  # outcome
  y <- int * b0 + t * delta + Z %*% gamma + e %>% matrix
  colnames(y) <- "y"
  colnames(Z) <- paste0("cov_", 1:n_cov)

  # design matrix
  X <- cbind(t, int, Z)

  return(cbind(y, X))
}

#' Function to simulate partial observations
#' from ANCOVA, i.e. based on an information fraction
#' and how many people have ANCOVA measurements.
#'
#' The assumption is that once you have information
#' on a subject, you will not lose information on them.
#' You do not necessarily have to have 100% of the information
#' for covariates by the end of the trial.
#'
#' @param t_k Information rates at each time point
#' @param p_k Covariate information rates at each time point
#' @param N Total number of subjects
#' @param ... Additional arguments passed to \code{sim.ancova}
#'
#' @examples
#' t_k <- 1:5/5
#' p_k <- (1:5/5)**2
#'
#' sim.ancova.partial(t_k, p_k, N=100,
#'                    delta=1, theta=1, gamma=c(0.5, 2),
#'                    sigma2=1, b0=0)
sim.ancova.partial <- function(t_k, p_k, N, ...){

  # Simulate ANCOVA data
  data <- sim.ancova(n=N, ...) %>% data.table
  id <- 1:N
  data[, id := id]
  cov_cols <- grep("cov", colnames(data))

  # Get number of subjects per stage
  # based on the information rate
  n_k <- t_k * N
  n_k_diff <- c(n_k[1], diff(n_k))

  c_k <- p_k * N

  final <- data.table()
  has_cov <- c()
  no_cov <- c()

  # Loop through stages and add data
  # in long format
  for(k in 1:length(n_k)){

    # Get the ids of those already in the trial
    past_ids <- c()
    if(nrow(final) > 0) past_ids <- final[, id] %>% unique()

    # Get the ids of all enrolled in trial at stage k
    all_ids <- id[1:n_k[k]]

    # Get ids of only those for this stage
    k_ids <- all_ids[!all_ids %in% past_ids]

    # Needs covariates
    needs_cov <- c(no_cov, k_ids)
    num_sample <- c_k[k] - length(has_cov)

    # Figure out whether or not they have covariates
    has_cov <- c(has_cov, sample(needs_cov, size=num_sample))
    no_cov <- all_ids[!all_ids %in% has_cov]

    # Subset out the data and
    # get rid of covariates if the
    # sample doesn't have covariates yet
    dat <- data[all_ids, ]
    dat[no_cov, cov_cols] <- NA
    dat[, look := k]

    final <- rbind(final, dat)
  }
  final[, enter := lapply(.SD, min), .SDcols="look", by="id"]
  return(final)
}

fit.anova <- function()
