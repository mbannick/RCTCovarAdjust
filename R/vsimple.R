rm(list=ls())
source("R/trial-funcs.R")
source("R/ancova.R")
source("R/boundaries.R")
source("R/trial-data.R")

try.trial <- function(K){

  t_k <- 1:K/K
  p_k <- t_k

  bound_obf <- c(1.9600, 2.7965, 3.4711, 4.0486, 4.5617,
                 5.0283, 5.5490, 5.8611, 6.2395, 6.5981,
                 6.9396, 7.2663, 7.5799, 7.8820, 8.1736)
  bound_obf <- bound_obf[K] / sqrt(rep(1:K))

  a.func.poc <- function(a, t) a * log(1 + (exp(1) - 1) * t)
  a.func.obf <- function(a, t) 4 * (1 - pnorm(qnorm(1-a/4)/sqrt(t)))

  run.full.trial <- function(theta, gamma,
                             N, n_sims, a.func,
                             estimate_sigma=TRUE, do_ancova=TRUE) {

    df <- sim.ancova.partial(t_k, p_k, N=N,
                             delta=1, theta=theta, gamma=gamma,
                             sigma2=1, b0=0)

    cat(".")
    bounds <- c()
    for(i in 1:K){

      ancova <- i == K
      ancova <- ancova & do_ancova
      bounds <- fit.stage(df, stage=i, bounds=bounds, ancova=ancova,
                          estimate_sigma=estimate_sigma,
                          a.func=a.func, a=0.05,
                          rates=t_k[1:i], N=N, n_sims=n_sims)
    }

    return(bounds)
  }

  set.seed(10)
  trial <- run.full.trial(gamma=c(1), theta=c(1),
                          estimate_sigma=TRUE,
                          do_ancova=FALSE,
                          N=10000,
                          n_sims=10000,
                          a.func=a.func.obf)
  trial.2 <- run.full.trial(gamma=c(1), theta=c(1),
                          estimate_sigma=TRUE,
                          do_ancova=TRUE,
                          N=10000,
                          n_sims=10000,
                          a.func=a.func.obf)
  print(bound_obf)
  print(trial)
  return(cbind(trial, trial.2))
}

TRIALS <- sapply(2:5, try.trial)
