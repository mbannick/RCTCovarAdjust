rm(list=ls())
source("R/trial-funcs.R")
source("R/ancova.R")
source("R/boundaries.R")
source("R/trial-data.R")

K <- 2

t_k <- 1:K/K
p_k <- (1:K/K)**2

a.func.poc <- function(a, t) a * log(1 + (exp(1) - 1) * t)
a.func.obf <- function(a, t) 4 * (1 - pnorm(qnorm(1-a/4)/sqrt(t)))

bound_poc <- rep(2.413, K) # THIS ISNT CORRECT
bound_obf <- c(1.9600, 2.7965, 3.4711, 4.0486, 4.5617,
               5.0283, 5.5490, 5.8611, 6.2395, 6.5981,
               6.9396, 7.2663, 7.5799, 7.8820, 8.1736)
bound_obf <- bound_obf[K] / sqrt(rep(1:K))

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

n_sims <- c(10000)
N <- c(1000)
gamma <- c(0.5, 1, 5)
theta <- c(1)

params <- expand.grid(n_sims=n_sims, N=N, gamma=gamma, theta=theta)

res.anova <- list()
res.ancova <- list()

for(i in 1:nrow(params)){
  parm <- params[i, ]
  print(parm)
  set.seed(10)
  reps.anova <- replicate(5, run.full.trial(gamma=c(parm$gamma), theta=c(parm$theta),
                                             estimate_sigma=TRUE,
                                             do_ancova=FALSE, N=parm$N, n_sims=parm$n_sims,
                                             a.func=a.func.obf))
  res.anova[[i]] <- reps.anova
  set.seed(10)
  reps.ancova <- replicate(5, run.full.trial(gamma=c(parm$gamma), theta=c(parm$theta),
                                             estimate_sigma=TRUE,
                                             do_ancova=TRUE, N=parm$N, n_sims=parm$n_sims,
                                             a.func=a.func.obf))
  res.ancova[[i]] <- reps.ancova

}

pdf("exmaple-obf-2.pdf", height=50, width=8)
par(mfrow=c(nrow(params), 2))
for(i in 1:nrow(params)){

  title <- paste0(colnames(params), ": ", params[i, ]) %>% paste0(collapse = " ")

  t <- 1:K

  plot(bound_obf ~ t, type='l', ylim=c(0, 6), main=paste0(title, "\nANOVA"))
  for(j in 1:dim(res.anova[[i]])[2]){
    lines(res.anova[[i]][, j] ~ t, col='lightblue')
  }
  lines(bound_obf ~ t)

  plot(bound_obf ~ t, type='l', ylim=c(0, 6), main=paste0(title, "\nANCOVA"))
  for(j in 1:dim(res.ancova[[i]])[2]){
    lines(res.ancova[[i]][, j] ~ t, col='lightblue')
  }
  lines(bound_obf ~ t)
}
dev.off()

res.anova <- list()
res.ancova <- list()

for(i in 1:nrow(params)){
  parm <- params[i, ]
  print(parm)
  set.seed(10)
  reps.anova <- replicate(10, run.full.trial(gamma=c(parm$gamma), theta=c(parm$theta),
                                             estimate_sigma=FALSE,
                                             do_ancova=FALSE, N=parm$N, n_sims=parm$n_sims,
                                             a.func=a.func.poc))
  res.anova[[i]] <- reps.anova
  set.seed(10)
  reps.ancova <- replicate(10, run.full.trial(gamma=c(parm$gamma), theta=c(parm$theta),
                                              estimate_sigma=FALSE,
                                              do_ancova=TRUE, N=parm$N, n_sims=parm$n_sims,
                                              a.func=a.func.poc))
  res.ancova[[i]] <- reps.ancova

}

pdf("exmaple-poc-2.pdf", height=50, width=8)
par(mfrow=c(nrow(params), 2))
for(i in 1:nrow(params)){

  title <- paste0(colnames(params), ": ", params[i, ]) %>% paste0(collapse = " ")

  t <- 1:5

  plot(bound_poc ~ t, type='l', ylim=c(0, 6), main=paste0(title, "\nANOVA"))
  for(j in 1:dim(res.anova[[i]])[2]){
    lines(res.anova[[i]][, j] ~ t, col='lightblue')
  }
  lines(bound_poc ~ t)

  plot(bound_poc ~ t, type='l', ylim=c(0, 6), main=paste0(title, "\nANCOVA"))
  for(j in 1:dim(res.ancova[[i]])[2]){
    lines(res.ancova[[i]][, j] ~ t, col='lightblue')
  }
  lines(bound_poc ~ t)
}
dev.off()

# SINGLE FOR DEBUGGING

set.seed(10)
trial <- run.full.trial(gamma=c(2), theta=c(1),
                        estimate_sigma=FALSE,
                        do_ancova=FALSE,
                        N=100,
                        n_sims=10000,
                        a.func=a.func.obf)

set.seed(10)
reps.anova <- replicate(5, run.full.trial(gamma=c(2), theta=c(1),
                                          estimate_sigma=TRUE,
                                          do_ancova=FALSE,
                                          N=100,
                                          n_sims=10000,
                                          a.func=a.func.obf))
set.seed(10)
reps.ancova <- replicate(5, run.full.trial(gamma=c(100), theta=c(10),
                                          estimate_sigma=TRUE,
                                          do_ancova=TRUE,
                                          N=100,
                                          n_sims=1000,
                                          a.func=a.func.obf))




#
#
## fit.stage(df, stage=5, bounds=c(4.562, 3.226, 2.634, 2.281), ancova=F,
#           a.func=a.func.obf, a=0.05,
#           rates=t_k, N=1000, n_sims=10000)
# set.seed(10)
# reps <- replicate(10, run.full.trial(gamma=c(1, 1), theta=c(1, 1), estimate_sigma=FALSE,
#                                      do_ancova=FALSE, N=100, n_sims=1000,
#                                      a.func=a.func.obf))
#
# run.trial(n=2, gamma=c(1, 1), theta=c(1, 1), estimate_sigma=FALSE,
#           do_ancova=FALSE, N=100, n_sims=1000,
#           a.func=a.func.obf)
#
# set.seed(10)
# reps.2 <- replicate(2, run.full.trial(estimate_sigma=FALSE,
#                                     do_ancova=TRUE, theta=c(1, 1),
#                                     gamma=c(10, 1)))
#
# set.seed(10)
# reps.3 <- replicate(2, run.full.trial(estimate_sigma=TRUE,
#                                     do_ancova=TRUE, theta=c(1, 1),
#                                     gamma=c(10, 1)))
#
# bound_obf <- c(1.9600, 2.7965, 3.4711, 4.0486, 4.5617,
#                5.0283, 5.5490, 5.8611, 6.2395, 6.5981,
#                6.9396, 7.2663, 7.5799, 7.8820, 8.1736)
# bound_obf <- bound_obf[5] / sqrt(rep(1:5))
#
# t <- 1:5
# par(mfrow=c(1, 3))
# plot(bound_obf ~ t, type='l', ylim=c(0, 6))
# for(i in 1:dim(reps.1)[2]){
#   lines(reps.1[, i] ~ t, col='lightblue')
# }
# lines(bound_obf ~ t)
# plot(bound_obf ~ t, type='l', ylim=c(0, 6))
# for(i in 1:dim(reps.2)[2]){
#   lines(reps.2[, i] ~ t, col='lightblue')
# }
# lines(bound_obf ~ t)
# plot(bound_obf ~ t, type='l', ylim=c(0, 6))
# for(i in 1:dim(reps.3)[2]){
#   lines(reps.3[, i] ~ t, col='lightblue')
# }
# lines(bound_obf ~ t)
