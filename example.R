source("R/trial-funcs.R")

t_k <- 1:5/5
p_k <- (1:5/5)**2

a.func.obf <- function(a, t) 4 * (1 - pnorm(qnorm(1-a/4)/sqrt(t)))
fit.stage(df, stage=5, bounds=c(4.562, 3.226, 2.634, 2.281), ancova=F,
          a.func=a.func.obf, a=0.05,
          rates=t_k, N=1000, n_sims=10000)

set.seed(10)
reps.1 <- replicate(2, run.full.trial(estimate_sigma=FALSE,
                                    do_ancova=FALSE, theta=c(1, 1),
                                    gamma=c(1, 1)))

set.seed(10)
reps.2 <- replicate(2, run.full.trial(estimate_sigma=FALSE,
                                    do_ancova=TRUE, theta=c(1, 1),
                                    gamma=c(10, 1)))

set.seed(10)
reps.3 <- replicate(2, run.full.trial(estimate_sigma=TRUE,
                                    do_ancova=TRUE, theta=c(1, 1),
                                    gamma=c(10, 1)))

bound_obf <- c(1.9600, 2.7965, 3.4711, 4.0486, 4.5617,
               5.0283, 5.5490, 5.8611, 6.2395, 6.5981,
               6.9396, 7.2663, 7.5799, 7.8820, 8.1736)
bound_obf <- bound_obf[5] / sqrt(rep(1:5))

t <- 1:5
par(mfrow=c(1, 3))
plot(bound_obf ~ t, type='l', ylim=c(0, 6))
for(i in 1:dim(reps.1)[2]){
  lines(reps.1[, i] ~ t, col='lightblue')
}
lines(bound_obf ~ t)
plot(bound_obf ~ t, type='l', ylim=c(0, 6))
for(i in 1:dim(reps.2)[2]){
  lines(reps.2[, i] ~ t, col='lightblue')
}
lines(bound_obf ~ t)
plot(bound_obf ~ t, type='l', ylim=c(0, 6))
for(i in 1:dim(reps.3)[2]){
  lines(reps.3[, i] ~ t, col='lightblue')
}
lines(bound_obf ~ t)
