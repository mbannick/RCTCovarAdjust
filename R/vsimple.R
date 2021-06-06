rm(list=ls())
source("R/trial-funcs.R")
source("R/ancova.R")
source("R/boundaries.R")
source("R/boundaries2.R")
source("R/boundaries3.R")
source("R/trial-data.R")
library(ggplot2)
library(RColorBrewer)

try.trial <- function(K, gamma, new=TRUE){

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
                             estimate_sigma=TRUE, do_ancova=TRUE, new=FALSE) {

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
                          rates=t_k[1:i], N=N, n_sims=n_sims, new=new)
    }

    return(bounds)
  }

  set.seed(10)
  trial <- run.full.trial(gamma=c(gamma), theta=c(1),
                          estimate_sigma=FALSE,
                          do_ancova=FALSE,
                          N=1000,
                          n_sims=100000,
                          a.func=a.func.obf, new=new)
  trial.2 <- run.full.trial(gamma=c(gamma), theta=c(1),
                          estimate_sigma=FALSE,
                          do_ancova=TRUE,
                          N=1000,
                          n_sims=100000,
                          a.func=a.func.obf, new=new)
  print(bound_obf)
  print(trial)
  return(cbind(trial, trial.2))
}

t1 <- try.trial(2, gamma=3)
t2 <- try.trial(2, gamma=3, new=TRUE)

TRIALS.1 <- sapply(2:5, try.trial, gamma=0.1)
TRIALS.2 <- sapply(2:5, try.trial, gamma=0.5)
TRIALS.3 <- sapply(2:5, try.trial, gamma=1)
TRIALS.4 <- sapply(2:5, try.trial, gamma=1.5)

base <- TRIALS.1
b <- c(base[[1]][2, 1], base[[2]][3, 1], base[[3]][4, 1], base[[4]][5, 1])
gammas <- c(0.0, 0.1, 0.5, 1, 1.5)
times <- c(2, 3, 4, 5)

data <- data.table(b=b, times=times, gamma=c(0))
for(i in 1:4){
  b.trial <- get(paste0("TRIALS.", i))
  bs <- c(b.trial[[1]][2, 2],
          b.trial[[2]][3, 2],
          b.trial[[3]][4, 2],
          b.trial[[4]][5, 2])
  df <- data.table(b=bs, times=times, gamma=gammas[i+1])
  data <- rbind(data, df)
}

colors <- data.table(color=colors, gamma=gammas)
data <- merge(data, colors, on='gamma')

data[, stop := times + 0.5]
data[, start := times - 0.5]

cols <- brewer.pal(10, "Spectral")[1:4]
cols <- c("black", cols)

pdf("critical-values-2.pdf", height=5, width=8)
ggplot(data=data) +
  geom_segment(aes(x=start, y=b, xend=stop, yend=b,
                   color=factor(gamma)), size=1) +
  theme_minimal() + xlab("Maximum Num. of Stages (K)") +
  ylab("Critical Value") +
  scale_color_manual(values=cols) +
  labs(color="Covariate \nCoefficient")
dev.off()


brewer.pal(5, "Set1")
