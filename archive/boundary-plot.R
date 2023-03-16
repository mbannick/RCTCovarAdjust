rm(list=ls())

source("R/ancova.R")
source("R/boundaries.R")
source("R/covariance.R")

library(ggplot2)
library(RColorBrewer)
library(mvtnorm)

BOUND_OBF <- c(1.9600, 2.7965, 3.4711, 4.0486, 4.5617,
               5.0283, 5.5490, 5.8611, 6.2395, 6.5981,
               6.9396, 7.2663, 7.5799, 7.8820, 8.1736)

make.boundaries <- function(K, gamma){

  t_k <- 1:K/K
  p_k <- t_k

  bound_obf <- BOUND_OBF[K] / sqrt(rep(1:K))
  a.func <- OBF.SPEND(0.05)
  rho <- 1/sqrt(1 + gamma**2)

  bounds1 <- get.boundaries(a.func=a.func,
                            rates=t_k, rho=1)
  bounds2 <- get.boundaries(a.func=a.func,
                            rates=t_k, rho=rho)

  return(cbind(bounds1, bounds2))
}

t1 <- make.boundaries(2, gamma=3)

TRIALS.1 <- sapply(2:5, make.boundaries, gamma=0.1)
TRIALS.2 <- sapply(2:5, make.boundaries, gamma=0.5)
TRIALS.3 <- sapply(2:5, make.boundaries, gamma=1)
TRIALS.4 <- sapply(2:5, make.boundaries, gamma=1.5)

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

pdf("critical-values-3.pdf", height=5, width=8)
ggplot(data=data) +
  geom_segment(aes(x=start, y=b, xend=stop, yend=b,
                   color=factor(gamma)), size=1) +
  theme_minimal() + xlab("Maximum Num. of Stages (K)") +
  ylab("Critical Value") +
  scale_color_manual(values=cols) +
  labs(color="Covariate \nCoefficient")
dev.off()
