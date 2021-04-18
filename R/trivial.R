rm(list=ls())

library(MASS)
library(magrittr)

# Information fraction
# and data generating process parameters
n_k <- c(100, 100)
t_k <- cumsum(n_k) / sum(n_k)
gamma <- 6; theta <- 1; sigma2 <- 1
nsim <- 10000

# OBF alpha spending function
a.func.obf <- function(a, t) 4 * (1 - pnorm(qnorm(1-a/4)/sqrt(t)))
a.spend <- a.func.obf(a=0.05, c(0, t_k)) %>% diff

# Get correlations for 1 v. 2 based on covariate adjustment or not
corr.1 <- sqrt(t_k[1])
corr.2 <- sqrt(t_k[1]) * (sigma2 / (sigma2 + (gamma * theta)**2))

# Create covariance matrices
covv.1 <- matrix(c(1, corr.1, corr.1, 1), nrow=2, ncol=2)
covv.2 <- matrix(c(1, corr.2, corr.2, 1), nrow=2, ncol=2)

set.seed(10)
t.1 <- mvrnorm(n=nsim, mu=c(0, 0), covv.1)
set.seed(10)
t.1.1 <- mvrnorm(n=nsim, mu=c(0), matrix(1))
set.seed(10)
t.2 <- mvrnorm(n=nsim, mu=c(0, 0), covv.2)

spend.1 <- function(u, a, vec) mean(abs(vec) >= u) - a
spend.2 <- function(u, a, vec, u.1, vec.1) mean(abs(vec) >= u & abs(vec.1) < u.1) - a

u.1 <- uniroot(function(u) spend.1(u, a.spend[1], t.1.1[, 1]), c(0, 100))$root

# I HAVE CHECKED THAT U.2.1 MATCHES WITH THE EXAMPLE
u.2.1 <- uniroot(function(u) spend.2(u, a.spend[2], t.1[, 2], u.1, t.1[, 1]), c(0, 100))$root
# THIS ONE DOES NOT MATCH WITH THE EXAMPLE; INSTEAD IT IS IDENTICAL TO ANOVA
u.2.2 <- uniroot(function(u) spend.2(u, a.spend[2], t.2[, 2], u.1, t.2[, 1]), c(0, 100))$root

# CORRECTION:
# YOU CAN ONLY SIMULATE WHAT YOU KNOW ABOUT AT THIS TIME.
# IF THAT MEANS ONLY ONE STAGE, THEN ONLY SIMULATE THAT ONE STAGE!
# FIX THE RUN.FULL.TRIAL FUNCTION
