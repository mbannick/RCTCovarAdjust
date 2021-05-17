rm(list=ls())
source("R/trial-funcs.R")
source("R/ancova.R")
source("R/boundaries.R")
source("R/boundaries2.R")
source("R/trial-data.R")

sigma2 <- 1
gamma <- 0
theta <- 1
rho <- sqrt(sigma2) / sqrt(sigma2 + (gamma*theta)**2)
n_k <- c(500, 500)
t_k <- cumsum(n_k) / sum(n_k)
u_1 <- 2.943058

# OBF alpha spending function
a.func.obf <- function(a, t) 4 * (1 - pnorm(qnorm(1-a/4)/sqrt(t)))
a.spend <- a.func.obf(a=0.05, c(0, t_k)) %>% diff
spend.2 <- function(u, a, vec, u.1, vec.1) mean(abs(vec) >= u & abs(vec.1) < u.1) - a

dens <- get.joint.density(n_k=n_k, lr_bounds=c(u_1),
                          delta=0, rho=rho, gridsize=1e5)

corr.2 <- sqrt(t_k[1]) * rho
covv.2 <- matrix(c(1, corr.2, corr.2, 1), nrow=2, ncol=2)
set.seed(10)
t.2 <- mvrnorm(n=1000000, mu=c(0, 0), covv.2)
u.2.2 <- uniroot(function(u) spend.2(u, a.spend[2], t.2[, 2], u_1, t.2[, 1]), c(0, 100))$root

u.2.2

## THE MULTIVARIATE NORMAL HAS THE SAME MARGINAL DISTRIBUTION FOR

# these are the simulated test statistics based on sim
mean(abs(t.2[, 2]) >= u.2.2 & abs(t.2[, 1]) < u_1) + mean(abs(t.2[, 1]) >= u_1)

u.2.2 * sqrt(2)

# these are the test statistics density based on rec
dens <- get.joint.density(n_k=n_k, lr_bounds=c(u_1),
                          gridsize=1e6, rho=rho)

sum(dens$dens[, 2][abs(dens$grid) >= u.2.2*sqrt(2)]) * 2e-5

tau <- get.tau(n_k)
tau_s <- cumsum(tau)
grid2 <- dens$grid / sqrt(tau_s[2])

A <- sqrt(tau_s[2])
D <- dens$dens[, 2]
X <- dens$grid

Y <- X / A
cumsum(D)
plot(cumsum(D * 2e-4)~X, type='l')

# DOING THE CONVERSION BETWEEN THE SCORE AND LIKELIHOOD RATIO TOTALLY WRONG
# AND FOR SOME REASON THAT'S MESSING UP THE DENSITY BECAUSE IT SHOULD BE A LOT
# LESS FAT.

hist(t.2[, 2], freq=FALSE, breaks=100, prob=TRUE, add=F)
new.sim <- t.2[, 2]*sqrt(tau_s[2])
hist(new.sim, freq=FALSE, breaks=100, prob=TRUE, add=T, col='blue')
lines(dens$dens[, 2] ~ dens$grid) # WRONG - SCORE SCALE
tot <- sum(dens$dens[, 2] * grid2)
lines(dens$dens[, 2] / tot ~ grid2, col='blue')
lines(density(t.2[, 2]), col='red')
samp <- replicate(10, sample(x=dens$grid, prob=dens$dens[, 2], replace=T))
hist(samp, freq=FALSE, breaks=100, prob=TRUE, add=T, col='green')


h <- hist(t.2[, 2], breaks=100, plot=F)


x <- seq(-5, 5, by=0.01)
y <- dnorm(x)
py <- pnorm(x)
plot(py~x, type='l')
lines(y~x, type='l')


samp <- replicate(10, sample(x=dens$grid, prob=dens$dens[, 2], replace=T)) %>% c
# hist(samp)
mu <- mean(samp)
sd <- var(samp)**0.5

new.samp <- samp/sqrt(tau_s[2])
hist(new.samp, add=T, col='pink', breaks=100, freq=FALSE, prob=TRUE, alpha=0.4)

lines(dnorm(x, mean=mu, sd=sd) ~ x, col='orange')
