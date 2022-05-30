
set.seed(100)
sim.data <- sim.data.closure(delta=0.5, n_cov=1, rho=0.9)
data <- sim.data(100)

anova <- fit.model.closure(ancova=FALSE)
ancova <- fit.model.closure(ancova=TRUE)

anova.res <- anova(data$X, data$y)
ancova.res <- ancova(data$X, data$y)

u_k <- 3
rho <- sqrt(ancova.res$variance / anova.res$variance)

# We have rejected at the first stage based on ANOVA and want to calculate
# a p-value, CI, and point estimate based on ANCOVA

# The distribution of the test statistics under the NULL hypothesis is a bivariate normal distribution
# with covariance matrix ((1, 0.9), (0.9, 1))

# The distribution of the test statistics under the ALTERNATIVE hypothesis with mean d is a
# bivariate normal distribution with means alt * sqrt(n) / \sigma

# First do just ANOVA p-value and confidence interval
obs <- anova.res$tstat
means <- function(alt) alt * 10 / anova.res$variance**0.5
lower <- pnorm(obs, lower.tail=TRUE)
upper <- pnorm(obs, lower.tail=FALSE)
pval <- 2*min(c(lower, upper))

ci.lower <- uniroot(function(x) 0.025 - pnorm(obs, mean=means(x), lower.tail=F), interval=c(-100, 100))$root
ci.upper <- uniroot(function(x) 0.025 - pnorm(obs, mean=means(x), lower.tail=T), interval=c(-100, 100))$root
med.esti <- uniroot(function(x) 0.500 - pnorm(obs, mean=means(x), lower.tail=F), interval=c(-100, 100))$root

c(ci.lower, med.esti, ci.upper)

alts <- seq(-3, 3, by=0.01)
vals.lower <- sapply(alts, function(x) pnorm(obs, mean=means(x), lower.tail=F))
vals.upper <- sapply(alts, function(x) pnorm(obs, mean=means(x), lower.tail=T))

plot(vals.lower ~ alts, type='l')
lines(vals.upper ~ alts, col='blue')

# Now do ANOVA --> ANCOVA p-value
obs <- ancova.res$tstat
means <- function(alt) alt * 10 / c(anova.res$variance, ancova.res$variance)**0.5

var <- matrix(c(1, rho, rho, 1), nrow=2, byrow=FALSE)

# get.lower1 <- function(x) pmvnorm(lower=c(-Inf, -Inf), upper=c(-u_k, obs), mean=means(x), corr=var)
# get.lower2 <- function(x) pmvnorm(lower=c(u_k, -Inf), upper=c(Inf, obs), mean=means(x), corr=var)
#
# get.upper1 <- function(x) pmvnorm(lower=c(-Inf, obs), upper=c(-u_k, Inf), mean=means(x), corr=var)
# get.upper2 <- function(x) pmvnorm(lower=c(u_k, obs), upper=c(Inf, Inf), mean=means(x), corr=var)

# THIS WAS FINE FOR SINGLE STAGE

# alt <- 0.0
# u_k <- matrix(c(-3, 3), nrow=1)
#
# # Probability of not stopping
# p.notstop <- pnorm(u_k[1, 2], mean=alt, lower.tail=T) -
#   pnorm(u_k[1, 1], mean=alt, lower.tail=T)
#
# # Stopping in this stage and more negative
# p.L.neg <- pnorm(min(obs, u_k[1, 1]), mean=alt, lower.tail=T)
# p.L.pos <- pnorm(max(obs, u_k[1, 2]), mean=alt, lower.tail=T) -
#   pnorm(u_k[1, 2], mean=alt, lower.tail=T)
#
# p.lower <- p.L.neg + p.L.pos + p.notstop
#
# # Stop in this stage, and more positive
# p.U.neg <- pnorm(min(obs, u_k[1, 1]), mean=alt, lower.tail=F) -
#   pnorm(u_k[1, 1], mean=alt, lower.tail=F)
# p.U.pos <- pnorm(max(obs, u_k[1, 2]), mean=alt, lower.tail=F)
#
# p.upper <- p.U.neg + p.U.pos

l.K <- u_k[K, 1] # lower bound at this stage
u.K <- u_k[K, 2] # upper bound at this stage

# Stop in this stage with MONITOR and more negative TEST, or do not stop
lower.cross <- c(-Inf, l.K)
upper.cross <- c(u.K, Inf)
no.cross <- c(l.K, u.K)

# More negative TEST
lower.tail <- c(-Inf, obs)
upper.tail <- c(obs, Inf)

lower.block.1 <- rbind(prev, lower.cross, lower.tail)
upper.block.1 <- rbind(prev, upper.cross, upper.tail)

lower.block.2 <- rbind(prev, upper.cross, lower.tail)
upper.block.2 <- rbind(prev, lower.cross, upper.tail)

lower.block.3 <- rbind(prev, no.cross, c(-Inf, Inf))

lower.blocks <- list(lower.block.1, lower.block.2, lower.block.3)
upper.blocks <- list(upper.block.1, upper.block.2)

p.lower <- .pmvnorm.list(blocks=lower.blocks, corr=corr, mean=alt)
p.upper <- .pmvnorm.list(blocks=upper.blocks, corr=corr, mean=alt)

get.probs <- function(alt){
  # lower1 <- pmvnorm(lower=c(-Inf, -Inf), upper=c(-u_k, obs), mean=means(alt), corr=var)
  # lower2 <- pmvnorm(lower=c(u_k, -Inf), upper=c(Inf, obs), mean=means(alt), corr=var)
  #
  # upper1 <- pmvnorm(lower=c(-Inf, obs), upper=c(-u_k, Inf), mean=means(alt), corr=var)
  # upper2 <- pmvnorm(lower=c(u_k, obs), upper=c(Inf, Inf), mean=means(alt), corr=var)
  #
  # lower <- max(lower1, 0) + max(lower2, 0)
  # upper <- max(upper1, 0) + max(upper2, 0)

  p.notstop <- pnorm(u_k[1, 2], mean=means(alt), lower.tail=T) -
    pnorm(u_k[1, 1], mean=means(alt), lower.tail=T)

  # Stopping in this stage and more negative
  p.L.neg <- pnorm(min(obs, u_k[1, 1]), mean=means(alt), lower.tail=T)
  p.L.pos <- pnorm(max(obs, u_k[1, 2]), mean=means(alt), lower.tail=T) -
    pnorm(u_k[1, 2], mean=means(alt), lower.tail=T)

  lower <- p.L.neg + p.L.pos + p.notstop

  # Stop in this stage, and more positive
  p.U.neg <- pnorm(min(obs, u_k[1, 1]), mean=means(alt), lower.tail=F) -
    pnorm(u_k[1, 1], mean=means(alt), lower.tail=F)
  p.U.pos <- pnorm(max(obs, u_k[1, 2]), mean=means(alt), lower.tail=F)

  upper <- p.U.neg + p.U.pos

  return(c(lower, upper))
}

alts <- seq(-3, 3, by=0.01)
vals <- t(sapply(alts, get.probs))

plot(vals[, 1] ~ alts, type='l')
lines(vals[, 2] ~ alts, col='blue')


lower1s <- sapply(alts, get.lower1)
lower2s <- sapply(alts, get.lower2)
upper1s <- sapply(alts, get.upper1)
upper2s <- sapply(alts, get.upper2)

par(mfrow=c(2, 2))
plot(lower1s ~ alts)
abline(v=ancova.res$delta)
plot(lower2s ~ alts)
abline(v=ancova.res$delta)
plot(upper1s ~ alts)
abline(v=ancova.res$delta)
plot(upper2s ~ alts)
abline(v=ancova.res$delta)

library(mnormt)

x <- seq(-10, 10, 0.25)
y <- seq(-10, 10, 0.25)

build.contour <- function(alt){
  f <- function(x, y) dmnorm(cbind(x, y), mean=means(alt), varcov=var)
  z <- outer(x, y, f)
  levels <- c(0.00001, 0.0001, 0.001, 0.01, 0.1)
  contour(x, y, z, levels=levels)
}

add.lines <- function(){
  abline(h=0, lty=2)
  abline(v=0, lty=2)
  abline(v=u_k, lty=2, col='red')
  abline(v=-u_k, lty=2, col='red')
  abline(h=obs, col='blue')
  rect(u_k, -10, 10, obs, density = 5, border = "red", lty =2)
}

par(mfrow=c(2, 4))
for(alt in seq(0, 0.7, by=0.1)){
  build.contour(alt)
  add.lines()
}

# build.contour(0)
# add.lines()
#
# build.contour(0.2)
# add.lines()
#
#
#
#
# xs <- seq(-3, 3, by=0.01)
# sd_K <- 0.9755748
# est <- 0.1970725
# n_k <- 330
# rho <- 1
# sd_k <- c(sd_K/rho, sd_K)
# get.z <- function(eff) sqrt(n_k) * (eff) / sd_k
# ps.L <- sapply(xs, function(x) get.pvalue.sw(3.66, alt=get.z(x),
#                                              u_k=matrix(c(-2.282511, 2.282511), nrow=1),
#                                              n_k=c(330), k_r=1, rho=rho,
#                                              ancova_monitor=F,
#                                              ancova_test=F,
#                                              last_stage=F, type="lower"))
# ps.U <- sapply(xs, function(x) get.pvalue.sw(3.66, alt=get.z(x),
#                                              u_k=matrix(c(-2.282511, 2.282511), nrow=1),
#                                              n_k=c(330), k_r=1, rho=rho,
#                                              ancova_monitor=F,
#                                              ancova_test=F,
#                                              last_stage=F, type="upper"))
# plot(ps.L ~ xs, type='l')
# lines(ps.U ~ xs, col='blue')
# abline(v=est, lty=2)
#
#
# # probability of rejecting in the first stage
# pnorm()
#
#
#
#
#
#
# get.pvalue.sw(3.66, alt=get.z(0.19),
#               u_k=matrix(c(-2.282511, 2.282511), nrow=1),
#               n_k=c(330), k_r=1, rho=0.99,
#               ancova_monitor=F,
#               ancova_test=T,
#               last_stage=F, type="lower")

# u_k2 <- matrix(c(-2.797, 2.797, -1.977, 1.977), nrow=2, byrow=T)





# STOP AT THE SECOND STAGE IN A 3-STAGE TRIAL, ANOVA-ANOVA --
# LOOKS WEIRD WHEN HAVE NEGATIVE OBSERVED Z STATISTICS

obs <- 4
var <- 1
n_k <- c(10, 20)
means <- function(alt) alt * sqrt(n_k) / var**0.5

u_k3 <- matrix(c(-3.471, 3.471, -2.454, 2.454, -2.004, 2.004), nrow=3, byrow=T)

K <- 2
corr <- matrix(c(1, sqrt(1/2), sqrt(1/2), 1), byrow=T, nrow=2)

prev <- u_k3[1:(K-1), ]

get.probs <- function(alt){
  # Stop at previous stage
  p.L.prev <- pnorm(u_k3[1, 1], mean=means(alt)[1], lower.tail=T)
  p.U.prev <- pnorm(u_k3[1, 2], mean=means(alt)[1], lower.tail=F)

  # Stop in this stage and more negative, or do not stop
  p.notstop <- rbind(prev, u_k3[K,])
  p.L.neg <- rbind(prev, c(-Inf, min(obs, u_k3[K, 1])))
  p.L.pos <- rbind(prev, c(u_k3[K, 2], max(obs, u_k3[K, 2])))

  # Stop in this stage, and more positive
  p.U.neg <- rbind(prev, c(min(u_k3[K, 1], obs), u_k3[K, 1]))
  p.U.pos <- rbind(prev, c(max(obs, u_k3[K, 2]), Inf))

  lower.blocks <- list(p.L.neg, p.L.pos, p.notstop)
  upper.blocks <- list(p.U.neg, p.U.pos)# p.notstop)

  p.lower <- .pmvnorm.list(blocks=lower.blocks, corr=corr, mean=means(alt))
  p.upper <- .pmvnorm.list(blocks=upper.blocks, corr=corr, mean=means(alt))

  lower <- p.lower + p.L.prev
  upper <- p.upper + p.U.prev
  return(c(lower, upper))
}

alts <- seq(-3, 3, by=0.01)
vals <- t(sapply(alts, get.probs))

plot(vals[, 1] ~ alts, type='l')
lines(vals[, 2] ~ alts, col='blue')

# STOP AT THE FIRST STAGE IN A 3-STAGE TRIAL, ANOVA-ANCOVA

obs <- -4
var <- c(1, 0.7)
n_k <- c(10, 10)
means <- function(alt) alt * sqrt(n_k) / var**0.5

u_k3 <- matrix(c(-3.471, 3.471, -2.454, 2.454, -2.004, 2.004), nrow=3, byrow=T)

K <- 1
corr <- matrix(c(1, sqrt(0.7), sqrt(0.7), 1), byrow=T, nrow=2)

prev <- NULL

get.probs <- function(alt){
  l.K <- u_k3[K, 1] # lower bound at this stage
  u.K <- u_k3[K, 2] # upper bound at this stage

  # Stop in this stage with MONITOR and more negative TEST, or do not stop
  lower.cross <- c(-Inf, l.K)
  upper.cross <- c(u.K, Inf)
  no.cross <- c(l.K, u.K)

  # More negative TEST
  lower.tail <- c(-Inf, obs)
  upper.tail <- c(obs, Inf)

  lower.block.1 <- rbind(prev, lower.cross, lower.tail)
  upper.block.1 <- rbind(prev, upper.cross, upper.tail)

  lower.block.2 <- rbind(prev, upper.cross, lower.tail)
  upper.block.2 <- rbind(prev, lower.cross, upper.tail)

  lower.block.3 <- rbind(prev, no.cross, c(-Inf, Inf))

  lower.blocks <- list(lower.block.1, lower.block.2, lower.block.3)
  upper.blocks <- list(upper.block.1, upper.block.2)

  p.lower <- .pmvnorm.list(blocks=lower.blocks, corr=corr, mean=means(alt))
  p.upper <- .pmvnorm.list(blocks=upper.blocks, corr=corr, mean=means(alt))
  return(c(p.lower, p.upper))
}

alts <- seq(-3, 3, by=0.01)
vals <- t(sapply(alts, get.probs))

plot(vals[, 1] ~ alts, type='l')
lines(vals[, 2] ~ alts, col='blue')







# NEW DEBUGGING AS OF APRIL 5 EVENING

obs <- 2.227904
var <- c(0.9442099, 0.9442099, 0.9393716)
rho <- sqrt(0.9393716 / 0.9442099)
n_k <- c(16.17, 32.83, 32.83)
means <- function(alt) alt * sqrt(n_k) / var**0.5

u_k3 <- matrix(c(-2.282511, 2.282511, -2.291232, 2.291232), nrow=2, byrow=T)

K <- 2
corr <- matrix(c(1, sqrt(16.17/32.83), sqrt(16.17/32.83)*rho,
                 sqrt(16.17/32.83), 1, rho,
                 sqrt(16.17/32.83)*rho, rho, 1), byrow=T, nrow=3)

prev <- NULL

get.probs <- function(alt){
  l.K <- u_k3[K, 1] # lower bound at this stage
  u.K <- u_k3[K, 2] # upper bound at this stage

  # Stop in this stage with MONITOR and more negative TEST, or do not stop
  lower.cross <- c(-Inf, l.K)
  upper.cross <- c(u.K, Inf)
  no.cross <- c(l.K, u.K)

  # More negative TEST
  lower.tail <- c(-Inf, obs)
  upper.tail <- c(obs, Inf)

  lower.block.1 <- rbind(prev, lower.cross, lower.tail)
  upper.block.1 <- rbind(prev, upper.cross, upper.tail)

  lower.block.2 <- rbind(prev, upper.cross, lower.tail)
  upper.block.2 <- rbind(prev, lower.cross, upper.tail)

  nocross.block <- rbind(prev, no.cross, c(-Inf, Inf))

  lower.blocks <- list(lower.block.1, lower.block.2)
  upper.blocks <- list(upper.block.1, upper.block.2)

  if(obs <= l.K){
    lower.blocks <- list(lower.block.1, lower.block.2)
    upper.blocks <- list(upper.block.1, upper.block.2, nocross.block)
  } else if(obs >= u.K) {
    lower.blocks <- list(lower.block.1, lower.block.2, nocross.block)
    upper.blocks <- list(upper.block.1, upper.block.2)
  } else {
    lower.blocks <- list(lower.block.1, lower.block.2)
    upper.blocks <- list(upper.block.1, upper.block.2)
  }

  p.lower <- .pmvnorm.list(blocks=lower.blocks, corr=corr, mean=means(alt))
  p.upper <- .pmvnorm.list(blocks=upper.blocks, corr=corr, mean=means(alt))
  return(c(p.lower, p.upper))
}

alts <- seq(-3, 3, by=0.01)
vals <- t(sapply(alts, get.probs))

plot(vals[, 1] ~ alts, type='l')
lines(vals[, 2] ~ alts, col='blue')
