
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

get.lower1 <- function(x) pmvnorm(lower=c(-Inf, -Inf), upper=c(-u_k, obs), mean=means(x), corr=var)
get.lower2 <- function(x) pmvnorm(lower=c(u_k, -Inf), upper=c(Inf, obs), mean=means(x), corr=var)

get.upper1 <- function(x) pmvnorm(lower=c(-Inf, obs), upper=c(-u_k, Inf), mean=means(x), corr=var)
get.upper2 <- function(x) pmvnorm(lower=c(u_k, obs), upper=c(Inf, Inf), mean=means(x), corr=var)

get.probs <- function(alt){
  lower1 <- pmvnorm(lower=c(-Inf, -Inf), upper=c(-u_k, obs), mean=means(alt), corr=var)
  lower2 <- pmvnorm(lower=c(u_k, -Inf), upper=c(Inf, obs), mean=means(alt), corr=var)

  upper1 <- pmvnorm(lower=c(-Inf, obs), upper=c(-u_k, Inf), mean=means(alt), corr=var)
  upper2 <- pmvnorm(lower=c(u_k, obs), upper=c(Inf, Inf), mean=means(alt), corr=var)

  lower <- max(lower1, 0) + max(lower2, 0)
  upper <- max(upper1, 0) + max(upper2, 0)

  return(c(lower, upper))
}

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
