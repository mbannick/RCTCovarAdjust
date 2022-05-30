library(data.table)
library(ggplot2)

obs <- -1.35
u_k <- matrix(c(-3.730665, 3.730665,
              -2.503871, 2.503871), nrow=2, byrow=T)
n_k <- c(16.17, 32.83)
rho <- 0.11
type <- "two-sided"

get.pvalue.sw(obs=1.96, u_k=matrix(c(-1.96, 1.96), nrow=1),
              n_k=c(100), k_r=1, rho=0.11,
              ancova_monitor=F, ancova_test=T, last_stage=F,
              type="two-sided")

get.pvalue.sw(obs=obs, u_k=u_k, n_k=n_k, k_r=2, rho=rho,
              ancova_monitor=F, ancova_test=T, last_stage=F,
              type="two-sided")

get.pvalue.sw(obs=obs, u_k=u_k, n_k=n_k, k_r=2, rho=rho,
              ancova_monitor=F, ancova_test=F, last_stage=F,
              type="two-sided")

get.confint.sw(est=obs, sd_K=1, u_k=u_k, n_k=n_k, k_r=2, rho=rho,
              ancova_monitor=F, ancova_test=F, last_stage=F, alpha=0.05)

# I think we're missing a probability to add on to this...
# The p-value above should be LARGER when we do two tests, not smaller...

obs <- 1.0
rho <- 0.8
bound <- 1.96

obs <-
get.pval <- function(x, rho) get.pvalue.sw(obs=x,
                                      u_k=matrix(c(-1.96, 1.96), nrow=1),
                                      n_k=c(100), k_r=1, rho=rho,
                                      ancova_monitor=F, ancova_test=T,
                                      last_stage=F,
                                      type="two-sided")

xs <- seq(-3, 3, by=0.5)
rhos <- seq(0, 0.9, by=0.1)
grid <- expand.grid(x=xs, rho=rhos)

res <- mapply(FUN=get.pval, x=grid$x, rho=grid$rho)
grid$pval <- res

grid <- data.table(grid)
ggplot(data=grid, aes(x=x, y=rho, fill=pval)) + geom_tile()

pcnorm <- function(x, z, rho) pnorm(x, mean=rho*z, (1-rho^2))

p.ANCOVA <- function(obs, rho, bound){

  corr <- matrix(c(1, rho, rho, 1), byrow=T, nrow=2)

  lower.UU <- c(bound, obs)
  upper.UU <- c(Inf, Inf)

  lower.UL <- c(bound, -Inf)
  upper.UL <- c(Inf, obs)

  lower.LL <- c(-Inf, -Inf)
  upper.LL <- c(-bound, obs)

  lower.LU <- c(-Inf, obs)
  upper.LU <- c(-bound, Inf)

  p.UU <- pmvnorm(lower=lower.UU, upper=upper.UU, corr=corr, algorithm=GenzBretz(abseps=1e-8))
  p.UL <- pmvnorm(lower=lower.UL, upper=upper.UL, corr=corr, algorithm=GenzBretz(abseps=1e-8))
  p.LL <- pmvnorm(lower=lower.LL, upper=upper.LL, corr=corr, algorithm=GenzBretz(abseps=1e-8))
  p.LU <- pmvnorm(lower=lower.LU, upper=upper.LU, corr=corr, algorithm=GenzBretz(abseps=1e-8))

  pval <- 2*min(c(p.UU + p.LU), c(p.LL + p.UL))

  return(pval)
}


p.ANCOVA <- function(rho) ((1-pcnorm(obs, z=bound, rho=rho)) * 2 +
            (1-pcnorm(obs, z=-bound, rho=rho)) * 2) * (1-pnorm(bound)) * 2

p.ANCOVA2 <- function(rho){
  corr <- matrix(c(1, rho, rho, 1), byrow=T, nrow=2)

  lower <- c(-Inf, -Inf)
  upper <- c(-bound, obs)
}

rhos <- seq(0, 1, by=0.05)
ps <- sapply(rhos, p.ANCOVA)
plot(ps ~ rhos)

# What would the p-value be for ANOVA?
p.ANOVA <- (1-pnorm(1.96)) * 2

p.ANCOVA <- 1-pmvnorm(lower=bounds[, 1], upper=bounds[, 2],
                    ..., algorithm=GenzBretz(abseps=1e-8))


u_k <- matrix(c(-3.730665, 3.730665,
                -2.503871, 2.503871), nrow=2, byrow=T)
n_k <- c(16.17, 32.83)
rho <- 0.11
type <- "two-sided"
