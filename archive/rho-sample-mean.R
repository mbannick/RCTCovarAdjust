rm(list=ls())

setwd("simulations/")
source("../R/constants.R")
source("../R/pvalues.R")
source("sim-utils.R")

# OBF boundaries
u_k <- c(2.963, 1.969)
u_k1 <- cbind(-u_k, u_k)

# This gives exactly alpha = 0.05 by putting in the last
# boundary.
n_K <- c(50, 100)

# 1 ------------------------------------------------------
# SAMPLE MEAN

anova <- get.pvalue.sm(obs=-1.969/sqrt(100), u_K=t(matrix(u_k1[1:1,])), n_K=n_K,
              ancova_monitor=F, ancova_test=F)
get.pvalue.sm(obs=-1.969/sqrt(100), u_K=t(matrix(u_k1[1:1,])), n_K=n_K,
              ancova_monitor=F, ancova_test=T, rho=0.5)

vals.sm <- c()
i.s <- seq(0.01, 0.99, 0.01)

for(i in i.s){
    cat(i, ": ")
    res <- get.pvalue.sm(obs=-1.969/sqrt(100),
                        u_K=t(matrix(u_k1[1:1,])), n_K=n_K,
                        ancova_monitor=F, ancova_test=T, rho=i)
    cat("\n")
    vals.sm <- c(vals.sm, res)
}

# STAGE-WISE

vals.sw <- c()
anova.sw <- get.pvalue.sw(obs=-1.969,
                     u_k=t(matrix(u_k1[1:1,])), n_k=n_K,
                     ancova_monitor=F, ancova_test=F, rho=i,
                     last_stage=T)
for(i in i.s){
  cat(i, ": ")
  res <- get.pvalue.sw(obs=-1.969,
                       u_k=t(matrix(u_k1[1:1,])), n_k=n_K,
                       ancova_monitor=F, ancova_test=T, rho=i,
                       last_stage=T)
  cat("\n")
  vals.sw <- c(vals.sw, res)
}

load("empirical-ps-2.RData")


plot(vals.sm ~ i.s, xlab="rho", ylab="p-value", pch=1, col='red',
     ylim=c(min(vals.sw, vals.sm, vals[, 3]), max(vals.sw, vals.sm, vals[, 3])))
points(anova ~ 1, col='black', pch=1)
abline(h=0.05, lty='dashed')
points(anova.sw ~ 1, pch=4)
points(vals.sw ~ i.s, col='red', pch=4)
points(vals[, 1], vals[, 3], pch=8, col='blue')

legend(x=0.6, y=0.062, col=c('red', 'black', 'red', 'black', 'blue'),
       pch=c(1, 1, 4, 4, 8), legend=c("ANOVA, Sample Mean",
                                     "ANCOVA, Sample Mean",
                                     "ANOVA, Stage-Wise",
                                     "ANCOVA, Stage-Wise",
                                     "Simulation, ANCOVA, Stage-Wise"),
       cex=0.6)

# 3 --------------------------------------------------------------

# Observed values at the end of second stage (but that doesn't matter)
obs <- seq(-5.0, 5.0, by=0.1)
anova <- get.pvalue.sm(obs=obs/sqrt(20), u_K=t(matrix(u_k1[1:1,])), n_K=n_K,
                       ancova_monitor=F, ancova_test=F)

rhos <- c(0.25, 0.5, 0.75)
res.anova <- c()
for(o in obs){
  res <- get.pvalue.sm(obs=o/sqrt(20), u_K=t(matrix(u_k1[1:1,])), n_K=n_K,
                       ancova_monitor=F, ancova_test=F, rho=1)
  res.anova <- c(res.anova, res)
}
res.ancova <- matrix(data=NA, nrow=length(obs), ncol=length(rhos))
for(i in 1:length(obs)){
  for(j in 1:length(rhos)){
    res <- get.pvalue.sm(obs=obs[i]/sqrt(20), u_K=t(matrix(u_k1[1:1,])),
                         n_K=n_K,
                         ancova_monitor=F, ancova_test=T, rho=rhos[j])
    res.ancova[i, j] <- res
  }
}

cols <- c("red", "green", "blue")

plot(res.anova ~ obs, col='black', type='l')
for(i in 1:length(rhos)){
  lines(res.ancova[, i] ~ obs, col=cols[i])
}

