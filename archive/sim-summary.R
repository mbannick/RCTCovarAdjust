rm(list=ls())

in_dir <- "~/Documents/FileZilla/rct/"
setwd(in_dir)

stds <- c(0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.115, 0.2, 0.3, 0.75)
tasks <- 80
N <- 800000

for(st in stds){
  cat("\n", st)
  var_stage <- sprintf("stage_%s", st)
  var_tstat <- sprintf("tstat_%s", st)
  var_smean <- sprintf("smean_%s", st)

  assign(var_stage, c())
  assign(var_tstat, c())
  assign(var_smean, c())

  for(i in 1:80){
    cat(".")
    load(sprintf("task_%s_std_%s.RData", i, st))
    early <- result2$l.bounds[, 2]
    stage <- rep(1, length(early))
    stage[early] <- 2

    tstat <- result2$tstat
    smean <- result2$smean

    assign(var_stage, c(get(var_stage), stage))
    assign(var_tstat, c(get(var_tstat), tstat))
    assign(var_smean, c(get(var_smean), smean))
  }
}

# This is how to get rho from the observation standard deviations
# sim_stds <- c(0.001, 0.115)
rho <- function(x) sqrt(x**2 / (0.1**2 + x**2))
sim_rho <- rho(stds)

# holders for p-value estimates
sim_pvals_smean <- c()
sim_pvals_tstat <- c()

# standard errors of p-values
sim_pvals_smean_se <- c()
sim_pvals_tstat_se <- c()

# observed values on the boundaries
obs_smean <- -1.969/sqrt(100)
obs_tstat <- -1.969

# loop through each of the observation standard deviations
for(st in stds){

  std_name <- sprintf("%03d", st*1000)

  # CALCULATE SAMPLE MEAN P-VALUE EST AND P-VALUE SE
  smean_sims <- get(sprintf("smean_%s", st))
  # Take fraction of time that simulation is more extreme than observed (two-sided)
  pval_smean <- mean((smean_sims <= -abs(obs_smean)) | (smean_sims >= abs(obs_smean)))
  pval_smean_se <- sqrt(pval_smean * (1 - pval_smean) / N)

  # CALCULATE STAGE-WISE MEAN P-VALUE EST AND P-VALUE SE
  tstat_sims <- get(sprintf("tstat_%s", st))
  stage_sims <- get(sprintf("stage_%s", st))
  # probability of rejecting early
  reject_early <- stage_sims == 1
  last_stage_vals <- tstat_sims[!reject_early]
  # p-value is prob. of rejecting early + more extreme in second stage
  pval_tstat <- mean(!stage_sims) +
    mean((last_stage_vals <= -abs(obs_tstat)) | (last_stage_vals >= abs(obs_tstat)))
  pval_tstat_se <- sqrt(pval_tstat * (1 - pval_tstat) / N)

  # save estimates
  sim_pvals_smean <- c(sim_pvals_smean, pval_smean)
  sim_pvals_tstat <- c(sim_pvals_tstat, pval_tstat)

  # save standard errors
  sim_pvals_smean_se <- c(sim_pvals_smean_se, pval_smean_se)
  sim_pvals_tstat_se <- c(sim_pvals_tstat_se, pval_tstat_se)
}

setwd("~/repos/RCTCovarAdjust/simulations/")
source("../R/constants.R")
source("../R/pvalues.R")
source("sim-utils.R")

# OBF boundaries
u_k <- c(2.963, 1.969)
u_k1 <- cbind(-u_k, u_k)

# This gives exactly alpha = 0.05 by putting in the last
# boundary.
n_K <- c(50, 100)

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

# PLOT
par(mfrow=c(1, 2))
# SAMPLE MEAN
plot(vals.sm ~ i.s, xlab="rho", ylab="p-value", type='l', col='red',
     ylim=c(0.04, 0.058),
     main="Sample Mean")
points(anova ~ 1, col='black', pch=16)
abline(h=0.05, lty='dashed')

# Simulations
points(sim_pvals_smean ~ sim_rho, pch=16, col='blue')
lower <- sim_pvals_smean - 1.96 * sim_pvals_smean_se
upper <- sim_pvals_smean + 1.96 * sim_pvals_smean_se
arrows(x0=sim_rho, y0=lower, x1=sim_rho, y1=upper, code=3, angle=90, length=0.05, col='blue')

legend(x=0.6, y=0.062, col=c('red', 'black', 'blue'),
       pch=c(NA, 16, 16), lty=c(1, NA, 1), legend=c("ANOVA",
                                                    "ANCOVA",
                                                    "Simulation"), cex=0.7)

# STAGEWISE
plot(vals.sw ~ i.s, xlab="rho", ylab="p-value", type='l', col='red',
     ylim=c(0.04, 0.058),
     main="Stage Wise")
points(anova.sw ~ 1, pch=16, col='black')
abline(h=0.05, lty='dashed')
points(sim_pvals_tstat ~ sim_rho, pch=16, col='blue')

# Simulations
lower <- sim_pvals_tstat - 1.96 * sim_pvals_tstat_se
upper <- sim_pvals_tstat + 1.96 * sim_pvals_tstat_se
arrows(x0=sim_rho, y0=lower, x1=sim_rho, y1=upper, code=3, angle=90, length=0.05, col='blue')

legend(x=0.6, y=0.062, col=c('red', 'black', 'blue'),
       pch=c(NA, 16, 16), lty=c(1, NA, 1), legend=c("ANOVA",
                                      "ANCOVA",
                                      "Simulation"), cex=0.7)

