rm(list=ls())

library(pbapply)
setwd("~/repos/RCTCovarAdjust/simulations/")

source("../R/sim-analysis.R")
source("../R/sim-data.R")
source("../R/sim-procedure.R")
source("../R/sim-bounds.R")
source("../R/constants.R")
source("sim-utils.R")

# SET REPRODUCIBLE SEED
set.seed(715)

# GET INFORMATION FRACTION AND ALPHA SPENDING
# FUNCTION BASED ON ARGUMENTS
rates <- c(0.5, 1)

a.func <- spend(
  a=0.05,
  type="obf")

stds <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.2, 0.3, 0.75)
stds <- c(0, stds)
rho <- function(x) sqrt(x**2 / (0.1**2 + x**2))
# sim.datas <- lapply(stds, function(x)
#   sim.data.closure(
#     delta=0,
#     beta=1,
#     b0=1,
#     cov_std=0.1,
#     obs_std=x)
# )

sim.data <- sim.data.closure(
  delta=0,
  beta=1,
  b0=1,
  cov_std=0.1,
  obs_std=gp("obs_std"))

sim.trial <- sim.trial.closure(
  N=100,
  rates=rates)

procedure <- procedure.closure(
  monitor=F,
  final=T,
  correct=F,
  rates=rates,
  a.func=a.func)

# RUN SIMULATION
trial_data <- replicate(1000, sim.trial(sim.data), simplify=F)
result <- pblapply(trial_data, procedure)

# SAVE RESULTS
result2 <- condense.output(result)

# COMMENTED OUT ON THURS JULY 29
# trial.data <- function(x) replicate(100, sim.trial(x), simplify=F)
# run.result <- function(x) {
#   trial_data <- trial.data(x)
#   result <- pblapply(trial_data, procedure)
#   result <- condense.output(result)
#   return(mean(result$reject))
# }
#
# empirical.ps <- lapply(sim.datas, run.result)
# ps <- unlist(empirical.ps)
# rhos <- rho(stds)
#
# vals <- cbind(rho=rhos, std=stds, ps=ps)
# print(vals)
# save(vals, file="empirical-ps-hack.RData")

# u_k <- c(2.963, 1.969)
# u_k1 <- cbind(-u_k, u_k)
#
# # This gives exactly alpha = 0.05 by putting in the last
# # boundary.
# n_k1 <- c(10, 20)

# get.pvalue.sw(obs=-1.969, u_k=t(matrix(u_k1[1:1,])), n_k=n_k1,
#               ancova_monitor=F, ancova_test=F, last_stage=T)
# # 0.05000002
# get.pvalue.sw(obs=-1.969, u_k=t(matrix(u_k1[1:1,])), n_k=n_k1, rho=0.99,
#               ancova_monitor=F, ancova_test=T, last_stage=T)
# # 0.06235481
