rm(list=ls())

library(pbapply)
source("../R/sim-analysis.R")
source("../R/sim-data.R")
source("../R/sim-procedure.R")
source("../R/sim-bounds.R")
source("../R/constants.R")
source("sim-utils.R")

# DEBUGGING AND TESTING SET TO FALSE
parallel <- TRUE

# GET TASK ID FROM SGE
if(parallel){
  TASKID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
  args <- commandArgs(trailingOnly=TRUE)
  OUT_DIR <- args[1]
  STD <- as.numeric(args[2])
  NSIMS <- as.integer(args[3])
} else {
  TASKID <- 70
  OUT_DIR <- "."
  STD <- 0.001
  NSIMS <- 10
}

# SET REPRODUCIBLE SEED
set.seed(0)
SEEDS <- sample(1:1000, size=1000, replace=FALSE)
set.seed(SEEDS[TASKID])

# GET INFORMATION FRACTION AND ALPHA SPENDING
# FUNCTION BASED ON ARGUMENTS
rates <- c(0.5, 1)

a.func <- spend(
  a=0.05,
  type="obf")

# CLOSURE FUNCTIONS
sim.data <- sim.data.closure(
  delta=0,
  beta=1,
  b0=1,
  cov_std=0.1,
  obs_std=STD)

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
trial_data <- replicate(NSIMS, sim.trial(sim.data), simplify=F)
result <- pblapply(trial_data, procedure)

# SAVE RESULTS
result2 <- condense.output(result)
save(result2, file=sprintf("%s/task_%s_std_%s.RData", OUT_DIR, TASKID, STD))
