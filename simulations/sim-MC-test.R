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
  STD <- 0.1
  NSIMS <- 200
}

# SET REPRODUCIBLE SEED
set.seed(0)
SEEDS <- sample(1:1000, size=1000, replace=FALSE)
set.seed(SEEDS[TASKID])

# GET INFORMATION FRACTION AND ALPHA SPENDING
# FUNCTION BASED ON ARGUMENTS
rates <- c(0.1, 0.5, 1)

a.func <- spend(
  a=0.05,
  type="obf")

BETA <- 1

# CLOSURE FUNCTIONS
sim.data <- sim.data.closure(
  delta=0,
  beta=10,
  b0=1,
  cov_std=0.1,
  obs_std=STD)

sim.trial <- sim.trial.closure(
  N=10000,
  rates=rates)

procedure <- procedure.closure(
  monitor="anova",
  final="ancova",
  correct=T,
  rates=rates,
  a.func=a.func)
  # sd_anova=(STD + 0.1**2*BETA**2)**0.5,
  # sd_ancova=STD)

# RUN SIMULATION
trial_data <- replicate(NSIMS, sim.trial(sim.data), simplify=F)
result <- pblapply(trial_data, procedure)

# SAVE RESULTS
result2 <- condense.output(result)
save(result2, file=sprintf("%s/task_%s_std_%s-known.RData", OUT_DIR, TASKID, STD))
