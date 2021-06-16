rm(list=ls())

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
} else {
  TASKID <- 1
}

# GET COMMAND LINE ARGS
args <- commandArgs(trailingOnly=TRUE)
OUT_DIR <- args[1]
N_SIMS <- as.integer(args[2])

# SET REPRODUCIBLE SEED
set.seed(715)

# STATIC PARAMS
N_REPS <- 10

# PARAMETER GRID
load(sprintf("%s/params.RData", OUT_DIR))

# PARAMETER GETTER FOR THIS TASK ID
param_grid <- expand.grid(params)
gp <- get.param.closure(TASKID, param_grid)

# GET INFORMATION FRACTION AND ALPHA SPENDING
# FUNCTION BASED ON ARGUMENTS
rates <- info.fractions(
  stages=gp("stages"),
  type=gp("ifracts"))

a.func <- spend(
  a=gp("alpha"),
  type=gp("afunc"))

# CLOSURE FUNCTIONS
sim.data <- sim.data.closure(
  delta=gp("delta"),
  beta=gp("beta"),
  b0=gp("intercept"),
  cov_std=gp("cov_std"),
  obs_std=gp("obs_std"))

sim.trial <- sim.trial.closure(
  N=gp("n"),
  rates=rates)

procedure <- procedure.closure(
  monitor=gp("monitor"),
  final=gp("final"),
  correct=gp("correct"),
  rates=rates,
  a.func=a.func)

# RUN SIMULATION
trial_data <- replicate(N_SIMS, sim.trial(sim.data), simplify=F)
result <- lapply(trial_data, procedure)

# SAVE RESULTS
result <- condense.output(result)
save(result, file=sprintf("%s/params_%s.RData", OUT_DIR, TASKID))
