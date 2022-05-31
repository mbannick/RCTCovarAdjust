rm(list=ls())

source("../R/sim-analysis.R")
source("../R/sim-data.R")
source("../R/sim-procedure.R")
source("../R/sim-bounds.R")
source("../R/constants.R")
source("../simulations/sim-utils.R")

# WHEN DEBUGGING AND TESTING SET TO FALSE
parallel <- TRUE

# GET TASK ID FROM SGE
if(parallel){
  args <- commandArgs(trailingOnly=TRUE)

  TASKID <- as.numeric(Sys.getenv("SGE_TASK_ID"))
  OUT_DIR <- args[1]
  N_SIMS <- as.integer(args[2])

  # Read in parameter grid
  param_grid <- read.csv(sprintf("%s/params.csv", OUT_DIR))
} else {
  TASKID <- 1
  OUT_DIR <- "."
  N_SIMS <- 50

  # Create a one-row parameter grid to test out examples
  param_grid <- data.frame(
    n=50,
    delta=0.0,
    n_cov=1,
    rho=0.25,
    monitor="anova",
    final="ancova",
    correct=TRUE,
    afunc="obf",
    stages=3,
    ifracts=2,
    alpha=0.05,
    est_var=TRUE,
    est_bounds=FALSE
  )
}

# SET REPRODUCIBLE SEED
set.seed(715)

# PARAMETER GETTER FOR THIS TASK ID
gp <- get.param.closure(TASKID, param_grid)

# GET INFORMATION FRACTION AND ALPHA SPENDING
# FUNCTION BASED ON ARGUMENTS
rates <- info.fractions(
  stages=gp("stages"),
  type=gp("ifracts"))$tk

a.func <- spend(
  a=gp("alpha"),
  type=gp("afunc"))

# BOUNDARY FUNCTION -- WILL BE ESTIMATED USING ALPHA-SPENDING IF EST_BOUNDS,
# OTHERWISE, IT'S BASED ON THE ALPHA-SPENDING FUNCTIONS IN WASSMER-BRANNATH
b.func <- get.boundary.closure(
  a.func=a.func,
  rates=rates,
  est.bounds=gp("est_bounds"),
  a.type=gp("afunc"),
  rho=gp("rho"),
  n=gp("n"))

# CLOSURE FUNCTIONS
sim.data <- sim.data.closure(
  delta=gp("delta"),
  rho=gp("rho"),
  n_cov=gp("n_cov"))

# PROCEDURE FOR SIMULATING THE TRIAL DATA BASED ON INFORMATION FRACTIONS
sim.trial <- sim.trial.closure(
  N=gp("n"),
  rates=rates)

# PROCEDURE FOR ESTIMATING THE VARIANCE, EITHER YOU HAVE KNOWN VALUES
# OR YOU RETURN NANS SO THAT THEY WILL BE ESTIMATED
v.func <- variance.closure(
  rho=gp("rho"),
  est_var=gp("est_var"))

procedure <- procedure.closure(
  monitor=gp("monitor"),
  final=gp("final"),
  correct=gp("correct"),
  est.bounds=gp("est_bounds"),
  rates=rates,
  a.func=a.func,
  v.func=v.func,
  b.func=b.func)

# RUN SIMULATION
trial_data <- replicate(N_SIMS, sim.trial(sim.data), simplify=F)
result <- lapply(trial_data, procedure)

# SAVE RESULTS
result <- condense.output(result)
save(result, file=sprintf("%s/params_%s.RData", OUT_DIR, TASKID))
