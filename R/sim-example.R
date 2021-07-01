library(pbapply)
source("../R/sim-analysis.R")
source("../R/sim-data.R")
source("../R/sim-procedure.R")
source("../R/sim-bounds.R")
source("../R/constants.R")
source("sim-utils.R")

set.seed(0)

# GET INFORMATION FRACTION AND ALPHA SPENDING
# FUNCTION BASED ON ARGUMENTS
rates <- info.fractions(
  stages=3,
  type=2)

a.func <- spend(
  a=0.05,
  type="obf")

# CLOSURE FUNCTIONS
sim.data <- sim.data.closure(
  delta=0,
  beta=1,
  b0=1,
  cov_std=1,
  obs_std=1)

sim.trial <- sim.trial.closure(
  N=500,
  rates=rates)

procedure <- procedure.closure(
  monitor="anova",
  final="ancova",
  correct=TRUE,
  rates=rates,
  a.func=a.func)

# RUN SIMULATION
trial_data <- replicate(10, sim.trial(sim.data), simplify=F)
result <- pblapply(trial_data, procedure)

# SAVE RESULTS
result <- condense.output(result)
