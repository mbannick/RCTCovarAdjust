source("R/sim-analysis.R")
source("R/sim-data.R")
source("R/sim-procedure.R")
source("R/constants.R")

set.seed(160)

# HOW MANY SIMULATIONS TO PERFORM
N_REPS <- 10

# DEFINE CLOSURE FUNCTIONS
a.func <- OBF.SPEND(0.05)
sim.data <- sim.data.closure(delta=0.3, beta=1, b0=1,
                             cov_std=1, obs_std=1)
sim.trial <- sim.trial.closure(N=100, rates=1:4/4)
procedure <- procedure.closure(monitor="ancova", final="ancova",
                               rates=1:4/4, a.func=a.func)

# RUN SIMULATION
trial_data <- replicate(N_REPS, sim.trial(sim.data), simplify=F)
result <- lapply(trial_data, procedure)
