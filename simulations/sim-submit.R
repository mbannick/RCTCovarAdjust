# SIMULATION SUBMISSION SCRIPT
# ----------------------------
library(magrittr)
library(data.table)

# Read in command line arguments
args <- commandArgs(trailingOnly=TRUE)

# Output directory
OUT_DIR <- args[1]
dir.create(OUT_DIR, recursive=TRUE)

ERROR <- sprintf("%s/output", OUT_DIR)
dir.create(ERROR, recursive=TRUE)

# Number of simulations to run
N_SIMS <- args[2]

# Parameter grid
params <- list(
  n=c(100),
  delta=c(0.0, 0.25),
  n_cov=c(1),
  rho=c(0.5),
  monitor=c("anova", "ancova"),
  final=c("anova", "ancova"),
  correct=c(FALSE, TRUE),
  afunc=c("pocock"),
  stages=c(4),
  ifracts=c(2),
  alpha=c(0.05),
  est_var=c(TRUE)
)

# Save parameter list and number of tasks
param_grid <- expand.grid(params)
param_grid <- data.table(param_grid)
param_grid <- param_grid[!(monitor == "anova" & final == "anova" & correct == T)]
param_grid <- param_grid[!(monitor == "ancova" & final == "ancova" & correct == T)]
param_grid <- param_grid[!(monitor == "ancova" & final == "anova")]

N_JOBS <- nrow(param_grid)
write.csv(param_grid, file=sprintf("%s/params.csv", OUT_DIR))

# Submit job arrays
command <- sprintf("qsub -cwd -j y -o %s -t 1-%s -q normal.q shell.sh sim.R %s %s",
                   ERROR, N_JOBS, OUT_DIR, N_SIMS)
print(command)
system(command)
