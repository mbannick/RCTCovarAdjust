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
  n=c(50, 100, 250, 1000),
  delta=c(0.0, 0.1),
  n_cov=c(1),
  rho=seq(0.1, 0.9, by=0.1),
  monitor=c("anova", "ancova"),
  final=c("anova", "ancova"),
  correct=c(FALSE, TRUE),
  afunc=c("obf", "pocock"),
  stages=c(3),
  ifracts=c(2),
  alpha=c(0.05),
  est_var=c(TRUE),
  est_bounds=c(FALSE, TRUE)
)

# Save parameter list and number of tasks
param_grid <- expand.grid(params)
param_grid <- data.table(param_grid)
param_grid <- param_grid[!(monitor == "anova" & final == "anova" & correct == T)]
param_grid <- param_grid[!(monitor == "ancova" & final == "ancova" & correct == T)]
param_grid <- param_grid[!(monitor == "ancova" & final == "anova")]
param_grid <- param_grid[!(monitor == "anova" & final == "anova" & est_bounds == FALSE)]
param_grid <- param_grid[!(monitor == "ancova" & final == "ancova" & est_bounds == FALSE)]
param_grid <- param_grid[!(monitor == "anova" & final == "ancova" & est_bounds == FALSE & correct == FALSE)]

N_JOBS <- nrow(param_grid)
write.csv(param_grid, file=sprintf("%s/params.csv", OUT_DIR))

# Submit job arrays
command <- sprintf("qsub -cwd -j y -o %s -t 1-%s -q normal.q shell.sh sim.R %s %s",
                   ERROR, N_JOBS, OUT_DIR, N_SIMS)
print(command)
system(command)
