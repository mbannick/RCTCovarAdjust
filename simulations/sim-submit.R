# SIMULATION SUBMISSION SCRIPT
# ----------------------------
library(magrittr)

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
  n=c(100, 500),
  delta=c(0.0, 0.1, 0.5),
  beta=c(0, 0.5, 1.5),
  monitor=c("anova", "ancova"),
  final=c("anova", "ancova"),
  correct=c(FALSE, TRUE),
  afunc=c("obf", "pocock"),
  stages=c(4),
  ifracts=c(3),
  alpha=c(0.05),
  cov_std=c(1),
  obs_std=c(1),
  intercept=c(1)
)

# Save parameter list and number of tasks
N_JOBS <- expand.grid(params) %>% nrow
save(params, file=sprintf("%s/params.RData", OUT_DIR))

# Submit job arrays
command <- sprintf("qsub -cwd -j y -o %s -t 1-%s -q normal.q shell.sh sim.R %s %s",
                   ERROR, N_JOBS, OUT_DIR, N_SIMS)
print(command)
system(command)
