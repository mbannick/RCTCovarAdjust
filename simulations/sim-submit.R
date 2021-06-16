# SIMULATION SUBMISSION SCRIPT
# ----------------------------

# Read in command line arguments
args <- commandArgs(trailingOnly=TRUE)

# Output directory
OUT_DIR <- args[1]
dir.create(OUT_DIR, recursive=TRUE)

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
n_jobs <- expand.grid(params) %>% nrow
save(params, row.names=F, sprintf("%s/params.RData", OUT_DIR))

# Submit job arrays
command <- sprintf("qsub -cwd -t 1-%s shell.R sim.R %s %s",
                   n_jobs, OUT_DIR, n_sims)
system(command)
