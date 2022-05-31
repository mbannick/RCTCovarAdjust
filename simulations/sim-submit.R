# SIMULATION SUBMISSION SCRIPT
# ----------------------------
library(magrittr)
library(data.table)
library(R.utils)

source("../simulations/parse-args.R")

# Read in command line arguments, with the following defaults
args <- commandArgs(
  trailingOnly=TRUE, asValues=TRUE,
  defaults=list(
    out_dir=".",
    n_sims="10",
    n="50,100,250,1000",
    delta="0.0,0.1",
    n_cov="1",
    rho="0.25,0.5",
    monitor="anova,ancova",
    final="anova,ancova",
    correct="FALSE,TRUE",
    afunc="obf,pocock",
    stages="3",
    ifracts="2",
    alpha="0.05",
    est_var="TRUE",
    est_bounds="FALSE,TRUE",
    desc="test"
  )
)

# Output directory
OUT_DIR <- args$out_dir
dir.create(OUT_DIR, recursive=TRUE)

ERROR <- sprintf("%s/output", OUT_DIR)
dir.create(ERROR, recursive=TRUE)

# Get the current git hash
hash <- system("git rev-parse HEAD", intern=TRUE)

# Write description to a text file
fileConn <- file(sprintf("%s/DESCRIPTION.txt", OUT_DIR))
writeLines(c(
  "DESCRIPTION: ",
  args$desc,
  "GIT HASH: ",
  hash,
  "DATE: ",
  format(Sys.time(), "%d-%m-%y-%H")
  ),
  fileConn)
close(fileConn)

# Number of simulations to run
N_SIMS <- as.integer(args$n_sims)

# Parameter grid from arguments
params <- list(
  n           = parse.args(args$n, as.integer),
  delta       = parse.args(args$delta, as.numeric),
  n_cov       = parse.args(args$n_cov, as.integer),
  rho         = parse.args(args$rho, as.numeric),
  monitor     = parse.args(args$monitor, as.character),
  final       = parse.args(args$final, as.character),
  correct     = parse.args(args$correct, as.logical),
  afunc       = parse.args(args$afunc, as.character),
  stages      = parse.args(args$stages, as.integer),
  ifracts     = parse.args(args$ifracts, as.integer),
  alpha       = parse.args(args$alpha, as.numeric),
  est_var     = parse.args(args$est_var, as.logical),
  est_bounds  = parse.args(args$est_bounds, as.logical)
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
