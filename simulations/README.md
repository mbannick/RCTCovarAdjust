## Run Simulations

To run the simulations for the preprint, run `sh submit-script.sh`. Specifically, this is what each file does:

-   `sim.R`: Runs `N_SIMS` number of simulations for one setting. Can be run in parallel or locally. For local runs, can adjust the settings in `param_grid`, and set `parallel=FALSE`. Uses `sim-utils.R` to get output from simulations into a saveable format.
-   `sim-submit.R`: This R script submits combinations of parameter settings as a job array, and each job calls `sim.R` with a different parameter setting. Uses a custom command-line argument parser in `parse-args.R`. Uses `shell.sh` to pass arguments to `sim.R`.
-   `submit-script.sh`: This shell script calls `sim-submit.R` several times with the parameter combinations required to reproduce the tables and figures in the main text and supplemental appendix of the preprint.
-   `sim-summary.R`: This R script summarizes the results from one run of `sim-submit.R`. It takes in a directory with `.RData` outputs from each time `sim.R` is called and outputs a `summary.csv` file.
-   `copy-files.sh`: This shell script copies the summary files from a remote server to a local device. The raw files are very large if thousands of simulations are run, so recommend only working with the summary files locally.
