## Create Tables and Figures

After running simulations using the files in `../simulations/`, this folder contains scripts to create the tables and figures. Each requires a versioned input file with a `summary.csv` file. If you run this, you will need to change the hard-coded base directory to your own device location.

-   `table-1.R`: Creates main text simulation table results (table 1), or any of the appendix tables A1-A3. Depends on which version of the simulations you pass in, and whether you select `obf` of `pocock` as the argument for the type of boundaries.

-   `figure-2-3.R`: Creates figures 2 and 3 in the preprint, not based on simulations.

-   `figure-4-5.R`: Creates the type I error and power simulation figures.

-   `figure-A1.R`: Creates appendix figure A1 on misspecification of $\rho$ from a simulation run.

-   `make-tables-figures.sh`: This file has the versions of simulations that were run for the preprint, and creates all of the above simulation tables and figures by calling R scripts (besides `figure-2-3.R`).
