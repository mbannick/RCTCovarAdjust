#$ -S /bin/sh

ARGS="-n_sims ${2} -out_dir ~/rct/${1}"

# Table 1 and A1
# R-4.0.1 --no-save --args ${ARGS} -desc "table 1 and A1 -- main results" < "sim-submit.R"

# Table A2
# R-4.0.1 --no-save --args ${ARGS} -desc "table A2, 2 stages + pocock" -stages 2 -afunc pocock < "sim-submit.R"

# Table A3
# R-4.0.1 --no-save --args ${ARGS} -desc "table A3, 2 covariates + pocock" -n_cov 2 -afunc pocock < "sim-submit.R"

# Figure 1
# EXTRA="-delta 0.0,0.1 -rho 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 -est_bounds FALSE,TRUE"
# R-4.0.1 --no-save --args ${ARGS} -desc "figure 1 -- main results" ${EXTRA} < "sim-submit.R"

# Figure A1
EXTRA="-delta 0.0,0.1 -rho 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 -design_rho 0.2,0.4,0.6,0.8,1.0 -est_bounds FALSE"
R-4.0.1 --no-save --args ${ARGS} -desc "figure A1 -- misspecified rho" ${EXTRA} < "sim-submit.R"
