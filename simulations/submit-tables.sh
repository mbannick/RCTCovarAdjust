#$ -S /bin/sh

ALL="-n_sims ${1}"

# Table 1 and A1
R-4.0.1 --no-save --args ${ARGS} -desc "table 1 and A1 -- main results" < "sim-submit.R"

# Table A2
R-4.0.1 --no-save --args ${ARGS} -desc "table A2, 2 stages + pocock" -stages 2 -afunc pocock < "sim-submit.R"

# Table A3
R-4.0.1 --no-save --args ${ARGS} -desc "table A3, 2 covariates + pocock" -n_cov 2 -afunc pocock < "sim-submit.R"

