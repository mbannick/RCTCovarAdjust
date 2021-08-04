#!/bin/sh

REPO="/home/users/mnorwood/repos/RCTCovarAdjust"
RESULTS="/home/users/mnorwood/rct"

OUT="${RESULTS}/sgeout"
SHELL="${REPO}/simulations/shell.sh"
SCRIPT="${REPO}/simulations/sim-MC-test.R"
OUTPUT="${RESULTS}/simulation-MC-test"

# 800,000 simulations
NSIMS=10000
TASKS=80

# Standard deviations that correspond to
# rho of 0.01 and 0.7546, which maximize the difference
# between the expected value and 0.05.

for STD in 0.01 0.02 0.03 0.04 0.05 0.06 0.08 0.1 0.2 0.3 0.75
do
    CONSTANTS="-N sim -j y -o ${OUT} -pe smp 1 -t 1-${TASKS} -q normal.q ${SHELL} ${SCRIPT}"
    ARGS="${OUTPUT} ${STD} ${NSIMS}"
    CMD="${CONSTANTS} ${ARGS}"

    echo qsub -cwd ${CMD}
    qsub -cwd ${CMD}
done

