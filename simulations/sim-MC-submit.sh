#!/bin/sh

REPO="/home/students/mnorwood/repos/RCTCovarAdjust"
RESULTS="/home/students/mnorwood/rct"

OUT="${RESULTS}/sgeout/"
SHELL="${REPO}/simulations/shell.sh"
SCRIPT="${REPO}/simulations/sim-MC-test.R"
OUTPUT="${RESULTS}/simulation-MC-test/"

# 800,000 simulations
NSIMS=10000
TASKS=80

# Standard deviations that correspond to
# rho of 0.01 and 0.7546, which maximize the difference
# between the expected value and 0.05.
STD1=0.001
STD2=0.115

CONSTANTS="-cwd -N sim -j y -o ${OUT} -pe smp 1 -t 1-${TASKS} -q normal.q ${SHELL} ${SCRIPT}"
ARGS1="${OUT} ${NSIMS} ${STD1}"
ARGS2="${OUT} ${NSIMS} ${STD2}"

CMD1="${CONSTANTS} ${ARGS1}"
CMD2="${CONSTANTS} ${ARGS2}"

echo qsub ${CMD1}
qsub ${CMD1}

echo qsub ${CMD2}
qsub ${CMD2}
