#$ -S /bin/sh

# Table 1
R --no-save --args run-23-05-22-1 pocock < "table-1.R"

# Table A1
R --no-save --args run-23-05-22-1 obf < "table-1.R"

# Table A2
# R --no-save --args . < "table-1.R"

# Table A3
# R --no-save --args . < "table-1.R"
