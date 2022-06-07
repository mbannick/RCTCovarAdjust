#$ -S /bin/sh

# Table 1
R --no-save --args run-23-05-22-1 pocock < "table-1.R"

# Figure 4 and 5
R --no-save --args run-28-05-22-2 < "figure-4-5.R"

# Table A1
R --no-save --args run-23-05-22-1 obf < "table-1.R"

# Table A2
R --no-save --args run-30-05-22-7 pocock < "table-1.R"

# Table A3
R --no-save --args run-30-05-22-8 pocock < "table-1.R"
