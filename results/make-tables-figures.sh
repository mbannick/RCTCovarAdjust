#$ -S /bin/sh

# Table 1
R --no-save --args run-25-08-22-1/table1A1 pocock < "table-1.R"

# Figure 4 and 5
R --no-save --args run-25-08-22-1/figure1 < "figure-4-5.R"

# Table A1
R --no-save --args run-25-08-22-1/table1A1 obf < "table-1.R"

# Table A2
R --no-save --args run-25-08-22-1/tableA2 pocock < "table-1.R"

# Table A3
R --no-save --args run-25-08-22-1/tableA3 pocock < "table-1.R"

# Figure A1
R --no-save --args run-25-08-22-1/figureA1 < "figure-A1.R"
