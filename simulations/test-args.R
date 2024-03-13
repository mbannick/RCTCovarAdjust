
library(R.utils)
args <- commandArgs(trailingOnly=TRUE, asValues=TRUE, adhoc=TRUE)
print(args)
print(args$n)
print(args$n + 1)
