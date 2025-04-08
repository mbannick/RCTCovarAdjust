rm(list=ls())
library(data.table)
library(ggplot2)
library(magrittr)

in_dir <- "/Users/marlena/Documents/FileZilla/rct/run-17-02-22-5"
setwd(in_dir)

# results <- sapply(files, get.power)
# results <- t(results)
# results <- data.table(results)
# setnames(results, c("param", "power"))
#
parameters <- fread(paste0(in_dir, "/params.csv"))
setnames(parameters, "V1", "param")

parms <- parameters[monitor == "ancova" &
                      final == "ancova" &
                      delta == 0.5 &
                      afunc == "pocock" &
                      n == 250 &
                      stages == 3 &
                      rho == 0.5]

results <- list()
for(p in 1:length(parms$param)){
  load(paste0("params_", parms$param[p], ".RData"))
  results[[p]] <- result$tstat
}

# DELTA = 0
# > parms$param
# [1] 487 496 505 514 523 532 541 550 559

# DELTA = 0.1
# > parms$param
# [1] 490 499 508 517 526 535 544 553 562
