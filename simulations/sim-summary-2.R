rm(list=ls())
library(data.table)
library(ggplot2)

args <- commandArgs(TRUE)
in_dir <- args[1]
setwd(in_dir)

files <- list.files(in_dir, ".RData")

get.power <- function(file){
  param <- gsub("params_", "", file)
  param <- gsub(".RData", "", param) %>% as.integer()
  print(param)

  load(file)
  power <- mean(result$reject)
  return(c(param, power))
}

results <- sapply(files, get.power)
results <- t(results)
results <- data.table(results)
setnames(results, c("param", "power"))

parameters <- fread("params.csv")
setnames(parameters, "V1", "param")

df <- merge(parameters, results, by="param")
df[, std := sqrt(power * (1 - power)) / sqrt(1000)]
df[, lower := power - qnorm(0.975) * std]
df[, upper := power + qnorm(0.975) * std]

df[monitor == "anova" & final == "anova", type := "un-adjusted"]
df[monitor == "ancova" & final == "ancova", type := "adjusted"]
df[monitor == "anova" & final == "ancova" & correct == TRUE, type := "inconsistent, corrected"]
df[monitor == "anova" & final == "ancova" & correct == FALSE, type := "inconsistent, naive"]
df[, type := factor(type, levels=c("un-adjusted", "adjusted",
                                   "inconsistent, naive", "inconsistent, corrected"))]

df[, n := as.factor(n)]
write.csv(df, file=paste0(in_dir, "/summary.csv"), row.names=F)

