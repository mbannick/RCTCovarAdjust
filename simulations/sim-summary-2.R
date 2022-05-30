rm(list=ls())
library(data.table)
library(ggplot2)
library(magrittr)

args <- commandArgs(TRUE)
in_dir <- args[1]
# in_dir <- "~/Documents/FileZilla/rct/run-10-04-22-2/"
# in_dir <- "~/Documents/FileZilla/rct/run-11-04-22-3/"
# in_dir <- "~/Documents/FileZilla/rct/run-08-05-22-2/"
# in_dir <- "~/Documents/FileZilla/rct/run-09-05-22-1/"
# in_dir <- "~/Documents/FileZilla/rct/run-09-05-22-1/"
# in_dir <- "~/Documents/FileZilla/rct/run-18-04-22-1"
# in_dir <- "~/Documents/FileZilla/rct/run-16-05-22-2/"
setwd(in_dir)

files <- list.files(in_dir, ".RData")

parameters <- fread("params.csv")
setnames(parameters, "V1", "param")

stage_num <- unique(parameters$stages)
if(length(stage_num) > 1) stop()
potential_stages <- paste0("stage_", 1:stage_num)

get.results <- function(file){
  param_id <- gsub("params_", "", file)
  param_id <- gsub(".RData", "", param_id) %>% as.integer()
  print(param_id)

  truth <- c(parameters[param == param_id, delta])

  load(file)
  power <- mean(result$reject)
  bias_naive <- mean(result$est - truth)
  bias_corr <- mean(result$point - truth)
  bias_naive_med <- median(result$est - truth)
  bias_corr_med <- median(result$point - truth)

  cover_naive <- mean((result$naive_ci[, 1] <= truth) & (result$naive_ci[, 2] >= truth))
  cover_corr <- mean(result$ci[, 1] <= truth & (result$ci[, 2] >= truth))

  bounds <- result$u.bounds
  numbounds <- function(vec) sum(vec != 0)
  stages <- apply(bounds, MARGIN=1, FUN=numbounds)

  stageprop <- function(num) mean(stages == num)
  stageprops <- sapply(1:stage_num, stageprop)

  return(c(param_id, power,
           bias_naive, bias_corr,
           bias_naive_med, bias_corr_med,
           cover_naive, cover_corr,
           stageprops))
}

results <- sapply(files, get.results)
results <- t(results)
results <- data.table(results)
setnames(results, c("param", "power", "bias_naive", "bias_corr",
                    "bias_naive_med", "bias_corr_med",
                    "cover_naive", "cover_corr",
                    potential_stages))

df <- merge(parameters, results, by="param")
# df[, std := sqrt(power * (1 - power)) / sqrt(1000)]
# df[, lower := power - qnorm(0.975) * std]
# df[, upper := power + qnorm(0.975) * std]

df[monitor == "anova" & final == "anova", type := "un-adjusted"]
df[monitor == "ancova" & final == "ancova", type := "adjusted"]
df[monitor == "anova" & final == "ancova" & correct == TRUE, type := "inconsistent, corrected"]
df[monitor == "anova" & final == "ancova" & correct == FALSE, type := "inconsistent, naive"]
df[, type := factor(type, levels=c("un-adjusted", "adjusted",
                                   "inconsistent, naive", "inconsistent, corrected"))]

df[, n := as.factor(n)]
write.csv(df, file=paste0(in_dir, "/summary.csv"), row.names=F)

