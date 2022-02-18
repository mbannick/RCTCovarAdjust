rm(list=ls())
library(data.table)
library(ggplot2)

in_dir <- "~/Documents/FileZilla/rct/run-17-02-22-5/"
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

# Type I Error Plot

pdf("~/OneDrive/Documents/2021-2022/BIOST 600-noah/type1error.pdf", height=8, width=14)
ggplot(data=df[delta == 0.0 & stages == 3], aes(x=rho, y=power, color=type, linetype=afunc)) +
  facet_grid(~ n) +
  geom_point() +
  geom_hline(yintercept=0.05, linetype="dashed", color="black") +
  geom_line() +
  theme_minimal() +
  scale_color_brewer(palette="Set1") +
  theme(legend.position="top")
  # geom_errorbar(aes(ymin=lower, ymax=upper))
dev.off()

ggplot(data=df[delta > 0 & stages == 3], aes(x=rho, y=power, color=n, linetype=afunc)) +
  facet_grid(type ~ delta) +
  geom_point() + geom_line()

pdf("~/OneDrive/Documents/2021-2022/BIOST 600-noah/power.pdf", height=8, width=14)
ggplot(data=df[delta == 0.1 & stages == 3], aes(x=rho, y=power, color=type, linetype=afunc)) +
  facet_grid( ~ n) +
  geom_point() + geom_line() +
  scale_color_brewer(palette="Set1") +
  theme_minimal() +
  theme(legend.position="top")
dev.off()
