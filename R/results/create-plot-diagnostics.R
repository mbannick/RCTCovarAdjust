rm(list=ls())

library(data.table)
library(ggplot2)
library(gridExtra)

VERSION <- "09-05-22-1"
in_dir <- sprintf("~/Documents/FileZilla/rct/run-%s/", VERSION)
out_dir <- paste0(in_dir, "/plot-diagnostics/", VERSION, "/")
dir.create(out_dir, showWarnings=FALSE, recursive=T)

params <- fread(paste0(in_dir, "params.csv"))

small_params <- params[est_var == TRUE]

# DELTA <- 0.0
# RHO <- 0.5
# N <- 1000
# AFUNC <- "pocock"
# CORRECT <- FALSE
# MONITOR <- "anova"
# FINAL <- "anova"

for(i in small_params$V1){
  print(i)
  subset <- small_params[V1 == i]
  DELTA <- subset$delta[1]
  RHO <- subset$rho[1]
  N <- subset$n[1]
  AFUNC <- subset$afunc[1]
  CORRECT <- subset$correct[1]
  MONITOR <- subset$monitor[1]
  FINAL <- subset$final[1]

  texts <- paste0(
    "delta: ", DELTA,
    ",\nrho: ", RHO,
    ",\nn: ", N,
    ",\n", AFUNC,
    ",\n", MONITOR, "-", FINAL
  )
  if(MONITOR != FINAL){
    if(CORRECT){
      texts <- paste0(
        texts,
        ", corrected"
      )
    } else {
      texts <- paste0(
        texts,
        ", uncorrected"
      )
    }
  }

  filename <- paste(
    paste0("delta", DELTA),
    paste0("rho", RHO),
    N, AFUNC, MONITOR, FINAL, CORRECT,
    sep="-"
  )
  outfile <- paste0(out_dir, filename, ".pdf")

  # subset <- subset[delta == DELTA & rho == RHO & n == N]
  # subset <- subset[afunc == AFUNC & correct == CORRECT]
  # subset <- subset[monitor == MONITOR & final == FINAL]
  # print(subset$V1)
  # if(length(subset$V1) > 1) stop()

  load(paste0(in_dir, sprintf("params_%s.RData", subset$V1)))

  bounds <- result$u.bounds

  # x <- c(1, 2, 0)
  numbounds <- function(vec) sum(vec != 0)
  stages <- apply(bounds, MARGIN=1, FUN=numbounds)

  result$reject
  stages[(result$reject == FALSE) & (stages == 3)] <- NA

  proportions(table(stages))

  df <- data.table(
    adj_estimate=result$point,
    estimate=result$est,
    final_tstat=result$tstat,
    stage=stages
  )
  df[, sim := .I]
  melt(df, id.vars=c("sim", "stage", "final_tstat"))

  p0 <- ggplot() +
    annotate("text", x = 1, y = 1, size = 8,
             label = texts) +
    theme_void()

  (p1 <- ggplot(data=df) +
    geom_histogram(aes(x=final_tstat), binwidth=0.15) +
    facet_wrap(~stage) +
    geom_vline(xintercept=0, col='red')
  )

  (p2 <- ggplot(data=df) +
    geom_density(aes(x=estimate), color="blue") +
    geom_density(aes(x=adj_estimate), color="red") +
    facet_wrap(~stage, scales="free") +
    geom_vline(xintercept=DELTA, col='black', linetype='dashed')
  )

  (p3 <- ggplot(data=df) +
    geom_density(aes(x=estimate), color="blue") +
    geom_density(aes(x=adj_estimate), color="red") +
    geom_vline(xintercept=DELTA, col='black', linetype='dashed') +
    geom_vline(xintercept=mean(df$estimate), color='blue', linetype='dashed') +
    geom_vline(xintercept=mean(df$adj_estimate), color='red', linetype='dashed')
  )

  pdf(outfile, height=6, width=9)
  p <- grid.arrange(p0, p1, p2, p3, nrow=2)
  p
  dev.off()
}
