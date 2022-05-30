# in_dir <- "~/Documents/FileZilla/rct/run-08-05-22-2/"
# in_dir <- "~/Documents/FileZilla/rct/run-09-05-22-1/"
in_dir <- "~/Documents/FileZilla/rct/run-18-04-22-1/"
df <- fread(paste0(in_dir, "summary.csv"))

setnames(df, "bias_naive_med", "med_naive")
setnames(df, "bias_corr_med", "med_corr")

df.sub <- df[afunc == "obf"]

df.sub[, bias_naive := bias_naive * 100]
df.sub[, bias_corr := bias_corr * 100]
df.sub[, med_naive := med_naive * 100]
df.sub[, med_corr := med_corr * 100]

df.sub <- df.sub[monitor == "anova" & final == "ancova"]
df.sub <- df.sub[est_var == TRUE]
df.sub[, est_var := NULL]

df.sub[, neworder := sqrt(n) * delta]


setorder(df.sub, delta, rho, n, type)

df.sub[, dev_naive := paste0(sprintf("%0.2f", med_naive), " (", sprintf("%0.2f", bias_naive), ")")]
df.sub[, dev_corr := paste0(sprintf("%0.2f", med_corr), " (", sprintf("%0.2f", bias_corr), ")")]

df.sub2 <- df.sub[, .(delta, rho, n, type, dev_naive, dev_corr, cover_naive, cover_corr)]
df.sub2 <- data.table(df.sub2)
df.sub2 <- dcast(delta + rho + n ~ type, value.var = c("dev_naive", "dev_corr", "cover_naive", "cover_corr"), data=df.sub2)
df.sub2 <- df.sub2[, .(delta, rho, n,
                       `dev_naive_inconsistent, naive`, `dev_corr_inconsistent, naive`, `dev_corr_inconsistent, corrected`,
                       `cover_naive_inconsistent, naive`, `cover_corr_inconsistent, naive`, `cover_corr_inconsistent, corrected`)]

setnames(df.sub2, c("delta", "rho", "n", "bias (0)", "bias (1)", "bias (2)", "cover (0)", "cover (1)", "cover (2)"))

addtorow <- list()
addtorow$pos <- seq(4, nrow(df.sub2), by=4) %>% as.list
addtorow$command <- rep("\\hline \n", length(addtorow$pos))

tab <- xtable(df.sub2, align=rep("c", 10), digits=2, caption="Simulation results.")
print(tab, include.rownames=FALSE,
      add.to.row = addtorow)

df.sub3 <- df.sub[, .(stage_1, delta, rho, n, type, med_naive, bias_naive, med_corr, bias_corr, cover_naive, cover_corr)]
df.sub3 <- data.table(df.sub3)
df.sub3 <- dcast(stage_1 + delta + rho + n ~ type, value.var = c("med_naive", "med_corr","bias_naive", "bias_corr", "cover_naive", "cover_corr"), data=df.sub3)
df.sub3 <- df.sub3[, .(stage_1, delta, rho, n,
                       `med_naive_inconsistent, naive`,
                       `med_corr_inconsistent, naive`,
                       `med_corr_inconsistent, corrected`,
                       `bias_naive_inconsistent, naive`,
                       `bias_corr_inconsistent, naive`,
                       `bias_corr_inconsistent, corrected`,
                       `cover_naive_inconsistent, naive`, `cover_corr_inconsistent, naive`, `cover_corr_inconsistent, corrected`)]

setnames(df.sub3, c("stage_1", "delta", "rho", "n",
                    "med (0)", "med (1)", "med (2)",
                    "bias (0)", "bias (1)", "bias (2)",
                    "cover (0)", "cover (1)", "cover (2)"))

df.sub3[, neworder := sqrt(n) * delta]

VAL <- "neworder"
df.sub3[, xcol := get(VAL)]

med3 <- df.sub3[, c("xcol", "med (0)", "med (1)", "med (2)")]
setnames(med3, c("xcol", "Simple", "GS", "GS + Adjust"))
bias3 <- df.sub3[, c("xcol", "bias (0)", "bias (1)", "bias (2)")]
setnames(bias3, c("xcol", "Simple", "GS", "GS + Adjust"))

med3 <- melt(med3, id.vars="xcol")
med3[, vartype := "median"]
bias3 <- melt(bias3, id.vars="xcol")
bias3[, vartype := "bias"]

total3 <- rbind(med3, bias3)
total3[, value := value / 100]

dir.create(paste0(in_dir, "/figures"), showWarnings=FALSE)
pdf(paste0(in_dir, "/figures/", "std-effect.pdf"), height=5, width=7)
ggplot(data=total3, aes(x=xcol, y=value, color=variable,
                        group=variable)) +
  facet_wrap(~ vartype) + geom_point(alpha=0.5)
dev.off()
