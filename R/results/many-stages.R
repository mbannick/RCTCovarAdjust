in_dir <- "~/Documents/FileZilla/rct/run-16-05-22-2/"

df <- fread(paste0(in_dir, "summary.csv"))

setnames(df, "bias_naive_med", "med_naive")
setnames(df, "bias_corr_med", "med_corr")

df[, dev_naive := paste0(sprintf("%0.2f", med_naive*100), " (", sprintf("%0.2f", bias_naive*100), ")")]
df[, dev_corr := paste0(sprintf("%0.2f", med_corr*100), " (", sprintf("%0.2f", bias_corr*100), ")")]

df.sub <- df[, .(delta, n, rho, type, stage_1, dev_naive, dev_corr, cover_naive, cover_corr)]

df.sub.1 <- df.sub[type %in% c("inconsistent, naive", "inconsistent, corrected")] %>% data.table()
df.sub.2 <- df.sub[type %in% c("un-adjusted", "adjusted")] %>% data.table()

df.sub.1 <- dcast(delta + n + rho + stage_1 ~ type, value.var = c("dev_naive", "dev_corr", "cover_naive", "cover_corr"), data=df.sub.1)
df.sub.1 <- df.sub.1[, .(delta, n, rho, stage_1,
                       `dev_naive_inconsistent, naive`, `dev_corr_inconsistent, naive`, `dev_corr_inconsistent, corrected`,
                       `cover_naive_inconsistent, naive`, `cover_corr_inconsistent, naive`, `cover_corr_inconsistent, corrected`)]


df.sub.2[, type := NULL]
setnames(df.sub.1, c("delta", "rho", "n", "stage_1", "bias (0)", "bias (1)", "bias (2)", "cover (0)", "cover (1)", "cover (2)"))
setnames(df.sub.2, c("delta", "rho", "n", "stage_1", "bias (0)", "bias (2)", "cover (0)", "cover (2)"))

df.sub.1[, type := "switch"]
df.sub.2$type <- rep(c("anova", "ancova"), each=2)

dff <- rbind(df.sub.1, df.sub.2, fill=TRUE)
dff <- dff[, c("type", "delta", "rho", "n", "stage_1",
               "bias (0)", "bias (1)", "bias (2)",
               "cover (0)", "cover (1)", "cover (2)")]
dff$type <- factor(dff$type, levels=c("anova", "ancova", "switch"))
setorder(dff, type, delta)

tab <- xtable(dff, align=rep("c", 12), digits=3,
              caption=paste0("Simulation results for ", a,
                             " boundaries and monitoring and final test statistic ", tn,
                             ", estimating variance, 1000 simulations."))
print(tab, include.rownames=FALSE)

for(a in c("pocock")){
  for(tn in c("ancova", "anova")){
    df.sub <- df[rho %in% c(0.5) & afunc == a]
    df.sub[, dev_naive := get(paste0(BTYPE, "_naive"))]
    df.sub[, dev_corr := get(paste0(BTYPE, "_corr"))]
    df.sub[, bias_naive := dev_naive * 100]
    df.sub[, bias_corr := dev_corr * 100]
    df.sub <- df.sub[monitor == tn & final == tn]

    df.sub <- df.sub[est_var == TRUE]

    setorder(df.sub, delta, n)

    df.sub2 <- df.sub[, .(delta, n, bias_naive, bias_corr, cover_naive, cover_corr)]
    df.sub2 <- data.table(df.sub2)

    tab <- xtable(df.sub2, align=rep("c", 7), digits=3,
                  caption=paste0("Simulation results for ", a,
                                 " boundaries and monitoring and final test statistic ", tn,
                                 ", estimating variance, 1000 simulations."))
    print(tab, include.rownames=FALSE)
  }
}

#HACK
df.sub <- df[rho %in% c(0.5) & afunc == "pocock"]
df.sub[, bias_naive := bias_naive * 100]
df.sub[, bias_corr := bias_corr * 100]
df.sub <- df.sub[monitor == "anova" & final == "ancova"]

df.sub <- df.sub[est_var == TRUE]
df.sub <- df.sub[correct == TRUE]

setorder(df.sub, delta, n)

df.sub2 <- df.sub[, .(delta, n, bias_naive, bias_corr, cover_naive, cover_corr)]
df.sub2 <- data.table(df.sub2)

addtorow <- list()
addtorow$pos <- seq(5, nrow(df.sub2), by=5) %>% as.list
addtorow$command <- rep("\\hline \n", length(addtorow$pos))

tab <- xtable(df.sub2, align=rep("c", 7), digits=3,
              caption=paste0("Simulation results for ", afunc,
                             " boundaries and monitoring and final test statistic ", tn,
                             ", estimating variance, 1000 simulations."))
print(tab, include.rownames=FALSE,
      add.to.row = addtorow)

