# in_dir <- "~/Documents/FileZilla/rct/run-08-05-22-2/"
in_dir <- "~/Documents/FileZilla/rct/run-09-05-22-1/"
in_dir <- "~/Documents/FileZilla/rct/run-18-04-22-1/"
df <- fread(paste0(in_dir, "summary.csv"))

BTYPE <- "bias"

setnames(df, "bias_naive_med", "med_naive")
setnames(df, "bias_corr_med", "med_corr")

df.sub <- df[afunc == "pocock"]
df.sub[, dev_naive := get(paste0(BTYPE, "_naive")) * 100]
df.sub[, dev_corr := get(paste0(BTYPE, "_corr")) * 100]
df.sub <- df.sub[monitor == "anova" & final == "ancova"]
df.sub <- df.sub[est_var == TRUE]
df.sub[, est_var := NULL]

setorder(df.sub, delta, rho, n, type)

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
