df.sub <- df[rho %in% c(0.1, 0.5) & afunc == "pocock"]
df.sub[, bias_naive := bias_naive * 100]
df.sub[, bias_corr := bias_corr * 100]
df.sub <- df.sub[monitor == "anova" & final == "ancova"]

setorder(df.sub, delta, rho, n, type)

df.sub2 <- df.sub[, .(delta, rho, n, type, bias_naive, bias_corr, cover_naive, cover_corr)]
df.sub2 <- data.table(df.sub2)
df.sub2 <- dcast(delta + rho + n ~ type, value.var = c("bias_naive", "bias_corr", "cover_naive", "cover_corr"), data=df.sub2)
df.sub2 <- df.sub2[, .(delta, rho, n,
                       `bias_naive_inconsistent, naive`, `bias_corr_inconsistent, naive`, `bias_corr_inconsistent, corrected`,
                       `cover_naive_inconsistent, naive`, `cover_corr_inconsistent, naive`, `cover_corr_inconsistent, corrected`)]

setnames(df.sub2, c("delta", "rho", "n", "bias (0)", "bias (1)", "bias (2)", "cover (0)", "cover (1)", "cover (2)"))

addtorow <- list()
addtorow$pos <- seq(4, nrow(df.sub2), by=4) %>% as.list
addtorow$command <- rep("\\hline \\\\\n", 6)

tab <- xtable(df.sub2, align=rep("c", 10), digits=2, caption="Simulation results.")
print(tab, include.rownames=FALSE,
      add.to.row = addtorow)
