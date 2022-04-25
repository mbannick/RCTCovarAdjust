in_dir <- "~/Documents/FileZilla/rct/run-10-04-22-2/"
in_dir <- "~/Documents/FileZilla/rct/run-11-04-22-3/"

df <- fread(paste0(in_dir, "summary.csv"))

df.sub <- df[rho %in% c(0.5) & afunc == "pocock"]
df.sub[, bias_naive := bias_naive * 100]
df.sub[, bias_corr := bias_corr * 100]
df.sub <- df.sub[monitor == "ancova" & final == "ancova"]

df.sub <- df.sub[est_var == TRUE]
df.sub <- df.sub[delta < 0.5]
df.sub <- df.sub[n < 10000]

setorder(df.sub, delta, n)

df.sub2 <- df.sub[, .(delta, n, bias_naive, bias_corr, cover_naive, cover_corr)]
df.sub2 <- data.table(df.sub2)

addtorow <- list()
addtorow$pos <- seq(4, nrow(df.sub2), by=4) %>% as.list
addtorow$command <- rep("\\hline \n", 3)

tab <- xtable(df.sub2, align=rep("c", 7), digits=3, caption="Simulation results.")
print(tab, include.rownames=FALSE,
      add.to.row = addtorow)
