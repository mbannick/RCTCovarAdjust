# in_dir <- "~/Documents/FileZilla/rct/run-10-04-22-2/"
# in_dir <- "~/Documents/FileZilla/rct/run-11-04-22-3/"
# in_dir <- "~/Documents/FileZilla/rct/run-08-05-22-2/"
in_dir <- "~/Documents/FileZilla/rct/run-18-04-22-1/"

df <- fread(paste0(in_dir, "summary.csv"))

BTYPE <- "med"

setnames(df, "bias_naive_med", "med_naive")
setnames(df, "bias_corr_med", "med_corr")

for(a in c("pocock", "obf")){
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

    addtorow <- list()
    addtorow$pos <- seq(3, nrow(df.sub2), by=3) %>% as.list
    addtorow$command <- rep("\\hline \n", length(addtorow$pos))

    tab <- xtable(df.sub2, align=rep("c", 7), digits=3,
                  caption=paste0("Simulation results for ", a,
                                 " boundaries and monitoring and final test statistic ", tn,
                                 ", estimating variance, 1000 simulations."))
    print(tab, include.rownames=FALSE,
          add.to.row = addtorow)
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

