# HACK

df.sub <- fread("summary.csv")
df.sub <- df.sub[n < 2000]

setorder(df.sub, delta, n)

df.sub2 <- df.sub[, .(delta, n, bias_naive, bias_corr, bias_naive_med, bias_corr_med, cover_naive, cover_corr)]
df.sub2 <- data.table(df.sub2)

tab <- xtable(df.sub2, align=rep("c", 9), digits=3)
print(tab, include.rownames=FALSE)
