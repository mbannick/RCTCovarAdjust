# HACK

df.sub <- fread("summary.csv")
df.sub <- df.sub[n < 2000]

setorder(df.sub, delta, n)

df.sub2 <- df.sub[, .(delta, n, bias_naive, bias_corr, bias_naive_med, bias_corr_med, cover_naive, cover_corr)]
df.sub2 <- data.table(df.sub2)

# addtorow <- list()
# addtorow$pos <- seq(5, nrow(df.sub2), by=5) %>% as.list
# addtorow$command <- rep("\\hline \n", length(addtorow$pos))

tab <- xtable(df.sub2, align=rep("c", 9), digits=3)
print(tab, include.rownames=FALSE)
# add.to.row = addtorow)
  }
}
