library(data.table)
library(ggplot2)

df <- fread("/Users/marlena/Documents/FileZilla/rct/run-24-02-22-2/summary.csv")
df <- fread("/Users/marlena/Documents/FileZilla/rct/run-10-04-22-2/summary.csv")

df[, type := factor(type, levels=c("un-adjusted", "adjusted",
                                   "inconsistent, naive", "inconsistent, corrected"))]

df[, n := as.factor(n)]
df[, std := sqrt(power * (1 - power)) / sqrt(10000)]
df[, lower := power - qnorm(0.975) * std]
df[, upper := power + qnorm(0.975) * std]

df[, "un-adjusted"]

# Type I Error Plot

pdf("~/OneDrive/Documents/2021-2022/BIOST 600-noah/type1error-2.pdf", height=8, width=14)
ggplot(data=df[delta == 0.0 & stages == 3], aes(x=1-rho, y=power, color=type, linetype=afunc)) +
  facet_grid(afunc ~ n) +
  geom_point() +
  geom_hline(yintercept=0.05, linetype="dashed", color="black") +
  geom_line() +
  theme_minimal() +
  scale_color_brewer(palette="Set1") +
  theme(legend.position="top")
dev.off()

ggplot(data=df[delta > 0 & stages == 3], aes(x=rho, y=power, color=n, linetype=afunc)) +
  facet_grid(type ~ delta) +
  geom_point() + geom_line()

pdf("~/OneDrive/Documents/2021-2022/BIOST 600-noah/power.pdf", height=8, width=14)
ggplot(data=df[delta == 0.1 & stages == 3], aes(x=1-rho, y=power, color=type, linetype=afunc)) +
  facet_grid(afunc ~ n) +
  geom_point() + geom_line() +
  scale_color_brewer(palette="Set1") +
  theme_minimal() +
  theme(legend.position="top")
dev.off()
