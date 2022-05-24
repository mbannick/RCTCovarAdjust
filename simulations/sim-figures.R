library(data.table)
library(ggplot2)

# df <- fread("/Users/marlena/Documents/FileZilla/rct/run-24-02-22-2/summary.csv")
# df0 <- fread("/Users/marlena/Documents/FileZilla/rct/run-10-04-22-2/summary.csv")
df <- fread("/Users/marlena/Documents/FileZilla/rct/run-16-05-22-5/summary.csv")

df <- df[((type %in% c("un-adjusted", "adjusted",
           "inconsistent, naive", "inconsistent, corrected")) & est_bounds == TRUE) |
           (type == "inconsistent, corrected" & est_bounds == FALSE)]

df[type == "un-adjusted", label := "(A)(i)"]
df[type == "adjusted", label := "(A)(ii)"]
df[type == "inconsistent, naive", label := "(B)(i)"]
df[type == "inconsistent, corrected" & est_bounds == FALSE, label := "(B)(ii)"]
df[type == "inconsistent, corrected" & est_bounds == TRUE, label := "(B)(iii)"]

df[, n := as.factor(n)]
df[, std := sqrt(power * (1 - power)) / sqrt(10000)]
df[, lower := power - qnorm(0.975) * std]
df[, upper := power + qnorm(0.975) * std]

df[afunc == "obf", bound_type := "OBF"]
df[afunc == "pocock", bound_type := "Pocock"]

# Type I Error Plot

pdf("~/repos/Group-Sequential-Trials-Paper/figures/type1error-2.pdf", height=5, width=10)
ggplot(data=df[delta == 0.0 & stages == 3], aes(x=1-rho, y=power, color=label)) +
  facet_grid(bound_type ~ n) +
  geom_point() +
  geom_hline(yintercept=0.05, linetype="dashed", color="black") +
  geom_line() +
  theme_minimal() +
  scale_color_brewer(palette="Set1") +
  theme(legend.position="top") +
  labs(color="Scenario", linetype="Boundary Type",
       y="Type I Error",
       x="Reduction in variance by using ANCOVA")
  # geom_errorbar(aes(ymin=lower, ymax=upper, x=1-rho),
  #               linetype='dashed', width=.1)
dev.off()

# ggplot(data=df[delta > 0 & stages == 3], aes(x=rho, y=power, color=n, linetype=afunc)) +
#   facet_grid(type ~ delta) +
#   geom_point() + geom_line()

pdf("~/repos/Group-Sequential-Trials-Paper/figures/power-2.pdf", height=5, width=10)
ggplot(data=df[delta == 0.1 & stages == 3], aes(x=1-rho, y=power, color=label)) +
  facet_grid(bound_type ~ n) +
  geom_point() + geom_line() +
  scale_color_brewer(palette="Set1") +
  theme_minimal() +
  theme(legend.position="top") +
  labs(color="Scenario", linetype="Boundary Type",
       x="Reduction in variance by using ANCOVA",
       y="Power")
dev.off()
