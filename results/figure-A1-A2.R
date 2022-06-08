library(data.table)
library(ggplot2)
library(ggstar)

# SET DIRECTORIES AND VERSION
args <- commandArgs(TRUE)
VERSION <- "run-07-06-22-8"
IN_DIR <- "/Users/marlena/Documents/FileZilla/rct/"

# READ IN VERSION
df <- fread(paste0(IN_DIR, VERSION, "/summary.csv"))

# MAKE SURE THE VERSION HAS THE INFORMATION WE WANT
df <- df[(type == "inconsistent, corrected" & est_bounds == FALSE)]
df[, n := as.factor(n)]
df[afunc == "obf", bound_type := "OBF"]
df[afunc == "pocock", bound_type := "Pocock"]

# FIGURE A1 -----------------------------------

pdf("~/repos/Group-Sequential-Trials-Paper/figures/type1error-A1.pdf", height=5, width=10)
ggplot(data=df[delta == 0.0], aes(x=1-rho, y=power, color=1-design_rho, group=1-design_rho)) +
  facet_grid(bound_type~n) +
  geom_point() +
  geom_hline(yintercept=0.05, linetype="dashed", color="black") +
  geom_line() +
  theme_minimal() +
  scale_color_continuous() +
  geom_star(data=df[rho == design_rho & delta == 0.0], size=3, aes(fill=1-design_rho)) +
  # scale_color_brewer(palette="Set1") +
  theme(legend.position="top") +
  labs(color="Anticipated Reduction in Variance",
       fill="Anticipated Reduction in Variance",
       linetype="Boundary Type",
       y="Type I Error",
       x="Anticipated reduction in variance by using ANCOVA")
dev.off()

# FIGURE A2 -----------------------------------

pdf("~/repos/Group-Sequential-Trials-Paper/figures/power-A2.pdf", height=5, width=10)
ggplot(data=df[delta == 0.1 & n == 250], aes(x=1-design_rho, y=power, group=design_rho)) +
  facet_grid(~bound_type) +
  geom_point() + geom_line() +
  scale_color_brewer(palette="Set1") +
  theme_minimal() +
  theme(legend.position="top") +
  labs(color="Scenario", linetype="Boundary Type",
       x="Anticipated reduction in variance by using ANCOVA",
       y="Power")
dev.off()
