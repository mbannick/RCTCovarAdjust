library(data.table)
library(ggplot2)

# SET DIRECTORIES AND VERSION
args <- commandArgs(TRUE)
VERSION <- args[1]
IN_DIR <- "/Users/marlena/Documents/FileZilla/rct/"

# READ IN VERSION
df <- fread(paste0(IN_DIR, VERSION, "/summary.csv"))

# MAKE SURE THE VERSION HAS THE INFORMATION WE WANT
df <- df[((type %in% c("un-adjusted", "adjusted",
           "inconsistent, naive", "inconsistent, corrected")) & est_bounds == TRUE) |
           (type == "inconsistent, corrected" & est_bounds == FALSE)]

# LABEL BASED ON THE NUMBERING FOR THE PAPER
df[type == "un-adjusted", label := "(A)(i)"]
df[type == "adjusted", label := "(A)(ii)"]
df[type == "inconsistent, naive", label := "(B)(i)"]
df[type == "inconsistent, corrected" & est_bounds == FALSE, label := "(B)(ii)"]
df[type == "inconsistent, corrected" & est_bounds == TRUE, label := "(B)(iii)"]

df[, n := as.factor(n)]
df[afunc == "obf", bound_type := "OBF"]
df[afunc == "pocock", bound_type := "Pocock"]

# FIGURE 4 -----------------------------------

pdf("~/repos/Group-Sequential-Trials-Paper/figures/type1error-2.pdf", height=5, width=10)
ggplot(data=df[delta == 0.0 & stages == 3], aes(x=1-rho**2, y=power, color=label)) +
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
dev.off()

# FIGURE 4 -----------------------------------

pdf("~/repos/Group-Sequential-Trials-Paper/figures/power-2.pdf", height=5, width=10)
ggplot(data=df[delta == 0.1 & stages == 3], aes(x=1-rho**2, y=power, color=label)) +
  facet_grid(bound_type ~ n) +
  geom_point() + geom_line() +
  scale_color_brewer(palette="Set1") +
  theme_minimal() +
  theme(legend.position="top") +
  labs(color="Scenario", linetype="Boundary Type",
       x="Reduction in variance by using ANCOVA",
       y="Power")
dev.off()
