library(data.table)
library(ggplot2)

# SET DIRECTORIES AND VERSION
args <- commandArgs(TRUE)
# VERSION <- args[1]
VERSION <- "run-25-08-22-1"
IN_DIR <- "/Users/marlena/Documents/FileZilla/rct/"

# READ IN VERSION
df <- fread(paste0(IN_DIR, VERSION, "/figure1/summary.csv"))

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

df[label == "(B)(i)", type2 := "Regular"]
df[label == "(B)(ii)", type2 := "Uniform-Inflated"]
df[label == "(B)(iii)", type2 := "End-Inflated"]
df[label == "(A)(i)", type2 := "Regular (ANOVA)"]
df[label == "(A)(ii)", type2 := "Regular (ANCOVA)"]

df[, type2 := factor(type2, levels=c("Regular (ANOVA)", "Regular (ANCOVA)", "Regular", "Uniform-Inflated", "End-Inflated"))]

df[, n := as.factor(n)]
df[afunc == "obf", bound_type := "OBF"]
df[afunc == "pocock", bound_type := "Pocock"]
df[, bound_type := factor(bound_type, levels=c("Pocock", "OBF"))]

# FIGURE 4 -----------------------------------

# pdf("~/repos/Group-Sequential-Trials-Paper/figures/type1error-2.pdf", height=5, width=10)
plot1 <- ggplot(data=df[delta == 0.0 & stages == 3], aes(x=1-rho**2, y=power, color=label)) +
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
ggsave("~/OneDrive/Documents/2023-2024/GST-SIM-Revision/figures/Figure4v2.tiff",
       height=5, width=10, plot=plot1, dpi=600, units="in", bg = "white")

# dev.off()

ggplot(data=df[delta == 0.0 & stages == 3 & label %in% c("(B)(i)", "(B)(ii)", "(B)(iii)") & n == 1000],
       aes(x=1-rho**2, y=power, color=type2)) +
  facet_grid(bound_type ~ .) +
  geom_point() +
  geom_hline(yintercept=0.05, linetype="dashed", color="black") +
  geom_line() +
  theme_minimal() +
  scale_color_manual(values=c("#4DAF4A", "#984EA3", "#FF7F00")) +
  theme(legend.position="right") +
  labs(color="", linetype="Boundary Type",
       y="Type I Error Rate",
       x="Reduction in variance by using ANCOVA") +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(breaks=c(0.05, 0.055, 0.06))

# FIGURE 4 -----------------------------------

# pdf("~/repos/Group-Sequential-Trials-Paper/figures/power-2.pdf", height=5, width=10)
ggplot(data=df[delta == 0.1 & stages == 3 &
                          label %in% c("(A)(i)", "(A)(ii)", "(B)(ii)", "(B)(iii)") & n == 250],
                aes(x=1-rho**2, y=power, color=label)) +
  facet_grid(bound_type ~ n) +
  geom_point() + geom_line() +
  scale_color_brewer(palette="Set1") +
  theme_minimal() +
  theme(legend.position="top") +
  labs(color="Scenario", linetype="Boundary Type",
       x="Reduction in variance by using ANCOVA",
       y="Power") +
  xlim(c(0.4, 0.8)) +
  scale_x_continuous(labels = scales::percent)


ggsave("~/OneDrive/Documents/2023-2024/GST-SIM-Revision/figures/Figure5.tiff",
       height=5, width=10, plot=plot2, dpi=600, units="in", bg = "white")

# dev.off()
