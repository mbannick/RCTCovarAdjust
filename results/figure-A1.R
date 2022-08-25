library(data.table)
library(ggplot2)
library(ggstar)
library(viridis)

# SET DIRECTORIES AND VERSION
args <- commandArgs(TRUE)
VERSION <- args[1]
IN_DIR <- "/Users/marlena/Documents/FileZilla/rct/"

# READ IN VERSION
df <- fread(paste0(IN_DIR, VERSION, "/summary.csv"))

# MAKE SURE THE VERSION HAS THE INFORMATION WE WANT
df <- df[(type == "inconsistent, corrected" & est_bounds == FALSE)]
df[, n := as.factor(n)]
df[afunc == "obf", bound_type := "OBF"]
df[afunc == "pocock", bound_type := "Pocock"]

# FIGURE A1 -----------------------------------

pdf("~/repos/Group-Sequential-Trials-Paper/figures/type1error-A1.pdf", height=8, width=10)
ann_text <- data.frame(rho=1-0.4, power=0.085,
                       lab="No anticipated reduction",
                       design_rho=1,
                       n=factor(50, levels=c(50, 100, 250, 1000)),
                       bound_type="OBF")
ann_text2 <- data.frame(rho=1-0.6, power=0.06,
                       lab="Anticipated and true\nreduction match",
                       design_rho=1,
                       n=factor(50, levels=c(50, 100, 250, 1000)),
                       bound_type="OBF")
extra_point <- data.frame(rho=1-0.35, power=0.06,
                          design_rho=1,
                          n=factor(50, levels=c(50, 100, 250, 1000)),
                          bound_type="OBF")
ggplot(data=df[delta == 0.0],
       aes(x=1-rho, y=power, color=1-design_rho, group=1-design_rho)) +
  facet_grid(bound_type~n) +
  geom_point(size=1, fill="white") +
  geom_hline(yintercept=0.05, linetype="dashed", color="black") +
  geom_line() +
  theme_minimal() +
  scale_color_viridis(option="D") +
  scale_fill_viridis(option="D", guide = 'none') +
  geom_star(data=df[rho == design_rho & delta == 0.0],
            size=3, aes(fill=1-design_rho)) +
  theme(legend.position="top") +
  labs(color="Anticipated Reduction in Variance",
       fill=NA,
       linetype="Boundary Type",
       y="Type I Error",
       x="Reduction in variance by using ANCOVA") +
  geom_text(data=ann_text, aes(label=lab), size=3) +
  geom_segment(
    x = 0.5, xend = 0.53, y = 0.083, yend = 0.0775,
    arrow = arrow(length = unit(5, "pt")),
    data=ann_text
  ) +
  geom_star(data=extra_point, size=2, fill="black") +
  geom_text(data=ann_text2, aes(label=lab), size=3, color="black")
dev.off()

