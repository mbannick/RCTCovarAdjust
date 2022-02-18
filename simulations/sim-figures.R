# Type I Error Plot

pdf("~/OneDrive/Documents/2021-2022/BIOST 600-noah/type1error.pdf", height=8, width=14)
ggplot(data=df[delta == 0.0 & stages == 3], aes(x=rho, y=power, color=type, linetype=afunc)) +
  facet_grid(~ n) +
  geom_point() +
  geom_hline(yintercept=0.05, linetype="dashed", color="black") +
  geom_line() +
  theme_minimal() +
  scale_color_brewer(palette="Set1") +
  theme(legend.position="top")
# geom_errorbar(aes(ymin=lower, ymax=upper))
dev.off()

ggplot(data=df[delta > 0 & stages == 3], aes(x=rho, y=power, color=n, linetype=afunc)) +
  facet_grid(type ~ delta) +
  geom_point() + geom_line()

pdf("~/OneDrive/Documents/2021-2022/BIOST 600-noah/power.pdf", height=8, width=14)
ggplot(data=df[delta == 0.1 & stages == 3], aes(x=rho, y=power, color=type, linetype=afunc)) +
  facet_grid( ~ n) +
  geom_point() + geom_line() +
  scale_color_brewer(palette="Set1") +
  theme_minimal() +
  theme(legend.position="top")
dev.off()
