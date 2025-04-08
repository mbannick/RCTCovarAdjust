bounds <- result$u.bounds

numbounds <- function(vec) sum(vec != 0)
stages <- apply(bounds, MARGIN=1, FUN=numbounds)

proportions(table(stages))

df <- data.table(
  adj_estimate=result$point,
  estimate=result$est,
  final_tstat=result$tstat,
  stage=stages
)
df[, sim := .I]
melt(df, id.vars=c("sim", "stage", "final_tstat"))

u_k <- rep(2.178273, 2)
n_k <- c(50, 100)
vars <- u_k / sqrt(n_k)

stagedf <- data.frame(
  stage=1:2,
  bound=vars
)
df <- merge(df, stagedf, by=c("stage"))
df[, stagename := paste0("Stage ", stage)]

library(ggplot2)
pdf("stages-4.pdf", height=3, width=7)
ggplot(df) + geom_histogram(aes(x=estimate), bins=50, fill='#348feb') +
  facet_wrap(~stagename, nrow=1) +
  geom_vline(xintercept=0.15, linetype='dashed') +
  geom_vline(aes(xintercept=bound), color="#f73942", linetype='dashed') +
  geom_vline(aes(xintercept=-bound), color="#f73942", linetype='dashed') +
  labs(x="Point Estimate")
dev.off()

pdf("dist-of-test-stat-blank.pdf", height=5, width=10)
par(mfrow=c(1, 2))
vals <- seq(-4, 4, by=0.01)
plot(vals, dnorm(vals), lty='dashed', type='l',
      main="Distribution of Test Statistics Under the Null",
     xlab="Test Statistic (10,000 simulations)", ylab="Sampling Density",
     ylim=c(0, 0.4))
vals2 <- seq(0, 1, by=0.01)
plot(vals2, dunif(vals2), col='black', lty='dashed', type='l',
      main="Distribution of P-values Under the Null", ylab="Sampling Density",
      xlab="P-value (10,000 simulations)",
     ylim=c(0, 2.0))
dev.off()

pdf("dist-of-test-stat.pdf", height=5, width=10)
par(mfrow=c(1, 2))
vals <- seq(-4, 4, by=0.01)
d <- density(result$tstat)
plot(d, main="Distribution of Test Statistics Under the Null",
     xlab="Test Statistic (10,000 simulations)", ylab="Sampling Density",
     ylim=c(0, 0.4))
lines(vals, dnorm(vals), col='black', lty='dashed')
abline(v=2.178273, col='red', lty='dashed')
abline(v=-2.178273, col='red', lty='dashed')
vals2 <- seq(0, 1, by=0.01)
hist(pnorm(abs(result$tstat), lower.tail=F)*2, breaks=40,
     main="Distribution of P-values Under the Null",
     xlab="P-value (10,000 simulations)", freq=F, ylab="Sampling Density",
     ylim=c(0, 2.0))
lines(vals2, dunif(vals2), col='black', lty='dashed')
dev.off()
# hist(result$tstat, breaks=50, freq=F)
vals <- seq(-4, 4, by=0.01)
lines(vals, dnorm(vals), col='red')

hist(result$tstat, breaks=75)
hist(pnorm(abs(result$tstat), lower.tail=F)*2, breaks=100)


lines(seq(-3, 3, by=0.01), dnorm(seq(-3, 3, by=0.01)))

tstat <- result$tstat
hist(result$pval, breaks=30)
