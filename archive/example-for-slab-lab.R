library(latex2exp)
source("~/repos/RCTCovarAdjust/R/covariance.R")

K <- 4
pocock_10_vec <- get.boundaries.design(obf=FALSE, rho=1, rates=(1:K/K))
obf_10_vec <- get.boundaries.design(obf=TRUE, rho=1, rates=(1:K/K))

pdf("~/repos/SLAB-Lab-WIN-2022/bounds.pdf", height=5, width=9)

par(mfrow=c(1, 2))
plot(1:K, pocock_10_vec, type='l', ylim=c(-4.15, 4.15),
     font.main=1, ylab="Critical Value", xlab="Stage",
     main="Pocock Bounds", xaxt="n")
points(1:K, pocock_10_vec, pch=16)
lines(1:K, -pocock_10_vec)
points(1:K, -pocock_10_vec, pch=16)
lines(1:K, rep(0, 4), lty=3)
axis(1, at = 1:K)

plot(1:K, obf_10_vec, type='l', ylim=c(-4.15, 4.15),
     font.main=1, ylab="Critical Value", xlab="Stage",
     main="O'Brien-Fleming Bounds", xaxt="n")
points(1:K, obf_10_vec, pch=16)
lines(1:K, -obf_10_vec)
points(1:K, -obf_10_vec, pch=16)
lines(1:K, rep(0, 4), lty=3)
axis(1, at=1:K)

dev.off()

t <- seq(0, 1, by=0.01)

pdf("~/repos/SLAB-Lab-WIN-2022/bounds-alpha.pdf", height=5, width=5)
plot(t, obf.spend(0.05)(t), type='l', main="Alpha-Spending Functions",
     font.main=1, ylab=TeX("$\\alpha^*(t)$"), xlab="Information Accrued, t")
lines(t, pocock.spend(0.05)(t), type='l', col='purple')
legend(x=0, y=0.05, legend=c("Pocock", "OBF"),
       lty=c(1, 1), col=c("purple", "black"),
       bty="n",
       cex=0.9)
dev.off()
