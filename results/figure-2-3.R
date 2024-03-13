rm(list=ls())
library(latex2exp)
library(mvtnorm)
source("~/repos/RCTCovarAdjust/R/boundaries.R")

K <- 4 # number of stages total
c <- c(0, 0, 0, 1) # the change vector for changing to ANCOVA at the last stage

# FIGURE 2 ----------------------------------

pocock_10_vec <- get.boundaries.design(obf=FALSE, rho=1, rates=(1:K/K), change=c)
pocock_0K_vec <- get.boundaries.design(obf=FALSE, rho=0.5, rates=(1:K/K), change=c)

obf_10_vec <- get.boundaries.design(obf=TRUE, rho=1, rates=(1:K/K), change=c)
obf_0K_vec <- get.boundaries.design(obf=TRUE, rho=0.5, rates=(1:K/K), change=c)

tiff("~/OneDrive/Documents/2023-2024/GST-SIM-Revision/figures/Figure2.tiff", height=5, width=9, units="in", res=600)
par(mfrow=c(1, 2))
plot(1:K, pocock_10_vec, type='l', ylim=c(2.0, 4.20),
     xaxt="n", xlab="Stage", ylab="Critical Value",
     main="Pocock Bounds",
     font.main=1)
points(1:K, pocock_10_vec, pch=16)
lines(1:K, pocock_0K_vec, col='blue')
points(1:K, pocock_0K_vec, pch=16, col='blue')
axis(1, at = 1:K)
legend(x=3, y=4.1, legend=c(1.0, 0.5),
       lty=c(1, 1), col=c("black", "blue"),
       title=expression(rho),
       bty="n",
       cex=0.9,
       pch=16)

plot(1:K, obf_10_vec, type='l', ylim=c(2.0, 4.20),
     xaxt="n", xlab="Stage", ylab="Critical Value",
     main="O'Brien-Fleming Bounds",
     font.main=1)
points(1:K, obf_10_vec, pch=16)
lines(1:K, obf_0K_vec, col='blue')
points(1:K, obf_0K_vec, pch=16, col='blue')
axis(1, at = 1:K)
legend(x=3, y=4.1, legend=c(1, 0.5),
       lty=c(1, 1), col=c("black", "blue"),
       title=expression(rho),
       bty="n",
       cex=0.9,
       pch=c(16, 16))
dev.off()

# FIGURE 3 ----------------------------------

pocock_10_alpha <- get.boundaries.aspend(a.func=pocock.spend(0.05), rates=(1:K)/K, rho=1, change=c)
obf_10_alpha <- get.boundaries.aspend(a.func=obf.spend(0.05), rates=(1:K)/K, rho=1, change=c)

pocock_05_alpha <- get.boundaries.aspend(a.func=pocock.spend(0.05), rates=(1:K)/K, rho=0.5, change=c)
obf_05_alpha <- get.boundaries.aspend(a.func=obf.spend(0.05), rates=(1:K)/K, rho=0.5, change=c)

tiff("~/OneDrive/Documents/2023-2024/GST-SIM-Revision/figures/Figure3.tiff", height=5, width=9, units="in", res=600)
par(mfrow=c(1, 2))
plot(1:K, pocock_10_alpha[, 2], type='l', ylim=c(2.0, 4.30),
     xaxt="n", xlab="Stage", ylab="Critical Value",
     main="Pocock Approximation",
     font.main=1)
points(1:K, pocock_10_alpha[, 2], pch=16)
lines(1:K, pocock_05_alpha[, 2], col='red')
points(1:K, pocock_05_alpha[, 2], pch=16, col='red')
axis(1, at = 1:K)
legend(x=3, y=4.1, legend=c(1.0, 0.5),
       lty=c(1, 1), col=c("black", "red"),
       title=expression(rho),
       bty="n",
       cex=0.9,
       pch=16)

plot(1:K, obf_10_alpha[, 2], type='l', ylim=c(2.0, 4.30),
     xaxt="n", xlab="Stage", ylab="Critical Value",
     main="O'Brien-Fleming Approximation",
     font.main=1)
points(1:K, obf_10_alpha[, 2], pch=16)
lines(1:K, obf_05_alpha[, 2], col='red')
points(1:K, obf_05_alpha[, 2], pch=16, col='red')
axis(1, at = 1:K)
legend(x=3, y=4.1, legend=c(1.0, 0.5),
       lty=c(1, 1), col=c("black", "red"),
       title=expression(rho),
       bty="n",
       cex=0.9,
       pch=c(16, 16))
dev.off()
