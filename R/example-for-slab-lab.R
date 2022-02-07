source("~/repos/RCTCovarAdjust/R/covariance.R")


K <- 2
cr <- qnorm(1-0.05/(K*2))
1 - pmvnorm(lower=-rep(cr, K),
            upper=rep(cr, K),
            mean=rep(0, K), corr=diag(K),
            algorithm=Miwa(steps=1000))

get.power <- function(c, K, obf=FALSE, rho=1){

  rates <- (1:K)/K
  u_k <- rep(c, K)

  if(obf) u_k <- u_k / sqrt(1:K)

  Sigma <- corr.mat(rates, rho=rho)
  print(Sigma)
  val <- 1 - pmvnorm(lower=-u_k,
                     upper=u_k,
                     mean=rep(0, K), corr=Sigma,
                     algorithm=Miwa(steps=1000))

  return(val)
}

get.power(2.361286, K=4, rho=0.05)
get.power(2.2895, K=3, rho=0.5)
get.power(2.178, K=2, rho=0.1)
get.power(4.048565, K=4, rho=0.001, obf=T)

get.bound <- function(K, obf=FALSE, rho=1, power=0.05){

  f <- function(x) get.power(x, K=K, obf=obf, rho=rho) - power
  s <- uniroot(f, interval=c(0, 100))
  return(s$root)
}

get.bound(3)
sapply(2:5, get.bound)
sapply(2:5, get.bound, rho=0.5)

sapply(2:5, get.bound, obf=TRUE)
sapply(2:5, get.bound, obf=TRUE, rho=0.5)

K <- 4
rho <- 0.01
pocock_10 <- get.bound(K=K, rho=1, obf=FALSE)
pocock_0K <- get.bound(K=K, rho=rho, obf=FALSE)

obf_10 <- get.bound(K=K, rho=1, obf=TRUE)
obf_0K <- get.bound(K=K, rho=rho, obf=TRUE)

pocock_10_vec <- rep(pocock_10, K)
pocock_0K_vec <- rep(pocock_0K, K)

obf_10_vec <- obf_10 / sqrt(1:K)
obf_0K_vec <- obf_0K / sqrt(1:K)

pocock_10_alpha <- get.boundaries(a.func=pocock.spend(0.05), rates=(1:K)/K, rho=1)
obf_10_alpha <- get.boundaries(a.func=obf.spend(0.05), rates=(1:K)/K, rho=1)

pocock_05_alpha <- get.boundaries(a.func=pocock.spend(0.05), rates=(1:K)/K, rho=0.5)
obf_05_alpha <- get.boundaries(a.func=obf.spend(0.05), rates=t, rho=0.5)

pdf("~/repos/SLAB-Lab-WIN-2022/bounds-adjusted.pdf", height=5, width=9)
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
legend(x=3, y=4.1, legend=c(1.0, 0.5),
       lty=c(1, 1), col=c("black", "blue"),
       title=expression(rho),
       bty="n",
       cex=0.9,
       pch=c(16, 16))
dev.off()


pdf("~/repos/SLAB-Lab-WIN-2022/bounds-adjusted-alpha-1.pdf", height=5, width=9)
par(mfrow=c(1, 2))
plot(1:K, pocock_10_alpha, type='l', ylim=c(2.0, 4.30),
     xaxt="n", xlab="Stage", ylab="Critical Value",
     main="Pocock Approximation",
     font.main=1)
points(1:K, pocock_10_alpha, pch=16)
axis(1, at = 1:K)
legend(x=3, y=4.1, legend=c(1.0),
       lty=c(1), col=c("black"),
       title=expression(rho),
       bty="n",
       cex=0.9,
       pch=16)

plot(1:K, obf_10_alpha, type='l', ylim=c(2.0, 4.30),
     xaxt="n", xlab="Stage", ylab="Critical Value",
     main="O'Brien-Fleming Approximation",
     font.main=1)
points(1:K, obf_10_alpha, pch=16)
axis(1, at = 1:K)
legend(x=3, y=4.1, legend=c(1.0),
       lty=c(1), col=c("black"),
       title=expression(rho),
       bty="n",
       cex=0.9,
       pch=c(16))
dev.off()

pdf("~/repos/SLAB-Lab-WIN-2022/bounds-adjusted-alpha-2.pdf", height=5, width=9)
par(mfrow=c(1, 2))
plot(1:K, pocock_10_alpha, type='l', ylim=c(2.0, 4.30),
     xaxt="n", xlab="Stage", ylab="Critical Value",
     main="Pocock Approximation",
     font.main=1)
points(1:K, pocock_10_alpha, pch=16)
lines(1:K, pocock_05_alpha, col='red')
points(1:K, pocock_05_alpha, pch=16, col='red')
axis(1, at = 1:K)
legend(x=3, y=4.1, legend=c(1.0, 0.5),
       lty=c(1, 1), col=c("black", "red"),
       title=expression(rho),
       bty="n",
       cex=0.9,
       pch=16)

plot(1:K, obf_10_alpha, type='l', ylim=c(2.0, 4.30),
     xaxt="n", xlab="Stage", ylab="Critical Value",
     main="O'Brien-Fleming Approximation",
     font.main=1)
points(1:K, obf_10_alpha, pch=16)
lines(1:K, obf_05_alpha, col='red')
points(1:K, obf_05_alpha, pch=16, col='red')
axis(1, at = 1:K)
legend(x=3, y=4.1, legend=c(1.0, 0.5),
       lty=c(1, 1), col=c("black", "red"),
       title=expression(rho),
       bty="n",
       cex=0.9,
       pch=c(16, 16))
dev.off()


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

library(latex2exp)
pdf("~/repos/SLAB-Lab-WIN-2022/bounds-alpha.pdf", height=5, width=5)
plot(t, obf.spend(0.05)(t), type='l', main="Alpha-Spending Functions",
     font.main=1, ylab=TeX("$\\alpha^*(t)$"), xlab="Information Accrued, t")
lines(t, pocock.spend(0.05)(t), type='l', col='purple')
legend(x=0, y=0.05, legend=c("Pocock", "OBF"),
       lty=c(1, 1), col=c("purple", "black"),
       bty="n",
       cex=0.9)
dev.off()
