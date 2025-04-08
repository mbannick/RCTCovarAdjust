
Sigma <- matrix(c(1, 1/sqrt(2), 1/sqrt(2), 1), nrow=2, byrow=TRUE)

u <- qnorm(1-0.05/4)
p1 <- (1-pnorm(u))*2
p2 <- pmvnorm(lower=c(-u, u), upper=c(u, Inf), mean=c(0, 0), sigma=Sigma)
p3 <- pmvnorm(lower=c(-u, -Inf), upper=c(u, -u), mean=c(0, 0), sigma=Sigma)

p1 + p2 + p3
