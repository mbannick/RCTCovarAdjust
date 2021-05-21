source("R/recursive.R")
library(MASS)

nsim <- 1000000
lower = -15
upper = 15
l.out = 1000000
grid_x <- seq(from = lower, to = upper, length.out = l.out)
delta = (upper - lower)/l.out
RHO <- 0

a <- 0.05
as <- a * log(1 + (exp(1) - 1)*c(0.5, 1))
as <- c(as[1], diff(as))

true <- c(2.157, 2.201)

corr.1 <- matrix(c(1, 1/sqrt(2), 1/sqrt(2), 1), nrow=2)
corr.1[2, 1] <- corr.1[2, 1] * RHO
corr.1[1, 2] <- corr.1[1, 2] * RHO
t1.sim <- rnorm(n=nsim, mean=0, sd=1)
u1.sim <- uniroot(function(u) mean(abs(t1.sim) >= u) - as[1],
                  c(0, 100))$root
t2.sim <- mvrnorm(n=nsim, mu=c(0, 0), Sigma=corr.1)
u2.sim <- uniroot(function(u) mean((abs(t2.sim[, 2]) >= u) & (abs(t2.sim[, 1]) < u1.sim)) - as[2],
                  c(0, 100))$root

f1 <- dnorm(grid_x)
u1.int <- uniroot(function(u) sum(f1[abs(grid_x) >= u]*delta) - as[1],
                  c(0, 100))$root
u1.int

# gridded_base_dist <- dnorm(grid_x)
# gridded_new_increment <- dnorm(grid_x)
#
# recursive_integ <- function(gridded_base_dist,
#                             gridded_new_increment,
#                             delta,
#                             cutoff){
#   gridded_base_dist_continue <- gridded_base_dist
#   gridded_base_dist_continue[which(abs(grid_x) >= cutoff)] <- 0
#   gridded_base_dist_continue <- gridded_base_dist_continue # / (sum(gridded_base_dist_continue)*delta)
#
#   est_flipped <- convolve(gridded_base_dist_continue,
#                           gridded_new_increment)  * delta
#   gridded_convolved_dist <- est_flipped[c(ceiling(length(grid_x)/2):length(grid_x),
#                                           1:(ceiling(length(grid_x)/2)-1))]
#   return(gridded_convolved_dist)
# }
# approx_dist <- recursive_integ(gridded_base_dist,
#                                gridded_new_increment,
#                                delta,
#                                u1.int)
#
# u2.int <- uniroot(function(u) sum(approx_dist[which(abs(grid_x) >= u)])*delta - as[2],
#                   c(0, 100))$root
# u2.int <- u2.int / sqrt(2)

f1.cont <- f1
f1.cont[which(abs(grid_x) >= u1.int)] <- 0

# calculate rho
v <- 1 + 2 * (1 - RHO) * 1

f2.inc <- dnorm(grid_x, sd=sqrt(v))
f2 <- convolve(f1.cont, f2.inc) * delta
f2 <- f2[c(ceiling(length(grid_x)/2):length(grid_x), 1:(ceiling(length(grid_x)/2)-1))]

u2.int <- uniroot(function(u) sum(f2[which(abs(grid_x) >= u)]*delta) - as[2],
                  c(0, 100))$root
u2.int / sqrt(2)



# dens <- get.joint.density(n_k=c(100, 100), lr_bounds=c(u1.int),
#                           gridsize=l.out, rho=RHO, lower=lower, upper=upper)

# f2 <- dens$dens[, 2]

par(mfrow=c(1, 2))
hist(t2.sim[, 2], prob=TRUE, freq=FALSE, breaks=100, col='green', add=F)
hist(t2.sim[, 2][which(t2.sim[, 1] < u1.sim)], prob=TRUE, freq=FALSE, breaks=100, col='pink', add=F)
