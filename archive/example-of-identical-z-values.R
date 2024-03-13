dg1 <- sim.data.closure(delta=0, n_cov=1, rho=0.2)
dg2 <- sim.data.closure(delta=0, n_cov=1, rho=0.4)

set.seed(0)
dat1 <- dg1(100)
set.seed(0)
dat2 <- dg2(100)

mod1 <- lm(dat1$y ~ 0 + dat1$X)
mod2 <- lm(dat2$y ~ 0 + dat2$X)

# Different estimate and standard error, but same z-value
summary(mod1)
summary(mod2)
