library(magrittr)
library(data.table)

t_k <- 1:5/5
p_k <- (1:5/5)**2

df <- sim.ancova.partial(t_k, p_k, N=100,
                         delta=1, theta=1, gamma=c(0.5, 2),
                         sigma2=1, b0=0)

fit.trial.data <- function(data, ancova=FALSE){

  data <- as.matrix(data)
  columns <- colnames(data)

  out.col <- which(columns == "y")
  t.col <- which(columns == "t")
  var.cols <- which(columns == "int")

  if(ancova){
    cov.cols <- which(grepl("cov_", columns))
    var.cols <- c(var.cols, cov.cols)
  }

  mod <- fit.regression(data=data, out.col=out.col,
                        t.col=t.col, var.cols=var.cols)
  return(mod)
}

ancova.stats <- function(n_sims, n_k, ){

  Sigma <- basic.cov(n_k)
  sims <- mvrnorm(n=n_sims, mu=rep(0, length(n_k)),
                   Sigma=basic.cov(n_k))

  return(sims)
}

fit.stage <- function(data, stage, bounds, ancova=TRUE, ...){

  # Get covariate columns
  columns <- colnames(data)
  cov.cols <- columns[grepl("cov_", columns)]

  # Get the covariance matrix for the sim.generator function
  stage.df <- as.matrix(data[look == stage])

  if(ancova){

    # If doing ANCOVA at this stage, first estimate the
    # gamma and theta parameters and then simulate previous
    # trial data to be used to assess the boundaries
    model <- fit.trial.data(data=stage.df, ancova=TRUE)
    gamma <- model$bhat[cov.cols, ] %>% c
    theta <- apply(stage.df[, cov.cols], MARGIN=2, FUN=sd)
    sigma2 <- model$shat

    # data.generator <- function(n) sim.ancova(n=n, delta=0, b0=0, theta=theta,
    #                                          gamma=gamma, sigma2=sigma2)
    sim.generator <- ancova.stats

  } else {

    # data.generator <- function(n) sim.ancova(n=n, delta=0, b0=0, theta=theta,
    #                                          gamma=rep(0, length(cov.cols)))
    sim.generator <- ancova.stats
  }

  bounds <- get.boundaries.aspend(u_k=bounds,
                                  sim.generator=sim.generator, ...)

  return(bounds)
}

a.func.obf <- function(a, t) 4 * (1 - pnorm(qnorm(1-a/4)/sqrt(t)))
fit.stage(df, stage=5, bounds=c(2, 2, 2, 2), ancova=TRUE,
          a.func=a.func.obf, a=0.05,
          rates=t_k, N=1000, n_sims=1000)

bounds <- c()
for(i in 1:5){

  ancova <- i == 5
  bounds <- fit.stage(df, stage=5, bounds=bounds, ancova=ancova,
                      a.func=a.func.obf, a=0.05,
                      rates=t_k, N=1000, n_sims=100000)
}

bounds
