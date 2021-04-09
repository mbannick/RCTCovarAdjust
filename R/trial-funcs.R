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

test.stat <- function(data, ancova=FALSE){
  mod <- fit.trial.data(data, ancova)
  return(mod$tstat)
}

t.anova <- function(data) test.stat(data, ancova=FALSE)
t.ancova <- function(data) test.stat(data, ancova=TRUE)

fit.stage <- function(data, stage, bounds, ancova=TRUE, ...){

  # Get covariate columns
  columns <- colnames(data)
  cov.cols <- columns[grepl("cov_", columns)]

  # Create list of test statistic functions
  t.stats <- list()
  for(i in 1:(stage-1)) t.stats[[i]] <- t.anova

  d.gens <- list()

  stage.df <- as.matrix(data[look == stage])

  if(ancova){

    # If doing ANCOVA at this stage, first estimate the
    # gamma and theta parameters and then simulate previous
    # trial data to be used to assess the boundaries.
    t.stats[[stage]] <- t.ancova
    model <- fit.trial.data(data=stage.df, ancova=TRUE)
    gamma <- model$bhat[cov.cols, ] %>% c
    theta <- apply(stage.df[, cov.cols], MARGIN=2, FUN=sd)
    sigma2 <- model$shat
    data.generator <- function(n) sim.ancova(n=n, delta=0, b0=0, theta=theta,
                                             gamma=gamma, sigma2=sigma2)

  } else {

    t.stats[[stage]] <- t.anova
    data.generator <- function(n) sim.ancova(n=n, delta=0, b0=0, theta=theta,
                                             gamma=rep(0, length(cov.cols)))
  }

  bounds <- get.boundaries.aspend(stat.func=t.stats,
                                  u_k=bounds,
                                  data.generator=data.generator, ...)

  return(bounds)
}


a.func.obf <- function(a, t) 4 * (1 - pnorm(qnorm(1-a/4)/sqrt(t)))
fit.stage(df, stage=5, bounds=c(2, 2, 2, 2), ancova=TRUE,
          a.func=a.func.obf, a=0.05,
          rates=t_k, N=1000, n_sims=1000)
