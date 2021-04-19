library(magrittr)
library(MASS)
library(data.table)

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

ancova.stats <- function(n_sims, n_k, gamma, theta, sigma2, ancova=TRUE){
  K <- length(n_k)
  co.var <- sqrt(sigma2) / sqrt(sigma2 + (gamma %*% theta)**2)
  Sigma <- basic.cov(n_k)
  if(ancova){
    Sigma[K, 1:(K-1)] <- Sigma[K, 1:(K-1)] * c(co.var)
    Sigma[1:(K-1), K] <- Sigma[1:(K-1), K] * c(co.var)
  }
  sims <- mvrnorm(n=n_sims, mu=rep(0, length(n_k)),
                  Sigma=Sigma)
  return(sims)
}

fit.stage <- function(data, stage, bounds, ancova=TRUE, estimate_sigma=TRUE, ...){
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
    gamma <- model$bhat[cov.cols, ]
    theta <- apply(matrix(stage.df[, cov.cols]), MARGIN=2, FUN=sd)
    if(estimate_sigma){
      sigma2 <- model$shat
    } else {
      sigma2 <- 1
    }

    sim.generator <- function(n_sims, n_k) ancova.stats(n_sims=n_sims, n_k=n_k,
                                                        gamma=gamma,
                                                        theta=theta,
                                                        sigma2=sigma2,
                                                        ancova=TRUE)

  } else {

    model <- fit.trial.data(data=stage.df, ancova=FALSE)
    if(estimate_sigma){
      sigma2 <- model$shat
    } else {
      sigma2 <- 1
    }
    sim.generator <- function(n_sims, n_k) ancova.stats(n_sims=n_sims, n_k=n_k,
                                                        gamma=0,
                                                        theta=0,
                                                        sigma2=sigma2,
                                                        ancova=FALSE)
  }

  bounds <- get.boundaries.aspend(u_k=bounds,
                                  sim.generator=sim.generator, ...)

  return(bounds)
}

