library(mvtnorm)

solve.boundary3 <- function(power, corr, u_k, ...){
  if(is.null(u_k)){
    bin <- function(x) 2 * pnorm(x, mean=0, sd=1, lower.tail=F) - power
  } else {
    bin <- function(x) 1 - pmvnorm(lower=c(-u_k, -x),
                                   upper=c(u_k, x),
                                   mean=rep(0, length(u_k)+1), corr=corr) - power
  }
  u <- uniroot(bin, c(0, 100))
  return(u$root)
}

get.boundaries.aspend3 <- function(a.func, a, rates, N,
                                   u_k=c(), rho=1, ...){

  # Get sample size increments
  K <- length(rates)
  n_k <- round(rates*N)
  n_k <- c(n_k[1], diff(n_k))

  # Number of *fixed* previous bounds
  K_prev <- length(u_k)
  tau <- get.tau(n_k)
  tau_s <- cumsum(tau)

  # Create alpha spending function
  a.spend <- function(t) a.func(a=a, t)

  # Append 0 onto the rates
  a.cuml <- a.spend(rates)

  bounds <- c()
  for(i in 1:K){
    if(i > K_prev){

      # Create covariance matrix
      Sigma <- basic.cov(n_k[1:i])
      if(i == K){
        Sigma[K, 1:(K-1)] <- Sigma[K, 1:(K-1)] * rho
        Sigma[1:(K-1), K] <- Sigma[1:(K-1), K] * rho
      }
      bound <- solve.boundary3(power=a.cuml[i], corr=Sigma,
                               u_k=bounds)
    } else {
      bound <- u_k[i]
    }
    bounds <- c(bounds, bound)
  }
  return(bounds)
}
