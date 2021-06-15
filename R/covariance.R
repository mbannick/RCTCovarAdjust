#' Generate basic information fraction covariance matrix
#'
#' @example
#' corr.mat(1:3/3)
#' corr.mat(c(0.1, 0.4, 0.5))
corr.mat <- function(t_k, rho=1, extra=FALSE){

  K <- length(t_k)

  dim <- K
  if(extra) dim <- dim + 1

  Sigma <- matrix(NA, nrow=dim, ncol=dim)
  for(i in 1:K){
    for(j in 1:K){
      Sigma[i, j] <- sqrt(min(t_k[i], t_k[j])/max(t_k[i], t_k[j]))
    }
  }
  if(extra){
    Sigma[K+1, 1:(K-1)] <- Sigma[K, 1:(K-1)] * rho
    Sigma[1:(K-1), K+1] <- Sigma[1:(K-1), K] * rho
    Sigma[K:(K+1), K:(K+1)] <- 1
    Sigma[K, K+1] <- rho
    Sigma[K+1, K] <- rho
  } else {
    Sigma[K, 1:(K-1)] <- Sigma[K, 1:(K-1)] * rho
    Sigma[1:(K-1), K] <- Sigma[1:(K-1), K] * rho
  }
  return(Sigma)
}
