#' Generate basic information fraction covariance matrix
#'
#' @example
#' corr.mat(1:3/3 * 30)
#' corr.mat(c(0.1, 0.4, 0.5)*50)
corr.mat <- function(n_k, rho=1, mis=F, sme=F){
  K <- length(n_k)

  # Use the mismatch vector and rho to determine
  # Correlation when switching between ANOVA and ANCOVA
  if(length(mis) == 1){
    mis <- rep(mis, K)
  } else {
    if(length(mis) != K) stop()
  }
  rho.mat <- matrix(1, nrow=K, ncol=K)

  for(i in 1:K){
    for(j in 1:K){
      if(mis[i] != mis[j]){
        rho.mat[i, j] <- rho
        rho.mat[j, i] <- rho
      }
    }
  }
  diag(rho.mat) <- 1

  # Use sample mean indicator (rather than likelihood ratio
  # statistic to get the correlation *and* variance when switching
  # to sample mean.)
  if(length(sme) == 1) sme <- rep(sme, K)
  sme.sd <- (1 / sqrt(n_k)) * (sme) + 1 * (!sme)
  sme.mat <- sme.sd %*% t(sme.sd)

  # Use sample sizes to get the correlation based on the
  # information fraction.
  Sigma <- matrix(NA, nrow=K, ncol=K)
  for(i in 1:K){
    for(j in 1:K){
      Sigma[i, j] <- sqrt(min(n_k[i], n_k[j])/max(n_k[i], n_k[j]))
    }
  }

  # Element-wise multiplication of the three components.
  mat <- Sigma * sme.mat * rho.mat
  return(mat)
}
