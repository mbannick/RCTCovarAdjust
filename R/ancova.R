fit.regression <- function(data, out.col, t.col, var.cols){

  # Get dimensions of problem
  n <- nrow(data)
  p <- length(var.cols) + 1 # plus one for treatment col

  # Get the outcome and design matrix
  # Design matrix is based on what is included in var.cols
  # If the covariate is in in var.cols, then it's ANCOVA otherwise it's ANOVA
  y <- data[, out.col]
  X <- data[, c(t.col, var.cols)]

  # Get the estimates
  esti <- solve(t(X) %*% X) %*% t(X) %*% y

  # Calculate residual variance
  residu <- y - X %*% esti

  # When we want to estimate the residual variance, this is what we'll do
  # but then we would need to use bounds for t-distribution
  shat <- sum(residu**2) / (n - p)

  # \Delta estimate
  delta <- esti[t.col-1]

  # Compute standardized test statistic
  tstat <- sqrt(n) * delta / sqrt(shat)

  # Variance-covariance of estimates
  vco <- shat * solve(t(X) %*% X)

  return(list(bhat=esti, vcov=vco, tstat=tstat, shat=shat))
}
