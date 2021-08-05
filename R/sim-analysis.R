#' Generate a function to fit either ANOVA
#' or ANCOVA based on outcome y and design matrix X.
#'
#' @examples
#' X <- rnorm(100) %>% matrix(ncol=4)
#' beta <- rnorm(4)
#' t <- rbinom(n=25, size=1, prob=0.5)
#' y <- 1 + 2 * t + X %*% beta + rnorm(25)
#' X <- cbind(t, int=1, X)
#'
#' anova <- fit.model.closure(ancova=F)
#' ancova <- fit.model.closure(ancova=T)
#'
#' anova(X, y)
#' ancova(X, y)
fit.model.closure <- function(ancova=FALSE, known_var=NA){
  fit.model <- function(X, y){
    if(ancova){
      X.sub <- X
    } else {
      X.sub <- X[, c("t", "int")]
    }
    n <- nrow(X.sub)
    p <- ncol(X.sub)

    est <- solve(t(X.sub) %*% X.sub) %*% t(X.sub) %*% y
    res <- y - X.sub %*% est

    if(is.na(known_var)){
      vhat <- sum(res**2) / (n - p)
    } else {
      vhat <- known_var
    }

    delta <- est["t",]
    tstat <- sqrt(n) * delta / sqrt(vhat)
    smean <- tstat / sqrt(n)

    return(list(delta=delta, tstat=tstat, variance=vhat,
                smean=smean))
  }
}

#' Estimate the variance reduction from ANOVA
#' to ANCOVA
#'
#' @examples
#' X <- rnorm(100) %>% matrix(ncol=4)
#' beta <- rnorm(4)
#' t <- rbinom(n=25, size=1, prob=0.5)
#' y <- 1 + 2 * t + X %*% beta + rnorm(25)
#' X <- cbind(t, int=1, X)
#'
#' estimate.rho(X, y)
estimate.rho <- function(X, y){
  anova <- fit.model.closure(ancova=FALSE)
  ancova <- fit.model.closure(ancova=TRUE)

  s.anova <- anova(X, y)$variance ** 0.5
  s.ancova <- ancova(X, y)$variance ** 0.5

  return(min(s.ancova / s.anova, 1.))
}
