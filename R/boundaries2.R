source("R/recursive.R")


solve.boundary2 <- function(power, dens, grid, lower=0, upper=50, ...){

  dx <- (max(grid) - min(grid)) / length(grid)
  get.power <- function(u) sum(dens[which(abs(grid) >= u)]) * dx
  u.search <- function(u) get.power(u) - power
  u <- uniroot(u.search, lower=0, upper=50)

  return(u$root)
}

get.boundaries.aspend2 <- function(a.func, a, rates, N,
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
  if(!0 %in% rates) rates <- c(0, rates)
  a.diff <- a.spend(rates) %>% diff

  # # Create a "previous" vector that we will use to
  # # indicate having *not* rejected at previous stages
  # bound <- solve.boundary2(power=a.diff[i],
  #                          dens=dens[, K_prev+1])

  bounds <- c()
  # browser()
  for(i in 1:K){
    if(i > K_prev){
      # Create the joint densities
      dens <- get.joint.density(n_k=n_k[1:i], lr_bounds=bounds,
                                delta=0, rho=rho, ...)
      bound <- solve.boundary2(power=a.diff[i],
                              dens=dens$dens[, i], grid=dens$grid)
      browser()
      bound <- bound / sqrt(tau_s[i])
      browser()
    } else {
      bound <- u_k[i]
    }
    bounds <- c(bounds, bound)
  }

  return(bounds)

}
