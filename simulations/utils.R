nk.from.rates <- function(N, rates){
  n_cuml <- floor(N * rates)
  n_k <- c(n_cuml[1], diff(n_cuml))
  return(n_k)
}
