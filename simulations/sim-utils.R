library(magrittr)

condense.output <- function(res){

  reject <- sapply(res, function(x) x$reject)
  est <- sapply(res, function(x) x$est)
  tstat <- sapply(res, function(x) x$tstat)
  ci <- sapply(res, function(x) x$ci) %>% t
  pval <- sapply(res, function(x) x$pval)

  bounds <- lapply(res, function(x) x$bounds)
  maxK <- sapply(bounds, length) %>% max
  v <- vector(mode="numeric", length=maxK)
  bounds <- sapply(bounds, function(x) c(x, v)[1:maxK]) %>% t

  return(list(
    reject=reject,
    est=est,
    tstat=tstat,
    ci=ci,
    pval=pval,
    bounds=bounds
  ))
}

get.param.closure <- function(i, grid){
  get.param <- function(param){
    return(grid[i, param][1])
  }
  return(get.param)
}
