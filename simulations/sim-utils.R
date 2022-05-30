library(magrittr)

condense.output <- function(res){

  reject <- sapply(res, function(x) x$reject)
  est <- sapply(res, function(x) x$est)
  tstat <- sapply(res, function(x) x$tstat)
  smean <- sapply(res, function(x) x$smean)
  ci <- sapply(res, function(x) x$ci) %>% t
  pval <- sapply(res, function(x) x$pval)
  point <- sapply(res, function(x) x$point)
  naive_ci <- sapply(res, function(x) x$naive_ci) %>% t

  l.bounds <- lapply(res, function(x) x$bounds[, 1])
  u.bounds <- lapply(res, function(x) x$bounds[, 2])
  maxK <- sapply(l.bounds, length) %>% max
  v <- vector(mode="numeric", length=maxK)
  l.bounds <- sapply(l.bounds, function(x) c(x, v)[1:maxK]) %>% t
  u.bounds <- sapply(u.bounds, function(x) c(x, v)[1:maxK]) %>% t

  return(list(
    reject=reject,
    est=est,
    tstat=tstat,
    smean=smean,
    point=point,
    ci=ci,
    naive_ci=naive_ci,
    pval=pval,
    l.bounds=l.bounds,
    u.bounds=u.bounds
  ))
}

get.param.closure <- function(i, grid){
  get.param <- function(param){
    return(grid[i, param][1])
  }
  return(get.param)
}
