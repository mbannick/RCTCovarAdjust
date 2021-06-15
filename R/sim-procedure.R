source("R/sim-data.R")
source("R/sim-analysis.R")

procedure.closure <- function(monitor, final, rates, a.func){

  bound.func <- get.boundary.closure(a.func=a.func, rates=rates)

  procedure <- function(data_list){

    monitor.func <- fit.model.closure(monitor == "ancova")
    final.func <- fit.model.closure(final == "ancova")

    total.alpha <- a.func(1)

    i <- 0
    bounds <- c()
    reject <- FALSE

    while(!reject & (i < length(data_list))){

      i <- i + 1
      X <- data_list[[i]]$X
      y <- data_list[[i]]$y

      corr <- corr.mat(rates[1:i])
      rho <- estimate.rho(X, y)

      end_stage <- i == length(data_list)
      match <- monitor == final

      if(end_stage){

        # Get the final test statistic
        test <- final.func(X, y)

        # Modify correlation matrix
        # to account for switch
        if(match){
          corr <- corr.mat(rates[1:i], 1)
        } else {
          corr <- corr.mat(rates[1:i], rho)
        }

        # Calculate the final boundary and perform test
        bound <- bound.func(prev_bounds=bounds, corr=corr)
        reject <- abs(test$tstat) >= bound

      } else {

        # Get the test statistic for monitoring
        test <- monitor.func(X, y)

        # Modify the correlation matrix without
        # any switching, since this is to derive the monitoring bound
        corr <- corr.mat(rates[1:i])
        bound <- bound.func(prev_bounds=bounds, corr=corr)

        # Perform hypothesis test
        reject <- abs(test$tstat) >= bound

        # If we reject in the interim stage, and there
        # is a switch, need to add another test statistic
        # to the stage-wise ordering
        if(!match & reject){
          corr <- corr.mat(rates[1:i], rho, extra=TRUE)
        }
      }
      bounds <- c(bounds, bound)
    }
    if(!match & reject & !end_stage){
      u_k <- bounds
    } else {
      u_k <- bounds[1:(i-1)]
    }

    final <- final.func(X, y)
    est <- final$delta
    ci <- get.confint.sw(est=est, sd_K=final$variance,
                         n_K=nrow(X), u_k=u_k,
                         alpha=total.alpha, corr=corr)
    pval <- get.pvalue.sw(obs=final$tstat, u_k=u_k,
                          corr=corr)

    return(list(reject=reject, est=est, tstat=final$tstat,
                ci=ci, pval=pval, bounds=bounds))
  }
  return(procedure)
}
