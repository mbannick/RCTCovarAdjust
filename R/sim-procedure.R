source("~/repos/RCTCovarAdjust/R/sim-data.R")
source("~/repos/RCTCovarAdjust/R/sim-analysis.R")
source("~/repos/RCTCovarAdjust/R/pvalues.R")
source("~/repos/RCTCovarAdjust/R/ci.R")

# Get variances for the monitor function and final function based
# on potentially known values.
variance.closure <- function(cov_std, obs_std, beta, est_var){

  get.variance <- function(ancova=FALSE){
    if(est_var){
      return(NA)
    } else {
      if(ancova){
        return(obs_std**2)
      } else {
        return(obs_std**2 + beta**2 * cov_std**2)
      }
    }
  }
}

procedure.closure <- function(monitor, final, correct, rates,
                              a.func, v.func, b.func){

  procedure <- function(data_list){

    monitor_var <- v.func(monitor == "ancova")
    final_var <- v.func(monitor == "ancova")

    monitor.func <- fit.model.closure(monitor == "ancova", known_var=monitor_var)
    final.func <- fit.model.closure(final == "ancova", known_var=final_var)

    total.alpha <- a.func(1)

    i <- 0
    bounds <- c()
    reject <- FALSE

    n_K <- nrow(data_list[[length(data_list)]]$X)
    n_k <- rates * n_K

    while(!reject & (i < length(data_list))){

      i <- i + 1
      X <- data_list[[i]]$X
      y <- data_list[[i]]$y

      rho <- estimate.rho(X, y)

      end_stage <- i == length(data_list)
      match <- monitor == final

      if(end_stage){

        # Get the final test statistic
        test <- final.func(X, y)

        # Modify correlation matrix
        # to account for switch
        if((match | !correct) & rho < 1.){
          corr <- corr.mat(n_k[1:i], 1)
        } else {
          corr <- corr.mat(n_k[1:i], rho=rho, mis=c(rep(F, i-1), T))
        }

        # Calculate the final boundary and perform test
        bound <- b.func(prev_bounds=bounds, corr=corr)
        # bound <- static.bounds[i]
        reject <- abs(test$tstat) >= bound

      } else {

        # Get the test statistic for monitoring
        test <- monitor.func(X, y)

        # Modify the correlation matrix without
        # any switching, since this is to derive the monitoring bound
        corr <- corr.mat(n_k[1:i])
        bound <- b.func(prev_bounds=bounds, corr=corr)
        # bound <- static.bounds[i]

        # Perform hypothesis test
        reject <- abs(test$tstat) >= bound
      }
      bounds <- rbind(bounds, c(-bound, bound))
    }
    if(i == 1){
      u_k <- c()
    } else {
      if(!match & reject & !end_stage & correct){
        u_k <- bounds
      } else {
        u_k <- matrix(bounds[1:(i-1),], ncol=2)
      }
    }
    final_est <- final.func(X, y)
    est <- final_est$delta

    ancova_monitor <- monitor == "ancova"
    ancova_test <- final == "ancova"

    if(!correct) ancova_test <- ancova_monitor

    pval <- get.pvalue.sw(obs=final_est$tstat, u_k=u_k,
                          n_k=n_k[1:i], rho=rho,
                          ancova_monitor=ancova_monitor,
                          ancova_test=ancova_test,
                          last_stage=end_stage)
    # ci <- get.confint.sw(est=est, sd_K=final_est$variance**0.5,
    #                      n_k=n_k[1:i], u_k=u_k,
    #                      alpha=total.alpha, rho=rho,
    #                      ancova_monitor=(monitor == "ancova"),
    #                      ancova_test=(final == "ancova" & correct),
    #                      last_stage=end_stage)
    ci <- c(0, 0)

    return(list(reject=reject, est=est, tstat=final_est$tstat,
                smean=final_est$smean,
                ci=ci, pval=pval, bounds=bounds))
  }
  return(procedure)
}
