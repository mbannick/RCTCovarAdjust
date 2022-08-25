source("~/repos/RCTCovarAdjust/R/sim-data.R")
source("~/repos/RCTCovarAdjust/R/sim-analysis.R")
source("~/repos/RCTCovarAdjust/R/pvalues.R")
source("~/repos/RCTCovarAdjust/R/ci.R")
source("~/repos/RCTCovarAdjust/R/point-estimate.R")

# Get variances for the monitor function and final function based
# on potentially known values.
variance.closure <- function(rho, est_var){

  get.variance <- function(ancova=FALSE){
    if(est_var){
      return(NA)
    } else {
      if(ancova){
        return(rho^2)
      } else {
        return(1)
      }
    }
  }
}

procedure.closure <- function(monitor, final, correct, rates,
                              a.func, v.func, b.func, est.bounds){

  procedure <- function(data_list){
    cat(".")

    # Get the variance functions for monitoring and final stage
    monitor_var <- v.func(monitor == "ancova")
    final_var <- v.func(final == "ancova")

    # Get the functions to get test statistic for monitoring and final stage
    monitor.func <- fit.model.closure(monitor == "ancova", known_var=monitor_var)
    final.func <- fit.model.closure(final == "ancova", known_var=final_var)

    # Total alpha to spend
    total.alpha <- a.func(1)

    i <- 0
    bounds <- c()
    reject <- FALSE

    # Get sample sizes for this simulation at each stage
    n_K <- nrow(data_list[[length(data_list)]]$X)
    n_k <- rates * n_K
    K <- length(n_k)

    # If we are not estimating boundaries
    # (so, getting bounds at the design stage)
    # potentially inflate all the monitoring boundaries
    # if "correct". Otherwise, boundaries are not corrected.
    browser()
    if(!est.bounds){
      inflate <- (monitor != final) & correct
      pre_bounds <- b.func(inflate)
    }

    # Loop through each of the stages of the trial.
    while(!reject & (i < length(data_list))){

      # Data list[i] contains all data collected up through stage i
      i <- i + 1
      X <- data_list[[i]]$X
      y <- data_list[[i]]$y

      if(length(unique(X[, "t"])) == 1){
        warning("Cannot perform analysis if all randomized to same group.")
        return()
      }

      # Get an estimate (or true value) of
      # rho (the reduction in variance using to ANCOVA).
      if(is.na(monitor_var)){
        rho <- estimate.rho(X, y)
      } else {
        rho <- sqrt(v.func(ancova=TRUE) / v.func(ancova=FALSE))
      }
      browser()

      # End stage tells us if we're at the end of the trial
      end_stage <- i == length(data_list)

      # Match tells us if there is no switch between methods
      # from monitoring to final stage inference
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
        if(est.bounds){
          bound <- b.func(prev_bounds=bounds, corr=corr)
        } else {
          bound <- pre_bounds[i]
        }

        reject <- abs(test$tstat) >= bound

      } else {

        # Get the test statistic for monitoring
        test <- monitor.func(X, y)

        # Modify the correlation matrix without
        # any switching, since this is to derive the monitoring bound
        corr <- corr.mat(n_k[1:i])
        if(est.bounds){
          bound <- b.func(prev_bounds=bounds, corr=corr)
        } else {
          bound <- pre_bounds[i]
        }

        # Perform hypothesis test
        reject <- abs(test$tstat) >= bound
        if(reject){
          crossed_lower <- (test$tstat <= bound)
        }
      }
      bounds <- rbind(bounds, c(-bound, bound))
    }

    final_est <- final.func(X, y)
    est <- final_est$delta

    ancova_monitor <- monitor == "ancova"
    ancova_test <- final == "ancova"

    if(!correct) ancova_test <- ancova_monitor

    pval <- get.pvalue.sw(obs=final_est$tstat, u_k=bounds,
                          n_k=n_k[1:i], k_r=i, rho=rho,
                          ancova_monitor=ancova_monitor,
                          ancova_test=ancova_test,
                          last_stage=end_stage,
                          crossed_lower=crossed_lower)
    ci <- get.confint.sw(est=est, sd_K=final_est$variance**0.5,
                         n_k=n_k[1:i], k_r=i, u_k=bounds,
                         alpha=total.alpha, rho=rho,
                         ancova_monitor=ancova_monitor,
                         ancova_test=ancova_test,
                         last_stage=end_stage,
                         crossed_lower=crossed_lower)
    point <- get.point.sw(est=est, sd_K=final_est$variance**0.5,
                          n_k=n_k[1:i], k_r=i, u_k=bounds,
                          rho=rho,
                          ancova_monitor=ancova_monitor,
                          ancova_test=ancova_test,
                          last_stage=end_stage,
                          crossed_lower=crossed_lower)

    return(list(reject=reject, est=est, tstat=final_est$tstat,
                smean=final_est$smean, point=point,
                ci=ci, pval=pval, bounds=bounds,
                naive_ci=final_est$ci))
  }
  return(procedure)
}
