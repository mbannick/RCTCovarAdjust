source("~/repos/RCTCovarAdjust/R/sim-data.R")
source("~/repos/RCTCovarAdjust/R/sim-analysis.R")
source("~/repos/RCTCovarAdjust/R/pvalues.R")
source("~/repos/RCTCovarAdjust/R/ci.R")

# Get variances for the monitor function and final function based
# on potentially known values.
get.variances <- function(monitor, final, sd_anova=NA, sd_ancova=NA){
  monitor_var <- NA
  final_var <- NA

  if(!is.na(sd_anova) | !is.na(sd_ancova)){
    if(!is.na(sd_anova) & !is.na(sd_ancova)){
      stop("Need to pass both as non-NA values, or none..")
    }
  }

  if(!is.na(sd_anova) & !is.na(sd_ancova)){
    if(monitor == "anova"){
      monitor_var <- sd_anova**2
    } else {
      monitor_var <- sd_ancova**2
    }
    if(final == "anova"){
      final_var <- sd_anova**2
    } else {
      final_var <- sd_ancova**2
    }
  }
  return(list(monitor=monitor_var, final=final_var))
}

procedure.closure <- function(monitor, final, correct, rates, a.func,
                              sd_anova=NA, sd_ancova=NA){

  bound.func <- get.boundary.closure(a.func=a.func, rates=rates)

  procedure <- function(data_list){

    vars <- get.variances(monitor=monitor, final=final,
                          sd_anova=sd_anova, sd_ancova=sd_ancova)

    monitor.func <- fit.model.closure(monitor == "ancova", known_var=vars$monitor)
    final.func <- fit.model.closure(final == "ancova", known_var=vars$final)

    total.alpha <- a.func(1)

    i <- 0

    bounds <- c()
    # TEMPORARY HACK
    static.bounds <- c(2.963, 1.969)
    reject <- FALSE

    n_K <- nrow(data_list[[length(data_list)]]$X)
    n_k <- rates * n_K

    while(!reject & (i < length(data_list))){

      i <- i + 1
      X <- data_list[[i]]$X
      y <- data_list[[i]]$y

      corr <- corr.mat(n_k[1:i])
      rho <- estimate.rho(X, y)

      end_stage <- i == length(data_list)
      match <- monitor == final

      if(end_stage){

        # Get the final test statistic
        test <- final.func(X, y)

        # Modify correlation matrix
        # to account for switch
        if(match | !correct){
          corr <- corr.mat(n_k[1:i], 1)
        } else {
          corr <- corr.mat(n_k[1:i], rho=rho, mis=c(rep(F, i-1), T))
        }

        # Calculate the final boundary and perform test
        # bound <- bound.func(prev_bounds=bounds, corr=corr)
        bound <- static.bounds[i]
        reject <- abs(test$tstat) >= bound

      } else {

        # Get the test statistic for monitoring
        test <- monitor.func(X, y)

        # Modify the correlation matrix without
        # any switching, since this is to derive the monitoring bound
        corr <- corr.mat(n_k[1:i])
        # bound <- bound.func(prev_bounds=bounds, corr=corr)
        bound <- static.bounds[i]

        # Perform hypothesis test
        reject <- abs(test$tstat) >= bound

        # If we reject in the interim stage, and there
        # is a switch, need to add another test statistic
        # to the stage-wise ordering
        if(!match & reject & correct){
          corr <- corr.mat(c(n_k[1:i], n_k[i]), rho=rho,
                           mis=c(rep(F, i), T))
        }
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
    ci <- get.confint.sw(est=est, sd_K=final_est$variance,
                         n_k=n_k[1:i], u_k=u_k,
                         alpha=total.alpha, rho=rho,
                         ancova_monitor=(monitor == "ancova"),
                         ancova_test=(final == "ancova"),
                         last_stage=end_stage)
    pval <- get.pvalue.sw(obs=final_est$tstat, u_k=u_k,
                          n_k=n_k[1:i], rho=rho,
                          ancova_monitor=(monitor == "ancova"),
                          ancova_test=(final == "ancova"),
                          last_stage=end_stage)

    return(list(reject=reject, est=est, tstat=final_est$tstat,
                smean=final_est$smean,
                ci=ci, pval=pval, bounds=bounds))
  }
  return(procedure)
}
