in_dir <- "~/Documents/FileZilla/rct/run-18-04-22-1/"
summ <- fread(paste0(in_dir, "params.csv"))

params <- summ[n == 50 & rho == 0.1 & delta == 0.5 & est_var == TRUE & afunc == "obf"]

load(paste0(in_dir, "params_301.RData"))
result_nocorr <- result

load(paste0(in_dir, "params_365.RData"))
result_corr <- result

# THE REASON THAT THESE THINGS ARE THE SAME IS BECAUSE WE ONLY EVER CHANGE
# THE LAST BOUNDARY, SO THE NAIVE ESTIMATES WILL BE THE SAME REGARDLESS
# OF WHETHER OR NOT WE CORRECT THE BOUNDS.
