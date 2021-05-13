# CUTOFFS ON THE LIKELIHOOD RATIO SCALE
rm(list=ls())

# SETUP --------------------------------

OBF <- FALSE
K <- 4

# OBrien and Fleming Critical Value Constants
bound_obf <- c(1.9600, 2.7965, 3.4711, 4.0486, 4.5617,
               5.0283, 5.5490, 5.8611, 6.2395, 6.5981,
               6.9396, 7.2663, 7.5799, 7.8820, 8.1736)

# Pocock Critical Value Constants
bound_poc <- c(1.9600, 2.1783, 2.2895, 2.3613, 2.4132,
               2.4532, 2.4855, 2.5123, 2.5352, 2.5550,
               2.5724, 2.5880, 2.6019, 2.6146, 2.6261)

# COMPUTE CRITICAL VALUES
if(OBF){
  u_k <- bound_obf[K] / sqrt(1:K)
} else {
  u_k <- bound_poc[K] / rep(1, K)
}

# HELPER FUNCTIONS --------------------------

# CONVOLUTION FUNCTION
conv <- function(f1, f2, grid, delta){
  flip <- convolve(f1, f2) * delta
  f.new <- flip[
    c(ceiling(length(grid)/2):length(grid),
    1:(ceiling(length(grid)/2)-1))
  ]
  return(f.new)
}

# GET THE INFORMATION FRACTION FUNCTION TAU
get.tau <- function(n_k){
  tau <- n_k / n_k[1]
  return(tau)
}

# CONVERT BOUNDS TO THE SCORE SCALE
convert.bounds <- function(bounds, n_k){
  tau <- get.tau(n_k)
  tau_s <- cumsum(tau)
  bounds <- bounds * sqrt(tau_s)
  return(bounds)
}

# MAIN FUNC -- GET THE JOINT DENSITY OF THE STATISTICS
# AND THE REJECTION STAGE
get.joint.density <- function(n_k, lr_bounds,
                              lower=-10, upper=10, gridsize=1e5){

  # QUALITY CHECK
  if(length(n_k) != length(lr_bounds)) stop("Different lengths for sample
                                            size and critical boundaries.")

  # NUMBER OF STAGES
  K <- length(n_k)

  # GET THE INFORMATION FRACTION
  tau <- get.tau(n_k)

  # CONVERT LIKELIHOOD RATIO SCALE BOUNDS
  # TO SCORE STATISTIC BOUNDS
  s_bounds <- convert.bounds(lr_bounds, n_k)

  # CREATE GRID FOR DENSITIES
  grid <- seq(lower, upper, length.out=gridsize)
  delta <- (upper - lower) / gridsize

  # DEFINE CONTINUATION AND REJECTION REGIONS
  con <- sapply(s_bounds, function(x) abs(grid) <= x)
  rej <- sapply(s_bounds, function(x) abs(grid) > x)

  # PLACEHOLDER FOR THE MATRIX OF DENSITIES FOR EACH
  # STAGE
  densities <- matrix(NA, nrow=gridsize, ncol=K)

  # FUNCTION TO TRUNCATE A DENSITY BASED
  # ON THE CONTINUATION OR REJECTION REGIONS
  get.truncations <- function(i, f){
    f.R <- f
    f.C <- f
    f.R[which(con[, i])] <- 0
    f.C[which(rej[, i])] <- 0
    return(list(R=f.R, C=f.C))
  }

  # DO THE FIRST DENSITY
  f1 <- dnorm(grid)
  f1 <- get.truncations(i=1, f=f1)
  densities[, 1] <- f1$R

  # FILL THE PREVIOUS DENSITY WITH THE CONTINUATION
  # OF THE FIRST DENSITY -- TO START THE CONVOLUTION
  f.prev <- f1$C

  for(i in 2:K){
    fi <- dnorm(grid, sd=sqrt(tau[i]))

    # GET THE CONVOLUTION WITH THE PREVIOUS INCREMENTS
    # AND THE INDEPENDENT INCREMENT
    fi.cuml <- conv(f.prev, fi, grid=grid, delta=delta)

    if(i == K){

      # IF IN THE LAST STAGE, FILL IN DENSITY WITH THE FULL CONVOLUTION
      # DOESN'T MATTER REJECTION / CONTINUATION B/C AT LAST STAGE
      densities[, i] <- fi.cuml

    } else {

      # TRUNCATE THE DENSITY BASED ON REJECTION REGION
      fi.trun <- get.truncations(i=i, f=fi.cuml)

      # FILL IN THE DENSITY FOR THIS STAGE REJECTION ONLY
      densities[, i] <- fi.trun$R

      # REPLACE THE PREVIOUS DENSITY WITH THIS ONE
      # IN THE CONTINUATION REGION
      f.prev <- fi.trun$C

    }
  }
  return(densities)
}

# TEST IT OUT -------------------------------

# THIS IS ACCURATE PROPORTIONALLY TO THE GRIDSIZE
dens <- get.joint.density(n_k=c(20, 30, 40, 20), lr_bounds=u_k,
                          gridsize=1e5)
sum(dens) * 2e-4
