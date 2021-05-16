# CUTOFFS ON THE LIKELIHOOD RATIO SCALE
rm(list=ls())

# SETUP --------------------------------

OBF <- TRUE
K <- 5

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
get.joint.density <- function(n_k, lr_bounds, delta=0, rho=1,
                              lower=-10, upper=10, gridsize=1e5){

  # QUALITY CHECK
  if(length(n_k)-1 != length(lr_bounds)) stop("Different lengths for sample
                                            size and critical boundaries.")

  # NUMBER OF STAGES
  K <- length(n_k)

  # GET THE INFORMATION FRACTION
  tau <- get.tau(n_k)

  # CREATE GRID FOR DENSITIES
  grid <- seq(lower, upper, length.out=gridsize)
  dx <- (upper - lower) / gridsize

  # CONVERT LIKELIHOOD RATIO SCALE BOUNDS
  # TO SCORE STATISTIC BOUNDS AND DEFINE CONTINUATION
  # AND REJECTION REGIONS
  if(K > 1){
    s_bounds <- convert.bounds(lr_bounds, n_k[1:(K-1)])
    con <- sapply(s_bounds, function(x) abs(grid) <= x)
    rej <- sapply(s_bounds, function(x) abs(grid) > x)
  }

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

  # GET THE EFFECT SIZE
  get.eff <- function(i) delta * sqrt(n_k[1]) * sum(tau[1:i])

  # DO THE FIRST DENSITY
  f1 <- dnorm(grid, mean=get.eff(1))
  if(K == 1){
    densities[, 1] <- f1
  } else {
    f1 <- get.truncations(i=1, f=f1)
    densities[, 1] <- f1$R

    # FILL THE PREVIOUS DENSITY WITH THE CONTINUATION
    # OF THE FIRST DENSITY -- TO START THE CONVOLUTION
    f.prev <- f1$C

    # FUNCTION TO GET THE VARIANCE OF THE INCREMENT
    # WHICH IS USUALLY JUST TAU_K BUT IF RHO != 1, WHICH WOULD MEAN
    # THAT THERE IS A PROGNOSTIC COVARIATE, AND IF ANCOVA IS BEING
    # PERFORMED AT THE LAST STAGE, THEN IT'S INFLATED
    get.rho.inflation <- function(i) 2 * (1 - rho) * sum(tau[1:(i-1)])

    for(i in 2:K){

      # COMPUTE VARIANCE
      v <- tau[i]
      if(i ==  K) v <- v + get.rho.inflation(i)

      # SPECIFY THE DENSITY
      m <- get.eff(i) - get.eff(i-1)
      fi <- dnorm(grid, mean=m, sd=sqrt(v))

      # GET THE CONVOLUTION WITH THE PREVIOUS INCREMENTS
      # AND THE INDEPENDENT INCREMENT
      fi.cuml <- conv(f.prev, fi, grid=grid, delta=dx)

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
  }
  return(list(grid=grid, dens=densities))
}

# TEST IT OUT -------------------------------

n_k <- rep(200, 5)
s_bounds <- convert.bounds(u_k, n_k)

# THIS IS ACCURATE PROPORTIONALLY TO THE GRIDSIZE
dens <- get.joint.density(n_k=n_k, lr_bounds=u_k[1:(K-1)],
                          gridsize=1e5)
sum(dens) * 2e-4

dens.ancova <- get.joint.density(n_k=n_k, lr_bounds=u_k, rho=0.9,
                                 gridsize=1e5)
sum(dens$dens) * 2e-4
#
# # CHECK FOR TYPE 1 ERROR
#
# grid <- seq(-10, 10, length.out=1e5)
# reject <- sapply(s_bounds, function(x) abs(grid) > x)
# sum(reject * dens) * 2e-4
# sum(reject * dens.ancova) * 2e-4 # type 1 error inflated
#
# # WHAT HAPPENS IF WE CHANGE THE EFFECT SIZE
# # THIS SHOULD BE A MUCH LARGER PROBABILITY OF REJECTION!
# dens.d <- get.joint.density(n_k=n_k, lr_bounds=u_k, delta=0.2,
#                           gridsize=1e5)
# sum(reject * dens.d) * 2e-4
#
# # NOTE: IF DELTA GETS TOO LARGE, YOU NEED TO CHANGE THE UPPER AND LOWER
# # GRID FOR X, BECAUSE THEN YOU RUN OUT OF ROOM IN THE GRID!
#
# # COMPARE TO SIMUALTED VERSION
# u_k1 <- c(4.468297, 3.381480, 2.710148, 2.305284, 2.031209)
# u_k2 <- c(4.468297, 3.381480, 2.710148, 2.305284, 2.033433)
#
# dens.1 <- get.joint.density(n_k=n_k, lr_bounds=u_k1, delta=0, rho=1,
#                             gridsize=1e5)
# dens.2 <- get.joint.density(n_k=n_k, lr_bounds=u_k1, delta=0, rho=0.9950372,
#                             gridsize=1e5)
#
# s_bounds1 <- convert.bounds(u_k1, n_k)
# s_bounds2 <- convert.bounds(u_k2, n_k)
# grid <- seq(-10, 10, length.out=1e5)
# reject1 <- sapply(s_bounds1, function(x) abs(grid) > x)
# reject2 <- sapply(s_bounds2, function(x) abs(grid) > x)
#
#
# sum(reject1 * dens.1) * 2e-4
# sum(reject1 * dens.2) * 2e-4
# sum(reject2 * dens.2) * 2e-4
