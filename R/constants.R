# Common alpha-spending functions.
# Use these to generate alpha spending functions
# for a fixed alpha level.

OBF.SPEND <- function(a){
  func <- function(t) 4 * (1 - pnorm(qnorm(1-a/4)/sqrt(t)))
  return(func)
}

POCOCK.SPEND <- function(a){
  func <- function(t) a * log(1 + (exp(1) - 1) * t)
  return(func)
}
