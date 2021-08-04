# Common alpha-spending functions.
# Use these to generate alpha spending functions
# for a fixed alpha level.

obf.spend <- function(a){
  func <- function(t) 4 * (1 - pnorm(qnorm(1-a/4)/sqrt(t)))
  return(func)
}

pocock.spend <- function(a){
  func <- function(t) a * log(1 + (exp(1) - 1) * t)
  return(func)
}

spend <- function(a, type){
  if(type == "obf") return(obf.spend(a))
  if(type == "pocock") return(pocock.spend(a))
}

# Common Information Fractions

info.fractions <- function(stages, type){

  msg <- "Unrecognized argument."

  if(stages == 2){
    if(type == 1){
      return(c(0.30, 1.00))
    } else if(type == 2){
      return(c(0.50, 1.00))
    } else if(type == 3){
      return(c(0.90, 1.00))
    } else {
      stop(msg)
    }
  } else if(stages == 3){
    if(type == 1){
      return(c(0.30, 0.90, 1.00))
    } else if(type == 2){
      return(c(0.33, 0.67, 1.00))
    } else if(type == 3){
      return(c(0.80, 0.90, 1.00))
    } else {
      stop(msg)
    }
  } else if(stages == 4){
    if(type == 1){
      return(c(0.20, 0.40, 0.90, 1.00))
    } else if(type == 2){
      return(c(0.25, 0.50, 0.75, 1.00))
    } else if(type == 3){
      return(c(0.30, 0.60, 0.90, 1.00))
    } else {
      stop(msg)
    }
  } else {
    stop(msg)
  }
}
