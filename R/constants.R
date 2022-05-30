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
      return(list(
        tk=c(0.30, 1.00),
        uk_obf=c(3.929, 1.960),
        uk_poc=c(2.312, 2.124)
      ))
    } else if(type == 2){
      return(list(
        tk=c(0.50, 1.00),
        uk_obf=c(2.963, 1.969),
        uk_poc=c(2.157, 2.201)
      ))
    } else if(type == 3){
      return(list(
        tk=c(0.90, 1.00),
        uk_obf=c(2.094, 2.053),
        uk_poc=c(1.989, 2.241)
      ))
    } else {
      stop(msg)
    }
  } else if(stages == 3){
    if(type == 1){
      return(list(
        tk=c(0.30, 0.90, 1.00),
        uk_obf=c(3.929, 2.094, 2.053),
        uk_poc=c(2.312, 2.162, 2.342)
      ))
    } else if(type == 2){
      return(list(
        tk=c(0.33, 0.67, 1.00),
        uk_obf=c(3.710, 2.511, 1.993),
        uk_poc=c(2.279, 2.295, 2.296)
      ))
    } else if(type == 3){
      return(list(
        tk=c(0.80, 0.90, 1.00),
        uk_obf=c(2.250, 2.177, 2.072),
        uk_poc=c(2.021, 2.271, 2.332)
      ))
    } else {
      stop(msg)
    }
  } else if(stages == 4){
    if(type == 1){
      return(list(
        tk=c(0.20, 0.40, 0.90, 1.00),
        uk_obf=c(4.877, 3.357, 2.097, 2.054),
        uk_poc=c(2.438, 2.427, 2.224, 2.376)
      ))
    } else if(type == 2){
      return(list(
        tk=c(0.25, 0.50, 0.75, 1.00),
        uk_obf=c(4.333, 2.963, 2.359, 2.014),
        uk_poc=c(2.368, 2.368, 2.358, 2.350)
      ))
    } else if(type == 3){
      return(list(
        tk=c(0.30, 0.60, 0.90, 1.00),
        uk_obf=c(3.929, 2.670, 2.121, 2.063),
        uk_poc=c(2.312, 2.321, 2.318, 2.412)
      ))
    } else {
      stop(msg)
    }
  } else {
    stop(msg)
  }
}
