parse.args <- function(arg, type_fun){
  vec <- strsplit(arg, ",")[[1]]
  vec <- type_fun(vec)
  return(vec)
}
