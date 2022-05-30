# DEBUG DIFFERENCE BETWEEN SIMULATION FIGURE VERSIONS

old <- "/Users/marlena/Documents/FileZilla/rct/run-17-02-22-6"
old <- "/Users/marlena/Documents/FileZilla/rct/run-10-04-22-2"
new <- "/Users/marlena/Documents/FileZilla/rct/run-16-05-22-5"

pold <- fread(paste0(old, "/params.csv"))
pnew <- fread(paste0(new, "/params.csv"))

cold <- pold[n == 250 & delta == 0.1 & rho == 0.1 & afunc == "obf" & stages == 3]
cnew <- pnew[n == 250 & delta == 0.1 & rho == 0.1 & afunc == "obf" & stages == 3]

setnames(cold, "V1", "old_param")
setnames(cnew, "V1", "new_param")

name <- names(cold)[!names(cold) %in% "old_param"]

df <- merge(cold, cnew, by=name)

load(paste0(old, "/params_7.RData"))
rold <- result

load(paste0(new, "/params_7.RData"))
rnew <- result



# load("/Users/marlena/Documents/FileZilla/rct/run-17-02-22-6/params_1.RData")
