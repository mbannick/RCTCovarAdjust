library(data.table)
library(magrittr)
library(xtable)

# SET DIRECTORIES AND VERSION
args <- commandArgs(trailingOnly=TRUE)

VERSION <- args[1]
AFUNC <- args[2]

IN_DIR <- "/Users/marlena/Documents/FileZilla/rct/"

# READ IN VERSION
df <- fread(paste0(IN_DIR, VERSION, "/summary.csv"))

# Rename the column for median bias
setnames(df, "bias_naive_med", "med_naive")
setnames(df, "bias_corr_med", "med_corr")

# Get the results for the particular alpha spending function
df.sub <- df[afunc == AFUNC]

# Scale the bias by 100
df.sub[, bias_naive := bias_naive * 100]
df.sub[, bias_corr := bias_corr * 100]
df.sub[, med_naive := med_naive * 100]
df.sub[, med_corr := med_corr * 100]

# If we also have results when not estimating the variance, get rid of those
# Also subset only to the switching from anova to ancova
df.sub <- df.sub[monitor == "anova" & final == "ancova"]
df.sub <- df.sub[est_var == TRUE]
df.sub[, est_var := NULL]

setorder(df.sub, delta, rho, n, type)

# Format the bias column to be median (mean)
df.sub[, dev_naive := paste0(sprintf("%0.2f", med_naive),
                             " (", sprintf("%0.2f", bias_naive), ")")]
df.sub[, dev_corr := paste0(sprintf("%0.2f", med_corr),
                            " (", sprintf("%0.2f", bias_corr), ")")]

# Subset and reformat data
df.sub2 <- df.sub[, .(delta, rho, n, type,
                      dev_naive, dev_corr,
                      cover_naive, cover_corr)]
df.sub2 <- data.table(df.sub2)
df.sub2 <- dcast(delta + rho + n ~ type,
                 value.var = c("dev_naive", "dev_corr",
                               "cover_naive", "cover_corr"), data=df.sub2)
df.sub2 <- df.sub2[, .(delta, rho, n,
                       `dev_naive_inconsistent, naive`,
                       `dev_corr_inconsistent, naive`,
                       `dev_corr_inconsistent, corrected`,
                       `cover_naive_inconsistent, naive`,
                       `cover_corr_inconsistent, naive`,
                       `cover_corr_inconsistent, corrected`)]

setnames(df.sub2, c("delta", "rho", "n",
                    "bias (0)", "bias (1)", "bias (2)",
                    "cover (0)", "cover (1)", "cover (2)"))

# Format table and print
addtorow <- list()
addtorow$pos <- seq(4, nrow(df.sub2), by=4) %>% as.list
addtorow$command <- rep("\\hline \n", length(addtorow$pos))

tab <- xtable(df.sub2, align=rep("c", 10), digits=2,
              caption=paste0("Simulation results for ", AFUNC,
                             " version: ", VERSION))
print(tab, include.rownames=FALSE,
      add.to.row = addtorow)
