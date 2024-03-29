% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boundaries.R
\name{get.boundaries.design}
\alias{get.boundaries.design}
\title{Get Pocock or OBF-style boundaries at the design stage}
\usage{
get.boundaries.design(
  rates,
  obf,
  unequal_type = FALSE,
  rho = 1,
  change = 0,
  algorithm = Miwa(steps = 1000)
)
}
\arguments{
\item{rates}{A vector of information rates (between 0 and 1)}

\item{obf}{Whether to use OBF (TRUE) or Pocock (FALSE)}

\item{unequal_type}{Correction for unequal sample sizes across stage}

\item{rho}{Fraction of variance explained by fitting ANCOVA.}

\item{change}{A vector indicating which stages use ANOVA v. ANCOVA.}
}
\description{
Get Pocock or OBF-style boundaries at the design stage
}
\examples{
# Information fractions
t <- 1:3/3

# OBF-type boundaries
get.boundaries.design(rates=t, obf=TRUE)
get.boundaries.design(rates=c(0.3, 0.9, 1.0), obf=TRUE)
get.boundaries.design(rates=c(0.3, 0.9, 1.0), obf=TRUE, unequal_type=TRUE)

# ANCOVA at last stage, R^2 = 0.5
get.boundaries.design(rates=c(0.3, 0.9, 1.0), obf=TRUE, unequal_type=TRUE, rho=sqrt(0.5), change=c(0, 0, 1))

# ANCOVA at last two stages, R^2 = 0.5
get.boundaries.design(rates=t, obf=TRUE, rho=sqrt(0.5), change=c(1, 0, 0))
}
