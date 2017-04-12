# CBM systems paper - extended simulation study for revision -
# Case D: fluctuating failure time behaviour (due to external factors)
# --------------------------------------------------------------------------- #

# packages, code files, example system
source("extsim-setup.R")

# simulate 100 25-cycle machines
set.seed(2411)
repetitions <- 100
opcycles <- 25
Dsimdata <- list()
mttffactor <- rep(sample(c(0.5, 1, 2), size = repetitions/2, replace = TRUE), each = 2)
for (i in 1:repetitions){
  Dsimdata[[i]] <- brWeibullData(opcycles, br1beta, br1mttf*mttffactor[i])
}

#