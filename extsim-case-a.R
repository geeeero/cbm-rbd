# --------------------------------------------------------------------------- #
# CBM systems paper - extended simulation study for revision -
# Case A: failure times as expected
# --------------------------------------------------------------------------- #

# packages, code files, example system
source("extsim-setup.R")

# simulate 100 25-cycle machines with priors as given
set.seed(2211)
br1sim5cycle20data1 <- list()
for (i in 1:20){
  br1sim5cycle20data1[[i]] <- brWeibullData(5, br1beta, br1mttf)
}


#