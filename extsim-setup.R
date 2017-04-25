# --------------------------------------------------------------------------- #
# CBM systems paper - extended simulation study for revision - simulation setup
# --------------------------------------------------------------------------- #

# install packages if necessary
install.packages(c("actuar", "igraph", "ReliabilityTheory", "reshape2"))
# on LISA, if they install "actuar", "igraph", and "mcmc",
# then I need to install only "ReliabilityTheory",
# which depends on the below other three packages:
install.packages("FRACTION")
install.packages("PhaseType")
install.packages("HI")
install.packages("ReliabilityTheory")



# load packages
library(actuar)
library(ReliabilityTheory)
library(reshape2) # might not be necessary
#library(plyr)

# code files
source("weibull-sys-functions.R") # contains also plotting functions that work only with ggplot2 and gridExtra
source("cbm-sim2.R")

# the example system
source("cbm-sys-other.R")

#