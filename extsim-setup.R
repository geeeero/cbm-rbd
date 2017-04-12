# --------------------------------------------------------------------------- #
# CBM systems paper - extended simulation study for revision - simulation setup
# --------------------------------------------------------------------------- #

# install packages if necessary
install.packages(c("actuar", "ReliabilityTheory", "reshape2"))

# load packages
library(actuar)
library(ReliabilityTheory)
library(reshape2)

# code files
source("weibull-sys-functions.R") # contains also plotting functions that work only with ggplot2 and gridExtra
source("cbm-sim2.R")

# the example system
source("cbm-sys-other.R")

#