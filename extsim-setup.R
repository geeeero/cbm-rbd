# --------------------------------------------------------------------------- #
# CBM systems paper - extended simulation study for revision - simulation setup
# --------------------------------------------------------------------------- #

# install packages if necessary
install.packages(c("actuar", "igraph", "ReliabilityTheory", "reshape2"))
# on LISA, need to install "actuar", "igraph", "ReliabilityTheory"

# load packages
library(actuar)
library(ReliabilityTheory)
library(reshape2) # might not be necessary

# code files
source("weibull-sys-functions.R") # contains also plotting functions that work only with ggplot2 and gridExtra
source("cbm-sim2.R")

# the example system
source("cbm-sys-other.R")

#