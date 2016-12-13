# --------------------------------------------------------------------------- #
# CBM systems paper - functions for simulation, version 2:
# sysrel re-evaluation also at exact failure times
# --------------------------------------------------------------------------- #

# function to simulate a single operational cycle, i.e. until the system fails,
# or it is repaired preventively
# see flow chart in paper
# sys      system reliability block diagram at t = 0
# ctypes   list giving the types of components in the system, same format as
#          used for setCompTypes(), but needs the same order as in survsign table (!!!)
# compfts  list with elements named like the vertices in sys, giving the (simulated)
#          failure time of this component -- assumed that no failure times == 0
# n0y0     list of K prior parameter pairs c(n0,y0)
# beta     vector of K fixed weibull shape parameters
# tnowstep time interval after which current system reliability is re-evaluated
#          when no component failures occur
# hor      how many time units from tnow into the future to calculate sysrelnow
# seqlen   at how many points to evaluate sysrelnow
# cu       cost of unplanned (corrective) repair action
# cp       cost of planned (preventive) repair action, cp < cu
# onecycle whether to use the one-cycle criterion or the renewal-based one
sim1cycle <- function(sys, ctypes, compfts, n0y0, beta, tnowstep, hor, seqlen = 101,
                      prior = FALSE, cu = 1, cp = 0.2, onecycle = TRUE){
  K <- length(ctypes)
  ftschron <- sort(unlist(compfts)) # when something fails
  # initializing: initial arguments for gnowhor
  tnow <- 0 # tnow: can be set to failure times
  m <- 0    # tnowstep counter, tnow grid = m * tnowstep
  sysnow <- sys
  signnow <- sign <- computeSystemSurvivalSignature(sysnow)
  wctnow <- 1:K       # which component types are now present in the system?
  ftsnow <- as.list(rep(list(NULL), K))
  failedcompsnow <- numeric(0)
  gotonext <- TRUE     # indicator that loop should go on
  failed <- FALSE      # indicator whether the system has failed
  res <- data.frame() # initialize results data frame
  while(gotonext){
    cat("tnow =", tnow, ": ") #"failedcompsnow =", failedcompsnow, "\n")
    if ((any(ftschron) > tnow) & (any(ftschron) < ((m + 1) * tnowstep))){ # failure before next tnow grid time
      tnow <- min(ftschron[ftschron > tnow])
      # check if system has failed with this component failure
    } else { # no failure before next tnow grid time
      
    }
  }



#