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
# cc       cost of corrective repair action
# cp       cost of preventive repair action, cp < cu
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
  # calculate tstarnow at time tnow = 0
  gnowvec <- gnowhor(signnow, n0y0[wctnow], beta[wctnow], ftsnow[wctnow], tnow, hor = hor,
                     seqlen = seqlen, prior = prior, cu = cu, cp = cp, onecycle = onecycle)
  taustarnow <- gnowvec$tau[which.min(gnowvec$gnow)] - tnow
  # initialize results data frame with entries for tnow = 0
  res <- data.frame(tnow = tnow, failed = failed, taustar = taustarnow)
  # stop already now if taustarnow smaller than tnowstep
  if (taustarnow < tnowstep){
    cost <- cp
    gotonext <- FALSE
  }
  while(gotonext){
    cat("tnow =", tnow, ": ") #"failedcompsnow =", failedcompsnow, "\n")
    nextgridtime <- (m + 1) * tnowstep
    if ((any(ftschron) > tnow) & (any(ftschron) < nextgridtime){ # failure before next tnow grid time
      tnow <- min(ftschron[ftschron > tnow])
      # update all gnowvec inputs to current tnow
      failedcompsnow <- ftschron[ftschron <= tnow]
      sysnow <- induced_subgraph(sys, vids=V(sys)[!(name %in% names(failedcompsnow))])
      # catch error when sysnow now contains vertices s and t only (suddenly all components fail)
      if(distances(sysnow, "s", "t") < Inf)
        signnow <- computeSystemSurvivalSignature(sysnow)
      else
        signnow <- data.frame(Probability = 0)
      # to which fts list element do the failure times in failedcompsnow belong?
      ftsnowindex <- sapply(ctypes, function(ctypesl) names(failedcompsnow) %in% ctypesl)
      if (length(failedcompsnow) == 1)
        ftsnowindex <- which(ftsnowindex)
      else
        ftsnowindex <- apply(ftsnowindex, 1, which)
      # write failure times in corresponding list element
      for (k in 1:K)
        ftsnow[[k]] <- failedcompsnow[ftsnowindex == k]
      # which component types are now present? (subset of 1:K)
      wctnow <- which(names(ctypes) %in% names(signnow))
      # check if system has failed with this component failure
      if (all(signnow$Probability == 0)){ # system has failed now, repair with cost cu
        failed <- TRUE
        cost <- cc
        taustarnow <- NA
        gotonext <- FALSE
      } else { # system has not failed, recalculate taustarnow at current tnow and compare with next grid time        
        gnowvec <- gnowhor(signnow, n0y0[wctnow], beta[wctnow], ftsnow[wctnow], tnow, hor = hor,
                           seqlen = seqlen, prior = prior, cu = cu, cp = cp, onecycle = onecycle)
        taustarnow <- gnowvec$tau[which.min(gnowvec$gnow)] - tnow
        if (taustarnow < nextgridtime) { 
          cost <- cc
          gotonext <- FALSE
        }
      }
    } else { # no failure before next tnow grid time
      
    }
  }



#