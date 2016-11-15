# --------------------------------------------------------------------------- #
# CBM systems paper - functions for simulation
# --------------------------------------------------------------------------- #

# function to simulate a single operational cycle, i.e. until the system fails,
# or it is repaired preventively
# sys      system reliability block diagram at t = 0
# ctypes   list giving the types of components in the system, same format as
#          used for setCompTypes(), but needs the same order as in survsign table (!!!)
# compfts  list with elements named like the vertices in sys, giving the (simulated)
#          failure time of this component -- assumed that no failure times == 0
# n0y0     list of K prior parameter pairs c(n0,y0)
# beta     vector of K fixed weibull shape parameters
# tnowstep time interval after which current system reliability is re-evaluated
# hor      how many time units from tnow into the future to calculate sysrelnow
# tprep    preparation time to do preventive maintenance
# trepa    time needed to repair the system
# seqlen   at how many points to evaluate sysrelnow
# cu       cost of unplanned (corrective) repair action
# cp       cost of planned (preventive) repair action, cp < cu
# onecycle whether to use the one-cycle criterion or the renewal-based one
sim1cycle <- function(sys, ctypes, compfts, n0y0, beta, fts, tnowstep, hor, tprep, trepa = 0,
                      seqlen = 101, prior = FALSE, cu = 1, cp = 0.2, onecycle = TRUE){
  K <- length(ctypes)
  ftschron <- sort(unlist(compfts)) # when something fails
  # initializing: initial arguments for gnowhor
  tnow <- 0
  sysnow <- sys
  signnow <- computeSystemSurvivalSignature(sysnow)
  wctnow <- 1:K # which component types are now present in the system?
  ftsnow <- as.list(rep(list(NULL), K))
  gotonext = TRUE # indicator that loop should go on
  failed = FALSE # indicator whether the system has failed
  repschedfor = Inf # for which time is repair scheduled?
  # initialize results data frame
  res <- data.frame()
  while(gotonext){
    if (!failed & repschedfor == Inf){ # not failed and no repair scheduled yet !!!still need to check if system fails if repair scheduled!!!
      if (all(signnow$Probability == 0)){ # system has failed, schedule repair for tnow + tprep
        failed <- TRUE
        taustarnow <- NA
        repschedfor <- tnow + tprep
      } else { # system has not failed
        # calculate current taustar
        gnowvec <- gnowhor(signnow, n0y0[wctnow], beta[wctnow], ftsnow[wctnow], tnow, hor = hor,
                           seqlen = seqlen, prior = prior, cu = cu, cp = cp, onecycle = onecycle)
        taustarnow <- gnowvec$tau[which.min(gnowvec$gnow)]
        if (taustarnow <= tprep){ # schedule repair for tnow + tprep
          repschedfor <- tnow + tprep
        } # else do nothing & go on to next loop (repschedfor is +Inf)
      }
    } else { # system was either in failed state already at start of loop or has repair already scheduled !!!!!
      taustarnow <- NA
    }
    # now we know if failed or not, if and when repair scheduled
    if (tnow >= repschedfor){ # time for repair has come, operational cycle ends
      gotonext <- FALSE
    }
    # write all current things in results data frame
    resnow <- data.frame(tnow = tnow, failed = failed, taustar = taustarnow, repschedfor = repschedfor)
    res <- rbind(res, resnow)
    # now prepare for next loop
    tnow <- tnow + tnowstep
    # update stuff if not failed (otherwise not needed except tnow update!)
    if (!failed){
      failedcompsnow <- ftschron[ftschron <= tnow]
      sysnow <- induced_subsystem(sysnow, vids=V(sys)[!(name %in% names(failedcompsnow))])
      signnow <- computeSystemSurvivalSignature(sysnow)
      # to which fts list element do the failure times in failedcompsnow belong?
      ftsnowindex <- sapply(ctypes, function(ctypesl) names(failedcompsnow) %in% ctypesl)
      ftsnowindex <- apply(ftsnowindex, 1, which)
      # write failure times in corresponding list element
      for (k in 1:K)
        ftsnow[[k]] <- failedcompsnow[ftsnowindex == k]
      # which component types are now present? (subset of 1:K)
      wctnow <- which(names(ctypes) %in% names(signnow))
    }
  } # end while loop
  # update Weibull parameters (all component types!) for next operational cycle
  # return res and all, tend = last.tnow + trepa, unit cost rate = (cp or cu) / last.tnow
  return(res)
}  








#