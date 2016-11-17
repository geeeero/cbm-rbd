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
sim1cycle <- function(sys, ctypes, compfts, n0y0, beta, tnowstep, hor, tprep, trepa = 0,
                      seqlen = 101, prior = FALSE, cu = 1, cp = 0.2, onecycle = TRUE){
  K <- length(ctypes)
  ftschron <- sort(unlist(compfts)) # when something fails
  # initializing: initial arguments for gnowhor
  tnow <- 0
  sysnow <- sys
  signnow <- computeSystemSurvivalSignature(sysnow)
  wctnow <- 1:K       # which component types are now present in the system?
  ftsnow <- as.list(rep(list(NULL), K))
  gotonext <- TRUE     # indicator that loop should go on
  failed <- FALSE      # indicator whether the system has failed
  repschedfor <- Inf   # for which time is repair scheduled?
  censtime <- Inf
  res <- data.frame() # initialize results data frame
  while(gotonext){
    cat("tnow =", tnow,"\n")
    if (!failed){ # if not failed already in previous loop, check...
      if (all(signnow$Probability == 0)){ # system has failed now, schedule repair for tnow + tprep if not scheduled already
        failed <- TRUE
        censtime <- tnow
        taustarnow <- NA
        if (repschedfor == Inf)
          repschedfor <- tnow + tprep
      } else { # system has not failed this time
        if (repschedfor == Inf){ # no repair scheduled yet
          # calculate current taustar
          gnowvec <- gnowhor(signnow, n0y0[wctnow], beta[wctnow], ftsnow[wctnow], tnow, hor = hor,
                             seqlen = seqlen, prior = prior, cu = cu, cp = cp, onecycle = onecycle)
          taustarnow <- gnowvec$tau[which.min(gnowvec$gnow)] - tnow
          if (taustarnow <= tprep){ # schedule repair for tnow + tprep if not scheduled already
            repschedfor <- tnow + tprep
            censtime <- tnow + tprep
          } # else do nothing & go on to next loop (repschedfor is still Inf)
        } else { # repair already scheduled: do nothing
          taustarnow <- NA
        }
      } # end "system has not failed this time"
    } else { # system was in failed state already at start of loop 
      taustarnow <- NA
    }
    # now we know if failed or not, if and when repair scheduled
    cat("failed =", failed, "taustar =", taustarnow, "repschedfor =", repschedfor, "\n")
    if (tnow >= repschedfor){ # time for repair has come, operational cycle ends
      gotonext <- FALSE
    }
    # write all current things in results data frame
    resnow <- data.frame(tnow = tnow, failed = failed, taustar = taustarnow, repschedfor = repschedfor, censtime = censtime)
    res <- rbind(res, resnow)
    # now prepare for next loop
    tnow <- tnow + tnowstep
    # update stuff if system not failed (otherwise not needed except tnow update!)
    if (!failed){
      failedcompsnow <- ftschron[ftschron <= tnow]
      if (length(failedcompsnow) > 0){ # update stuff only when the first component has failed
        # recalculates everything from scratch (do this in a more clever way?)
        sysnow <- induced_subgraph(sys, vids=V(sys)[!(name %in% names(failedcompsnow))])
        # TODO: catch error when sysnow now contains vertices s and t only (suddenly all components fail) 
        signnow <- computeSystemSurvivalSignature(sysnow)
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
      }
    } # end update if system not failed
  } # end while loop
  # update Weibull parameters (all component types!) for next operational cycle using nnynlist() (to be defined)
  # return res and all, tend = last.tnow + trepa, unit cost rate = (cp or cu) / last.tnow
  return(res) # TODO: return the other stuff as well
}  








#