# --------------------------------------------------------------------------- #
# CBM systems paper - functions for simulation with preparation time
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
# timeround to how many digits to round the time points
sim1cycleprep <- function(sys, ctypes, compfts, n0y0, beta, tnowstep, hor, tprep, trepa = 0,
                          seqlen = 101, prior = FALSE, cu = 1, cp = 0.2, onecycle = TRUE, timeround = 5){
  K <- length(ctypes)
  ftschron <- sort(unlist(compfts)) # when something fails
  # initializing: initial arguments for gnowhor
  tnow <- 0
  sysnow <- sys
  signnow <- sign <- computeSystemSurvivalSignature(sysnow)
  wctnow <- 1:K       # which component types are now present in the system?
  ftsnow <- as.list(rep(list(NULL), K))
  failedcompsnow <- numeric(0)
  gotonext <- TRUE     # indicator that loop should go on
  failed <- FALSE      # indicator whether the system has failed
  repschedfor <- Inf   # for which time is repair scheduled?
  res <- data.frame() # initialize results data frame
  while(gotonext){
    cat("tnow =", tnow, ":") #"failedcompsnow =", failedcompsnow, "\n")
    if (!failed){ # if not failed already in previous loop, check...
      if (all(signnow$Probability == 0)){ # system has failed now, schedule repair for tnow + tprep if not scheduled already
        failed <- TRUE
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
    resnow <- data.frame(tnow = tnow, failed = failed, taustar = taustarnow, repschedfor = repschedfor)
    res <- rbind(res, resnow)
    # now prepare for next loop
    tnow <- round(tnow + tnowstep, timeround)
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
  # update Weibull parameters (all component types!) for next operational cycle
  # get failure times from compftslist until time of repair repschedfor
  ftslist <- as.list(rep(list(NULL), K))
  ftsfinal <- ftschron[ftschron <= repschedfor]
  ftsindex <- sapply(ctypes, function(ctypesl) names(ftsfinal) %in% ctypesl)
  if (length(ftsfinal) == 1)
    ftsindex <- which(ftsindex)
  else
    ftsindex <- apply(ftsindex, 1, which)
  for (k in 1:K)
    ftslist[[k]] <- ftsfinal[ftsindex == k]
  # censoring time is last repschedfor, number of censored components from
  censlist <- as.list(rep(list(NULL), K))
  ek <- sapply(ftslist, length)
  Nk <- apply(sign, 2, max)
  Nk <- Nk[-length(Nk)]
  ck <- Nk - ek
  for (k in 1:K)
    censlist[[k]] <- rep(repschedfor, ck[k])
  nnyn <- nnynlist(n0y0, ftslist, censlist, beta) # updated parameters at end of cycle
  # time the whole cycle took
  tend <- repschedfor + trepa
  # time the system functioned
  tfunc <- max(res$tnow[!(res$failed)])
  # downtime and unit cost rate
  if (failed){
    downtime <- repschedfor - min(res$tnow[res$failed]) # only in tnowstep resolution!
    costrate <- cu / tend # *** or / tfunc? or penalty cost for downtime?
  } else {
    downtime <- 0
    costrate <- cp / tend # *** or / tfunc? or penalty cost for downtime?
  }
  # return res and all
  list(res = res, nnyn = nnyn, tend = tend, tfunc = tfunc, downtime = downtime, costrate = costrate)
}  


# function to simulate N operational cycles for one machine,
# i.e. updated parameters from one cycle are used as prior parameters in next cycle
# sys      system reliability block diagram at t = 0
# ctypes   list giving the types of components in the system, same format as
#          used for setCompTypes(), but needs the same order as in survsign table (!!!)
# compfts  list with elements named like the vertices in sys, giving the N (simulated)
#          failure times of this component -- assumed that no failure times == 0
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
# timeround to how many digits to round the time points
simNcycleprep <- function(sys, ctypes, compfts, n0y0, beta, tnowstep, hor, tprep, trepa = 0, seqlen = 101,
                          prior = FALSE, cycleupdate = TRUE, cu = 1, cp = 0.2, onecycle = TRUE, timeround = 5){
  N <- length(compfts[[1]])
  if (any(sapply(compfts, length) != N))
    stop("each element of compfts must contain the same number of failure times")
  compftsi <- lapply(compfts, function(x) x[1])
  n0y0i <- n0y0
  res <- res2 <- nnyn <- list()
  for (i in 1:N){
    cat("Operational cycle", i, "\n")
    res[[i]] <- sim1cycleprep(sys = sys, ctypes = ctypes, compfts = compftsi, n0y0 = n0y0i, beta = beta,
                              tnowstep = tnowstep, hor = hor, tprep = tprep, trepa = trepa, seqlen = seqlen,
                          prior = prior, cu = cu, cp = cp, onecycle = onecycle, timeround = timeround)
    res2 <- rbind(res2, data.frame(cycle = i, res[[i]]$res))
    nnynlist <- list(nnyn, res[[i]]$nnyn)
    if (i < N){ # update stuff for next cycle
      compftsi <- lapply(compfts, function(x) x[i + 1])
      if (cycleupdate)
        n0y0i <- res[[i]]$nnyn
    }
  }
  res2$cycle <- as.factor(res2$cycle)
  tendvec <- sapply(res, function(resi) resi$tend)
  tfuncvec <- sapply(res, function(resi) resi$tfunc)
  downtimevec <- sapply(res, function(resi) resi$downtime)
  costratevec <- sapply(res, function(resi) resi$costrate)
  return(list(res = res2, nnyn = nnynlist, tend = tendvec, tfunc = tfuncvec, downtime = downtimevec, costrate = costratevec))
}


#