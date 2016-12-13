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
# thresh   deprecated, equal to tnowstep
# seqlen   at how many points to evaluate sysrelnow
# cc       cost of corrective repair action
# cp       cost of preventive repair action, cp < cc
# onecycle whether to use the one-cycle criterion or the renewal-based one
sim1cycle <- function(sys, ctypes, compfts, n0y0, beta, tnowstep, hor, thresh = tnowstep,
                      seqlen = 101, prior = FALSE, cc = 1, cp = 0.2, onecycle = TRUE){
  K <- length(ctypes)
  ftschron <- sort(unlist(compfts)) # when something fails
  # initializing: initial arguments for gnowhor
  tnow <- 0 # tnow: can be set to failure times
  cat("tnow =", tnow, ": ") #"failedcompsnow =", failedcompsnow, "\n")
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
                     seqlen = seqlen, prior = prior, cu = cc, cp = cp, onecycle = onecycle)
  taustarnow <- gnowvec$tau[which.min(gnowvec$gnow)] - tnow
  # initialize results data frame with entries for tnow = 0
  res <- data.frame(tnow = tnow, failed = failed, taustar = taustarnow)
  cat("failed =", failed, "taustar =", taustarnow, "\n")
  # stop already now if taustarnow smaller than tnowstep
  if (taustarnow < tnowstep){
    cost <- cp
    gotonext <- FALSE
  }
  while(gotonext){
    nextgridtime <- (m + 1) * tnowstep
    if (any((ftschron > tnow) & (ftschron <= nextgridtime))){ # failure before next tnow grid time
      tnow <- min(ftschron[ftschron > tnow])
      cat("tnow =", tnow, ": ") #"failedcompsnow =", failedcompsnow, "\n")
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
      if (all(signnow$Probability == 0)){ # system has failed now, corrective repair with cost cc
        failed <- TRUE
        cost <- cc
        taustarnow <- NA
        gotonext <- FALSE
      } else { # system has not failed, recalculate taustarnow at current tnow and compare with next grid time        
        gnowvec <- gnowhor(signnow, n0y0[wctnow], beta[wctnow], ftsnow[wctnow], tnow, hor = hor,
                           seqlen = seqlen, prior = prior, cu = cc, cp = cp, onecycle = onecycle)
        taustarnow <- gnowvec$tau[which.min(gnowvec$gnow)] - tnow
        if (taustarnow < (nextgridtime - tnow)) { # do preventive repair 
          cost <- cp
          gotonext <- FALSE
        } # else go to next while loop, but do not increment m as there could be other failures before next grid time
      } # end failure before next grid time
    } else { # no failure before next grid time: go to this time and calculate taustarnow
      tnow <- nextgridtime
      cat("tnow =", tnow, ": ") #"failedcompsnow =", failedcompsnow, "\n")
      # no need to update gnowvec inputs to current tnow except tnow itself,
      # as no new failure has happened, and can use previous gnowvec inputs
      gnowvec <- gnowhor(signnow, n0y0[wctnow], beta[wctnow], ftsnow[wctnow], tnow, hor = hor,
                         seqlen = seqlen, prior = prior, cu = cc, cp = cp, onecycle = onecycle)
      taustarnow <- gnowvec$tau[which.min(gnowvec$gnow)] - tnow
      if (taustarnow < tnowstep) { # do preventive repair 
        cost <- cp
        gotonext <- FALSE
      } else { # go to next while loop and increment m
        m <- m + 1
      }
    } # end no failure before next grid time
    # write things from current while loop into results data.frame
    resnow <- data.frame(tnow = tnow, failed = failed, taustar = taustarnow)
    res <- rbind(res, resnow)
    cat("failed =", failed, "taustar =", taustarnow, "\n")
  } # end while loop
  # update Weibull parameters (all component types!) for next operational cycle
  # get failure times from compftslist until tnow (time of planned / unplanned repair)
  ftslist <- as.list(rep(list(NULL), K))
  ftsfinal <- ftschron[ftschron <= tnow]
  ftsindex <- sapply(ctypes, function(ctypesl) names(ftsfinal) %in% ctypesl)
  if (length(ftsfinal) >= 1){
    if (length(ftsfinal) == 1){
      ftsindex <- which(ftsindex)
    } else {
      ftsindex <- apply(ftsindex, 1, which)
    }
    for (k in 1:K)
      ftslist[[k]] <- ftsfinal[ftsindex == k]
  }
  # censoring time is last tnow, number of censored components = all components - failed components
  censlist <- as.list(rep(list(NULL), K))
  ek <- sapply(ftslist, length)
  Nk <- apply(sign, 2, max)
  Nk <- Nk[-length(Nk)]
  ck <- Nk - ek
  for (k in 1:K)
    censlist[[k]] <- rep(tnow, ck[k])
  nnyn <- nnynlist(n0y0, ftslist, censlist, beta) # updated parameters at end of cycle
  # lower bound for the time the system has functioned ( = last tnow where not failed)
  # (tend = tnow gives time of end of cycle, functioning ended immediately before that)
  tfunc <- max(res$tnow[!(res$failed)])
  # realized unit cost rate
  costrate <- cost / tnow
  # return res and all
  list(res = res, nnyn = nnyn, tfunc = tfunc, tend = tnow, costrate = costrate)
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
# thresh   deprecated, equal to tnowstep
# seqlen   at how many points to evaluate sysrelnow
# cu       cost of unplanned (corrective) repair action
# cp       cost of planned (preventive) repair action, cp < cc
# onecycle whether to use the one-cycle criterion or the renewal-based one
simNcycle <- function(sys, ctypes, compfts, n0y0, beta, tnowstep, hor, thresh = tnowstep, seqlen = 101,
                      prior = FALSE, cycleupdate = TRUE, cc = 1, cp = 0.2, onecycle = TRUE){
  N <- length(compfts[[1]])
  if (any(sapply(compfts, length) != N))
    stop("each element of compfts must contain the same number of failure times")
  compftsi <- lapply(compfts, function(x) x[1])
  n0y0i <- n0y0
  res <- res2 <- nnyn <- list()
  for (i in 1:N){
    cat("Operational cycle", i, "\n")
    res[[i]] <- sim1cycle(sys = sys, ctypes = ctypes, compfts = compftsi, n0y0 = n0y0i, beta = beta,
                          tnowstep = tnowstep, hor = hor, thresh = thresh, seqlen = seqlen, prior = prior,
                          cc = cc, cp = cp, onecycle = onecycle)
    res2 <- rbind(res2, data.frame(cycle = i, res[[i]]$res))
    nnynlist <- list(nnyn, res[[i]]$nnyn)
    if (i < N){ # update stuff for next cycle
      compftsi <- lapply(compfts, function(x) x[i + 1])
      if (cycleupdate)
        n0y0i <- res[[i]]$nnyn
    }
  }
  res2$cycle <- as.factor(res2$cycle)
  tfuncvec <- sapply(res, function(resi) resi$tfunc)
  tendvec <- sapply(res, function(resi) resi$tend)
  costratevec <- sapply(res, function(resi) resi$costrate)
  return(list(res = res2, nnyn = nnynlist, tfunc = tfuncvec, tend = tendvec, costrate = costratevec))
}


# simulate one operational cycle with age-based maintenance policy
# sys      system reliability block diagram at t = 0
# ctypes   list giving the types of components in the system, same format as
#          used for setCompTypes(), but needs the same order as in survsign table (!!!)
# compfts  list with elements named like the vertices in sys, giving the (simulated)
#          failure time of this component -- assumed that no failure times == 0
# n0y0     list of K prior parameter pairs c(n0,y0)
# beta     vector of K fixed weibull shape parameters
# tnowstep time interval after which system functioning is re-evaluated
# hor      how many time units from tnow into the future to calculate sysrelnow
# seqlen   at how many points to evaluate sysrelnow
# cc       cost of unplanned (corrective) repair action
# cp       cost of planned (preventive) repair action, cp < cc
# onecycle whether to use the one-cycle criterion or the renewal-based one
# timeround to how many digits to round the time points
sim1cycleAgebased <- function(sys, ctypes, compfts, n0y0, beta, tnowstep, hor, seqlen = 101,
                              cc = 1, cp = 0.2, onecycle = TRUE, timeround = 5){
  K <- length(ctypes)
  ftschron <- sort(unlist(compfts)) # when something fails
  ftschroni <- 1
  # get taustar: initial arguments for gnowhor
  cat("Calculating maintenance interval ... ")
  sign <- computeSystemSurvivalSignature(sys)
  fts0 <- as.list(rep(list(NULL), K))
  gvec <- gnowhor(sign, n0y0, beta, fts0, tnow = 0, hor = hor, seqlen = seqlen, cu = cc, cp = cp, onecycle = onecycle)
  taustar <- gvec$tau[which.min(gvec$gnow)]
  cat("taustar =", taustar, "\n")
  # set up loop over component failure times: assumes that system is not failed initially (before first component failure)
  gotonext <- TRUE
  failed <- FALSE
  res <- data.frame(tnow = 0, failed = FALSE) # initialize results data frame (assumed that initial system functions)
  while(gotonext){
    # when is the next component failure?
    failedcompsnow <- ftschron[1:ftschroni]
    tnow <- ftschron[ftschroni]
    if(taustar < tnow){ # taustar is before next failure: repair preventively
      tnow <- taustar # moment of repair = taustar
      gotonext <- FALSE
      cost <- cp
    } else { # taustar is after next failure: check if system failed, if not go on
      sysnow <- induced_subgraph(sys, vids=V(sys)[!(name %in% names(failedcompsnow))])
      # catch error when sysnow now contains vertices s and t only (suddenly all components fail)
      if(distances(sysnow, "s", "t") < Inf)
        signnow <- computeSystemSurvivalSignature(sysnow)
      else
        signnow <- data.frame(Probability = 0)
      # now check if system failed
      if (all(signnow$Probability == 0)){ # system has failed now, repair with cost cc
        failed <- TRUE
        gotonext <- FALSE
        cost <- cc
      } else { # system has not failed, go to next ftschron entry
        ftschroni <- ftschroni + 1
      }      
    }
    cat("tnow =", tnow, ": failed =", failed, "\n")
    # write all current things in results data frame
    resnow <- data.frame(tnow = tnow, failed = failed)
    res <- rbind(res, resnow) 
  } # end while loop
  # update Weibull parameters (all component types!) for next operational cycle
  # get failure times from compftslist until tnow (time of planned / unplanned repair)
  ftslist <- as.list(rep(list(NULL), K))
  ftsfinal <- ftschron[ftschron <= tnow]
  ftsindex <- sapply(ctypes, function(ctypesl) names(ftsfinal) %in% ctypesl)
  if (length(ftsfinal) >= 1){
    if (length(ftsfinal) == 1){
      ftsindex <- which(ftsindex)
    } else {
      ftsindex <- apply(ftsindex, 1, which)
    }
    for (k in 1:K)
      ftslist[[k]] <- ftsfinal[ftsindex == k]
  }
  # censoring time is last tnow, number of censored components from
  censlist <- as.list(rep(list(NULL), K))
  ek <- sapply(ftslist, length)
  Nk <- apply(sign, 2, max)
  Nk <- Nk[-length(Nk)]
  ck <- Nk - ek
  for (k in 1:K)
    censlist[[k]] <- rep(tnow, ck[k])
  nnyn <- nnynlist(n0y0, ftslist, censlist, beta) # updated parameters at end of cycle
  # time the system functioned ( = last tnow where not failed)
  names(tnow) <- NULL
  # time the system functioned
  if (failed) # last time on grid before tnow
    tfunc <- floor(round((tnow / tnowstep), timeround)) * tnowstep
  else # last tnow
    tfunc <- tnow
  # realized unit cost rate
  costrate <- cost / tnow
  # return res and all
  list(res = res, nnyn = nnyn, tfunc = tfunc, tend = tnow, costrate = costrate)
}  


# simulate N operational cycles with age-based maintenance policy
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
# seqlen   at how many points to evaluate sysrelnow
# cc       cost of unplanned (corrective) repair action
# cp       cost of planned (preventive) repair action, cp < cc
# onecycle whether to use the one-cycle criterion or the renewal-based one
# timeround to how many digits to round the time points
simNcycleAgebased <- function(sys, ctypes, compfts, n0y0, beta, tnowstep, hor, seqlen = 101, 
                              cycleupdate = TRUE, cc = 1, cp = 0.2, onecycle = TRUE, timeround = 5){
  N <- length(compfts[[1]])
  if (any(sapply(compfts, length) != N))
    stop("each element of compfts must contain the same number of failure times")
  compftsi <- lapply(compfts, function(x) x[1])
  n0y0i <- n0y0
  res <- res2 <- nnyn <- list()
  for (i in 1:N){
    cat("Operational cycle", i, "\n")
    res[[i]] <- sim1cycleAgebased(sys = sys, ctypes = ctypes, compfts = compftsi, n0y0 = n0y0i, beta = beta, tnowstep = tnowstep,
                                  hor = hor, seqlen = seqlen, cc = cc, cp = cp, onecycle = onecycle, timeround = timeround)
    res2 <- rbind(res2, data.frame(cycle = i, res[[i]]$res))
    nnynlist <- list(nnyn, res[[i]]$nnyn)
    if (i < N){ # update stuff for next cycle
      compftsi <- lapply(compfts, function(x) x[i + 1])
      if (cycleupdate)
        n0y0i <- res[[i]]$nnyn
    }
  }
  res2$cycle <- as.factor(res2$cycle)
  tfuncvec <- sapply(res, function(resi) resi$tfunc)
  tendvec <- sapply(res, function(resi) resi$tend)
  costratevec <- sapply(res, function(resi) resi$costrate)
  return(list(res = res2, nnyn = nnynlist, tfunc = tfuncvec, tend = tendvec, costrate = costratevec))
}



#