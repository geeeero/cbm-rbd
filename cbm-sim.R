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
# thresh   threshold time below which to do preventive maintenance
# seqlen   at how many points to evaluate sysrelnow
# cu       cost of unplanned (corrective) repair action
# cp       cost of planned (preventive) repair action, cp < cu
# onecycle whether to use the one-cycle criterion or the renewal-based one
# timeround to how many digits to round the time points
sim1cycle <- function(sys, ctypes, compfts, n0y0, beta, tnowstep, hor, thresh, seqlen = 101,
                      prior = FALSE, cu = 1, cp = 0.2, onecycle = TRUE, timeround = 5){
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
  res <- data.frame() # initialize results data frame
  while(gotonext){
    cat("tnow =", tnow, ": ") #"failedcompsnow =", failedcompsnow, "\n")
    if (all(signnow$Probability == 0)){ # system has failed now, repair with cost cu
      failed <- TRUE
      cost <- cu
      taustarnow <- NA
      gotonext <- FALSE
    } else { # system has not failed: calculate taustarnow and check whether to repair preventively
      gnowvec <- gnowhor(signnow, n0y0[wctnow], beta[wctnow], ftsnow[wctnow], tnow, hor = hor,
                         seqlen = seqlen, prior = prior, cu = cu, cp = cp, onecycle = onecycle)
      taustarnow <- gnowvec$tau[which.min(gnowvec$gnow)] - tnow
      if (taustarnow <= thresh){ # repair preventively if taustar below theshold
        cost <- cp
        gotonext <- FALSE
      } # else do nothing & go on to next loop      
    } # end "system has not failed this time"
    # now we know if failed or not, if repair done
    cat("failed =", failed, "taustar =", taustarnow, "\n")
    # write all current things in results data frame
    resnow <- data.frame(tnow = tnow, failed = failed, taustar = taustarnow)
    res <- rbind(res, resnow)
    # now prepare for next loop: update stuff only if needed
    if (gotonext){
      tnow <- round(tnow + tnowstep, timeround)
      failedcompsnow <- ftschron[ftschron <= tnow]
      if (length(failedcompsnow) > 0){ # update stuff only when the first component has failed
        # recalculates everything from scratch (do this in a more clever way?)
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
      }
    } # end update if gotonext = TRUE
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
# thresh   threshold time below which to do preventive maintenance
# seqlen   at how many points to evaluate sysrelnow
# cu       cost of unplanned (corrective) repair action
# cp       cost of planned (preventive) repair action, cp < cu
# onecycle whether to use the one-cycle criterion or the renewal-based one
# timeround to how many digits to round the time points
simNcycle <- function(sys, ctypes, compfts, n0y0, beta, tnowstep, hor, thresh, seqlen = 101,
                      prior = FALSE, cycleupdate = TRUE, cu = 1, cp = 0.2, onecycle = TRUE, timeround = 5){
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
                          cu = cu, cp = cp, onecycle = onecycle, timeround = timeround)
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
# cu       cost of unplanned (corrective) repair action
# cp       cost of planned (preventive) repair action, cp < cu
# onecycle whether to use the one-cycle criterion or the renewal-based one
sim1cycleAgebased <- function(sys, ctypes, compfts, n0y0, beta, tnowstep, hor, seqlen = 101,
                              cu = 1, cp = 0.2, onecycle = TRUE){
  K <- length(ctypes)
  ftschron <- sort(unlist(compfts)) # when something fails
  ftschroni <- 1
  # get taustar: initial arguments for gnowhor
  cat("Calculating maintenance interval ... ")
  sign <- computeSystemSurvivalSignature(sys)
  fts0 <- as.list(rep(list(NULL), K))
  gvec <- gnowhor(sign, n0y0, beta, fts0, tnow = 0, hor = hor, seqlen = seqlen, cu = cu, cp = cp, onecycle = onecycle)
  taustar <- gvec$tau[which.min(gvec$gnow)]
  cat("taustar =", taustar, "\n")
  # set up loop over component failure times: assumes that system is not failed initially (before first component failure)
  gotonext <- TRUE
  failed <- FALSE
  res <- data.frame(tnow = 0, failed = FALSE) # initialize results data frame (assumed that initial system functions)
  while(gotonext){
    # when is the next component failure?
    failedcompsnow <- ftschron[1:ftschroni]
    tnow <- ceiling(ftschron[ftschroni] / tnowstep) * tnowstep
    # have any more components failed by this tnow?
    failedcompsnow <- ftschron[ftschron <= tnow]
    ftschroni <- length(failedcompsnow)
    if(taustar < tnow){ # taustar is before next failure: repair preventively
      tnow <- ceiling(taustar / tnowstep) * tnowstep # moment of repair in tnowstep resolution
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
      if (all(signnow$Probability == 0)){ # system has failed now, repair with cost cu
        failed <- TRUE
        gotonext <- FALSE
        cost <- cu
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
  if (failed) # last tnow minus tnowstep)
    tfunc <- tnow - tnowstep
  else # last tnow
    tfunc <- tnow
  # realized unit cost rate
  costrate <- cost / tnow
  # return res and all
  list(res = res, nnyn = nnyn, tfunc = tfunc, tend = tnow, costrate = costrate)
}  

# function to calculate tend, tfunc, costrate for age-based policy from simNcycle object
#simNcycleAgebased <- function(simNcycleobj, thresh, cu = 1, cp = 0.2){
#  
#}


# function to calculate tend, costrate for a single cycle using the corrective policy (cp is ignored)
# sys      system reliability block diagram at t = 0
# ctypes   list giving the types of components in the system, same format as
#          used for setCompTypes(), but needs the same order as in survsign table (!!!)
# compfts  list with elements named like the vertices in sys, giving the N (simulated)
#          failure times of this component -- assumed that no failure times == 0
# tnowstep time interval after which current system reliability is re-evaluated
# cu       cost of unplanned (corrective) repair action
# cp       cost of planned (preventive) repair action, cp < cu
sim1cycleCorrective <- function(sys, ctypes, compfts, tnowstep, cu = 1, cp = 0.2){
  K <- length(ctypes)
  ftschron <- sort(unlist(compfts)) # when something happens
  ftschroni <- 1
  gotonext <- TRUE
  res <- data.frame(tnow = 0, failed = FALSE) # initialize results data frame (assumed that initial system functions)
  while(gotonext){
    failedcompsnow <- ftschron[1:ftschroni]
    tnow <- ceiling(ftschron[ftschroni] / tnowstep) * tnowstep
    sysnow <- induced_subgraph(sys, vids=V(sys)[!(name %in% names(failedcompsnow))])
    # catch error when sysnow now contains vertices s and t only (suddenly all components fail)
    if(distances(sysnow, "s", "t") < Inf)
      signnow <- computeSystemSurvivalSignature(sysnow)
    else
      signnow <- data.frame(Probability = 0)
    # now check if system failed
    if (all(signnow$Probability == 0)){ # system has failed now, repair with cost cu
      failed <- TRUE
      gotonext <- FALSE
    } else { # system has not failed, go to next ftschron entry
      failed <- FALSE
      ftschroni <- ftschroni + 1
    }
    # write all current things in results data frame
    resnow <- data.frame(tnow = tnow, failed = failed)
    res <- rbind(res, resnow) 
  }
  names(tnow) <- NULL
  # time the system functioned ( = last tnow minus tnowstep)
  tfunc <- tnow - tnowstep
  # realized unit cost rate
  costrate <- cu / tnow
  # return res and all
  list(res = res, tfunc = tfunc, tend = tnow, costrate = costrate) 
}

# function to calculate tend, costrate for N operational cycles using the corrective policy (cp is ignored)
# sys      system reliability block diagram at t = 0
# ctypes   list giving the types of components in the system, same format as
#          used for setCompTypes(), but needs the same order as in survsign table (!!!)
# compfts  list with elements named like the vertices in sys, giving the N (simulated)
#          failure times of this component -- assumed that no failure times == 0
# tnowstep time interval after which current system reliability is re-evaluated
# cu       cost of unplanned (corrective) repair action
# cp       cost of planned (preventive) repair action, cp < cu
simNcycleCorrective <- function(sys, ctypes, compfts, tnowstep, cu = 1, cp = 0.2){
  N <- length(compfts[[1]])
  if (any(sapply(compfts, length) != N))
    stop("each element of compfts must contain the same number of failure times")
  compftsi <- lapply(compfts, function(x) x[1])
  res <- res2 <- list()
  for (i in 1:N){
    cat("Operational cycle", i, "\n")
    res[[i]] <- sim1cycleCorrective(sys = sys, ctypes = ctypes, compfts = compftsi, tnowstep = tnowstep,
                                    cu = cu, cp = cp)
    res2 <- rbind(res2, data.frame(cycle = i, res[[i]]$res))
    if (i < N){ # update stuff for next cycle
      compftsi <- lapply(compfts, function(x) x[i + 1])
    }
  }
  res2$cycle <- as.factor(res2$cycle)
  rownames(res2) <- NULL
  tfuncvec <- sapply(res, function(resi) resi$tfunc)
  tendvec <- sapply(res, function(resi) resi$tend)
  costratevec <- sapply(res, function(resi) resi$costrate)
  return(list(res = res2, tfunc = tfuncvec, tend = tendvec, costrate = costratevec))
}

# function to simulate Weibull failure times according to prior parameter choices (note different parametrization!)
# ncycles  how many failure times to simulate for each component
brWeibullData <- function(ncycles, beta, mttf){
  C1sim1 <- rweibull(ncycles, shape = beta[1], scale = (failuretolambda(mttf[1], beta[1]))^(1/beta[1]))
  C2sim1 <- rweibull(ncycles, shape = beta[1], scale = (failuretolambda(mttf[1], beta[1]))^(1/beta[1]))
  C3sim1 <- rweibull(ncycles, shape = beta[1], scale = (failuretolambda(mttf[1], beta[1]))^(1/beta[1]))
  C4sim1 <- rweibull(ncycles, shape = beta[1], scale = (failuretolambda(mttf[1], beta[1]))^(1/beta[1]))
  Hsim1  <- rweibull(ncycles, shape = beta[2], scale = (failuretolambda(mttf[2], beta[2]))^(1/beta[2]))
  Msim1  <- rweibull(ncycles, shape = beta[3], scale = (failuretolambda(mttf[3], beta[3]))^(1/beta[3]))
  P1sim1 <- rweibull(ncycles, shape = beta[4], scale = (failuretolambda(mttf[4], beta[4]))^(1/beta[4]))
  P2sim1 <- rweibull(ncycles, shape = beta[4], scale = (failuretolambda(mttf[4], beta[4]))^(1/beta[4]))
  P3sim1 <- rweibull(ncycles, shape = beta[4], scale = (failuretolambda(mttf[4], beta[4]))^(1/beta[4]))
  P4sim1 <- rweibull(ncycles, shape = beta[4], scale = (failuretolambda(mttf[4], beta[4]))^(1/beta[4]))
  list(C1 = C1sim1, C2 = C2sim1, C3 = C3sim1, C4 = C4sim1, H = Hsim1,
       M = Msim1, P1 = P1sim1, P2 = P2sim1, P3 = P3sim1, P4 = P4sim1)  
}


#