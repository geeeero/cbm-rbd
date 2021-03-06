###############################################################################
#       Functions to calculate and plot system reliability functions          #
#       for systems consisting of Weibull components with fixed shape         #
#       Author: Gero Walter (g.m.walter@tue.nl, gero@gmx.at)                  #
###############################################################################

# load libraries if not yet loaded
# library(actuar) # provides dinvgamma(x, shape, scale) where shape = alpha and scale = beta
# library(ReliabilityTheory)
# library(ggplot2)
# library(reshape2)
# library(gridExtra)

# translate lambda parameter of weibull to expected failure time
# lambda - scale parameter of weibull
# beta   - shape parameter of weibull
lambdatofailure <- function(lambda, beta = 2){
  lambda^(1/beta) * gamma(1 + 1/beta)
}

# translate failure time to lambda parameter of weibull
# ft   - expected failure time
# beta - shape parameter of weibull
failuretolambda <- function(ft, beta = 2){
  (ft/gamma(1 + 1/beta))^beta
}

# inverse gamma pdf and cdf with canonical parameters n0 and y0
cdinvgamma <- function(x, n0, y0, ...)
  dinvgamma(x, shape = n0+1, scale = n0*y0, ...)

cpinvgamma <- function(x, n0, y0, ...)
  pinvgamma(x, shape = n0+1, scale = n0*y0, ...)

cpigforoptim <- function(ny, t, ...)
  cpinvgamma(t, ny[1], ny[2], ...)

cpigfornoptim <- function(n, y, t, ...)
  cpinvgamma(t, n, y, ...)

# function to calculate P(C_t = l | n^(0), y^(0), data),
# the posterior predictive probability that l components function at time t
# in the system observed until t_now (notion "of type k" dropped here)
# n0y0    pair c(n0,y0) of prior parameters (start of cycle values)
# t       time t for which to calculate P(C_t), t > t_now
# l       number of functioning components, \in {0, 1, ..., c = n-e}
# n       number of components in the system
# tnow    time until the system is observed
# fts     vector of length e giving the observed failure times,
#         or NULL if no failures observed
# beta    fixed weibull shape parameter
# prior   whether the prior predictive should be calculated instead
#         of the posterior predictive; when TRUE, fts is ignored 
postpredC <- function(n0y0, beta, n, fts, tnow, t, l, prior = FALSE){
  if (t < tnow)
    stop("t must be larger than tnow")
  e <- length(fts)
  if (n < e)
    stop("too many elements in fts, there can be at most n failure times!")
  if (l > n-e | l < 0)
    stop("l must be in [0, n-e]")
  if(!prior){
    nn <- n0y0[1] + e
    nnyn <- n0y0[1]*n0y0[2] + (n-e)*(tnow^beta) + sum(fts^beta)
  } else {
    nn <- n0y0[1]
    nnyn <- n0y0[1]*n0y0[2]
  }
  j <- seq(0, n-e-l)
  choose(n-e, l) * sum( (-1)^j * choose(n-e-l, j) * (nnyn/(nnyn + (l+j)*(t^beta - tnow^beta)))^(nn + 1) )
}

# calculates the system reliability for fixed future time t > tnow
# n0y0, beta, fts may contain only types that are still in the system (according to survsign),
# so K is here only the number of types still in the system at time tnow
# type order in n0y0, beta, fts must be the same as in (reduced) survsign table!
#
# survsign data frame with the survival signature as output by computeSystemSurvivalSignature()
#          this must be signature for the reduced system if components have failed!
#          -> use induced_subgraph()
# n0y0     list of K prior parameter pairs c(n0,y0)
# beta     vector of K fixed weibull shape parameters
# fts      list of K vectors giving the observed failure times for the compents of type 1,...,K;
#          the list element should be NULL if no failure has been observed for type k.
#          All list element entries must be <= tnow
# tnow     time until the system is observed
# t        time t for which to calculate P(T_sys > t), t > t_now
# table    if table of posterior predictive probabilities should be given along with the reliability
# prior    whether the prior system relibability should be calculated (= no updating with things observed until tnow)
sysrelnow <- function(survsign, n0y0, beta, fts, tnow, t, prior = FALSE, table = FALSE){
  K <- length(n0y0)
  #knames <- names(n0y0)
  # ck from survsign (number of functioning components per type)
  ck <- apply(survsign, 2, max)
  ck <- ck[-length(ck)]
  # e_k from fts (number of failed components per type)
  ek <- sapply(fts, length)
  Nk <- ck + ek
  if(length(fts) != K | length(beta) != K | length(Nk) != K)
    stop("Check the length of arguments n0y0, beta, fts.\n
         They need to contain the same number of elements as there are component types in survsign.")
  # keep only those rows where survsign > 0
  if (dim(survsign)[1] > 1) # but only if there is more than one row left
    survsign <- survsign[survsign$Probability > 0,]
  # for each component type, create a column containing 
  # the component posterior predictive probabilities
  for (k in 1:K){
    PCtk <- sapply(survsign[,k], function(l){
      postpredC(n0y0 = n0y0[[k]], beta = beta[k], n = Nk[k],
                fts = fts[[k]], tnow = tnow, t = t, l = l, prior = prior)
    })
    ssnames <- names(survsign)
    survsign <- data.frame(survsign, PCtk)
    names(survsign) <- c(ssnames, paste("PCt", ssnames[k], sep=""))
  }
  # joint posterior predictive probabilities
  survsign$summand <- apply(survsign[,-(1:K)], 1, prod)
  res <- sum(survsign$summand)
  if(table)
    return(list(rel = res, table = survsign))
  else
    return(res)
}

# calculates the system reliability for a vector of future times t > tnow,
# resp. hor time units into the future from tnow evaluated at seqlen points
# n0y0, beta, fts may contain only types that are still in the system (according to survsign),
# so K is here only the number of types still in the system at time tnow
# type order in n0y0, beta, fts must be the same as in (reduced) survsign table!
#
# survsign data frame with the survival signature as output by computeSystemSurvivalSignature()
#          this must be signature for the reduced system if components have failed!
#          -> use induced_subgraph()
# n0y0     list of K prior parameter pairs c(n0,y0)
# beta     vector of K fixed weibull shape parameters
# fts      list of K vectors giving the observed failure times for the compents of type 1,...,K;
#          the list element should be NULL if no failure has been observed for type k.
#          All list element entries must be <= tnow
# tnow     time until the system is observed
# prior    whether the prior system relibability should be calculated (= no updating with things observed until tnow)
# tvec     vector of times t for which to calculate P(T_sys > t), t > t_now
# hor      how many time units from tnow into the future to calculate sysrelnow
# seqlen   at how many points to evaluate sysrelnow
sysrelnowvec <- function(survsign, n0y0, beta, fts, tnow, tvec, prior = FALSE){
  if(min(tvec) < tnow)
    stop("Vector of times to evaluate sysrelnow must contain values > tnow")
  sapply(tvec, function(t){
    sysrelnow(survsign = survsign, n0y0 = n0y0, beta = beta, fts = fts, tnow = tnow, t = t, prior = prior)
  })
}

sysrelnowhor <- function(survsign, n0y0, beta, fts, tnow, hor, seqlen=101, prior = FALSE){
  tvec <- seq(tnow, tnow + hor, length.out = seqlen)
  rel <- sapply(tvec, function(t){
    sysrelnow(survsign = survsign, n0y0 = n0y0, beta = beta, fts = fts, tnow = tnow, t = t, prior = prior)
  })
  data.frame(t = tvec, rel = rel)
}


# calculate sysrelnow at repeated times tnow for given failure history
#
# sys      system reliability block diagram
# ctypes   list giving the types of components in the system, same format as
#          used for setCompTypes(), but needs the same order as in survsign table (!!!)
# compfts  list with elements named like the vertices in sys, giving the failure
#          time of this component; can be NA if component does not fail
# n0y0     list of K prior parameter pairs c(n0,y0)
# beta     vector of K fixed weibull shape parameters
# tnowvec  vector of times for which to evaluate sysrelnow
# hor      how many time units from tnow into the future to calculate sysrelnow
# seqlen   at how many points to evaluate sysrelnow
# prior    whether the prior system relibability should be calculated (= no updating with things observed until tnow)
# verbose  if true, prints some progress information to the console
sysrelnowhist <- function(sys, ctypes, compfts, n0y0, beta, tnowvec, hor, seqlen = 101, prior = FALSE, verbose = TRUE){
  K <- length(ctypes)
  # calculate survsign only when something fails
  # create fts only the same number of times
  ftschron <- c(0, sort(unlist(compfts))) # when something happens
  # create list with ftschron-th entry the survsign / fts valid from time in corresponding ftschron element;
  # first element corresponds to situation at system startup
  siglist <- as.list(rep(list(NULL), length(ftschron)))
  siglist[[1]] <- computeSystemSurvivalSignature(sys)
  ftslist <- as.list(rep(NA, length(ftschron)))
  ftslist[[1]] <- as.list(rep(list(NULL), K))
  # which component types are still in the system?
  presentctypes <- as.list(rep(NA, length(ftschron)))
  presentctypes[[1]] <- 1:K
  # now go through ftschron (first element already done) and calculate survsign and fts
  failedcomps <- numeric(0)
  if (verbose)
    cat("calculating survival signature(s) ...\n")
  for (i in 2:length(ftschron)){
    failedcomps <- c(failedcomps, names(ftschron[i]))
    reducedsys <- induced_subgraph(sys, vids=V(sys)[!(name %in% failedcomps)])
    siglist[[i]] <- computeSystemSurvivalSignature(reducedsys)
    ftslisti <- ftslist[[i-1]]
    ftslistiindex <- which(sapply(ctypes, function(ctypesl) names(ftschron[i]) %in% ctypesl))
    ftslisti[[ftslistiindex]] <- c(ftslisti[[ftslistiindex]], ftschron[i])
    ftslist[[i]] <- ftslisti
    presentctypes[[i]] <- which(names(ctypes) %in% names(siglist[[i]]))
  }
  #list(siglist=siglist, ftslist=ftslist, prctypes = presentctypes)
  # now get reliability curves for each tnow
  #res <- data.frame()
  reslen <- length(tnowvec) * seqlen
  res <- data.frame(tnow = rep(tnowvec, each = seqlen), t = rep(NA, reslen), rel = rep(NA, reslen))
  for (tnow in tnowvec){
    # index in siglist, ftslist, presentctypes for current (tnow) status
    tnowi <- max(which(ftschron <= tnow))
    prctps <- presentctypes[[tnowi]]
    if (verbose)
      cat("calculating sysrel for tnow =", tnow, "...\n")
    current <- data.frame(tnow = tnow,
                          # add nn, ny
                          sysrelnowhor(survsign = siglist[[tnowi]], n0y0 = n0y0[prctps],
                                       beta = beta[prctps], fts = ftslist[[tnowi]][prctps],
                                       tnow = tnow, hor = hor, seqlen = seqlen, prior = prior))
    #res <- rbind(res, current)
    res[res$tnow == tnow,] <- current
  }
  return(res)
}

# transform failure history argument for sysrelnowhist() into data.frame for plotting the history
# compfts  list with elements named like the vertices in sys, giving the failure
#          time of this component; can be NA if component does not fail
# maxtnow  censoring time for the non-failed components;
#          equal to last failure time if not given.
#          If set smaller than last failure time, maxtnow is the censoring time
#          for all failure times after maxtnow
compfts2df <- function(compfts, maxtnow = NA){
  compftsdf <- t(as.data.frame(compfts))
  compftsdf <- data.frame(compftsdf, Components = names(compfts))
  names(compftsdf)[1] <- "t"
  compftsdf <- data.frame(compftsdf, cens = is.na(compftsdf$t), id = 1:length(compfts))
  if(is.na(maxtnow)){
    maxtnow <- max(na.omit(compftsdf$t))
    compftsdf$t[is.na(compftsdf$t)] <- maxtnow
  } else {
    if(maxtnow >= max(na.omit(compftsdf$t))){
      compftsdf$t[is.na(compftsdf$t)] <- maxtnow
    } else {
      compftsdf$t[is.na(compftsdf$t)] <- maxtnow
      compftsdf$cens[compftsdf$t > maxtnow] <- TRUE
      compftsdf$t[compftsdf$t > maxtnow] <- maxtnow
    }
  }
  return(compftsdf)
}

# creates a ggplot2 plot object for plotting the failure history.
# To get the plot use print(obj)
# compfts  list with elements named like the vertices in sys, giving the failure
#          time of this component; can be NA if component does not fail
# maxtnow  censoring time for the non-failed components;
#          equal to last failure time if not given.
#          If set smaller than last failure time, maxtnow is the censoring time
#          for all failure times after maxtnow
plotfts <- function(compfts, maxtnow = NA, pointsize = 1, tlabelvjust = -1, labelround = 2){
  compftsdf <- compfts2df(compfts = compfts, maxtnow = maxtnow)
  ggplot(compftsdf, aes(x = Components, y = t, ymin = rep(0, dim(compftsdf)[1]), ymax = t)) + geom_linerange() +
    geom_pointrange(aes(shape = factor(cens)), size = pointsize) + geom_text(aes(label = round(t, labelround), vjust = tlabelvjust)) + 
    coord_flip() + scale_shape_manual(values = c(19, 1), label = c("Yes", "No"), name = "Failure")
}



# calculate g^(now)(tau) for given vectors of tau, R_sys and f_sys
# tau      vector of time points (prospective time, needs to start close to 0)
# rel      vector of reliability function values corresponding to tau
# f        vector of density function values corresponding to tau
# cu       cost of unplanned (corrective) repair action
# cp       cost of planned (preventive) repair action, cp < cu
# onecycle which one-cycle cost criterion to use
#          (was: whether to use the one-cycle criterion or the renewal-based one)
# tnow     one-cycle total cost criterion needs current tnow
gnow <- function(tau, rel, f, cu = 1, cp = 0.2, onecycle = TRUE, tnow = 0){
  if (onecycle){
    # one-cycle prospective cost criterion 
    gnow <- (cp / tau) * rel + cu * sapply(tau, function(ctau){
      sum(f[tau <= ctau] / tau[tau <= ctau])
    })
  } else {
    # one-cycle total cost criterion 
    gnow <- (cp / (tau + tnow)) * rel + cu * sapply(tau, function(ctau){
      sum(f[tau <= ctau] / (tau[tau <= ctau] + tnow))
    })
    # renewal reward criterion
    #gnow <- (cp * rel + cu * (1 - rel)) / ( tau * rel + sapply(tau, function(ctau){
    #  sum(rel[tau <= ctau])
    #}))
  }
  data.frame(tau = tau, rel = rel, gnow = gnow)  
}


# calculate g^(now)(tau) for current tnow based on sysrelnowhor()
# output is in prospective time, i.e. g(tau) cives the cost rate
# for repair action in tau time units from now.
# takes same arguments as sysrelnowhor(), and additionally
# cu       cost of unplanned (corrective) repair action
# cp       cost of planned (preventive) repair action, cp < cu
# onecycle whether to use the one-cycle criterion or the renewal-based one
gnowhor <- function(survsign, n0y0, beta, fts, tnow, hor, seqlen=101, prior = FALSE,
                    cu = 1, cp = 0.2, onecycle = TRUE){
  trel <- sysrelnowhor(survsign = survsign, n0y0 = n0y0, beta = beta, fts = fts,
                       tnow = tnow, hor = hor, seqlen = seqlen, prior = prior)
  f <- -diff(trel$rel)
  tau <- trel$t[-1] - tnow # shift to start directly after 0
  rel <- trel$rel[-1]
  res <- gnow(tau = tau, rel = rel, f = f, cu = cu, cp = cp, onecycle = onecycle, tnow = tnow)
  res$tau <- tau + tnow
  return(res)
}



# calculate gnow(tau) at repeated times tnow for given failure history
# uses sysrelnowhist()
#
# sys      system reliability block diagram
# ctypes   list giving the types of components in the system, same format as
#          used for setCompTypes(), but needs the same order as in survsign table (!!!)
# compfts  list with elements named like the vertices in sys, giving the failure
#          time of this component; can be NULL if component does not fail
# n0y0     list of K prior parameter pairs c(n0,y0)
# beta     vector of K fixed weibull shape parameters
# tnowvec  vector of times for which to evaluate sysrelnow
# hor      how many time units from tnow into the future to calculate sysrelnow
# seqlen   at how many points to evaluate sysrelnow
# cu       cost of unplanned (corrective) repair action
# cp       cost of planned (preventive) repair action, cp < cu
# onecycle whether to use the one-cycle criterion or the renewal-based one
# verbose  if true, prints some progress information to the console
gnowhist <- function(sys, ctypes, compfts, n0y0, beta, tnowvec, hor, seqlen = 101, prior = FALSE,
                     cu = 1, cp = 0.2, onecycle = TRUE, verbose = TRUE){  
  sysrels <- sysrelnowhist(sys = sys, ctypes = ctypes, compfts = compfts, n0y0 = n0y0, beta = beta,
                          tnowvec = tnowvec, hor = hor, seqlen = seqlen, prior = prior, verbose = verbose)
  # output: length(tnowvec) * seqlen x 3 data.frame with vars tnow, t, rel
  tnowlen <- length(tnowvec)
  gnowvec <- numeric(0)
  if (verbose)
    cat("calculating gnow function(s) ...\n")
  for (tnowi in tnowvec){
    trel <- subset(sysrels, tnow == tnowi)
    fi <- -diff(trel$rel)
    taui <- trel$t[-1] - tnowi # shift to start directly after 0
    reli <- trel$rel[-1]
    gnowvec <- c(gnowvec, NA, gnow(tau = taui, rel = reli, f = fi, cu = cu, cp = cp,
                                    onecycle = onecycle, tnow = tnowi)$gnow)
  }
  data.frame(sysrels, gnow = gnowvec)
}
  
# calculate optimal moment of maintenance history for given failure history
# uses gnowhist(), has same arguments and one extra argument (tool)
# sys      system reliability block diagram
# ctypes   list giving the types of components in the system, same format as
#          used for setCompTypes(), but needs the same order as in survsign table (!!!)
# compfts  list with elements named like the vertices in sys, giving the failure
#          time of this component; can be NULL if component does not fail
# n0y0     list of K prior parameter pairs c(n0,y0)
# beta     vector of K fixed weibull shape parameters
# tnowvec  vector of times for which to evaluate sysrelnow
# hor      how many time units from tnow into the future to calculate sysrelnow
# seqlen   at how many points to evaluate sysrelnow
# cu       cost of unplanned (corrective) repair action
# cp       cost of planned (preventive) repair action, cp < cu
# onecycle whether to use the one-cycle criterion or the renewal-based one
# verbose  if true, prints some progress information to the console
# tool     if true, returns also the gnowhist output, defaults to FALSE
taustarhist <- function(sys, ctypes, compfts, n0y0, beta, tnowvec, hor, seqlen = 101, prior = FALSE,
                        cu = 1, cp = 0.2, onecycle = TRUE, verbose = TRUE, tool = FALSE){
  res <- gnowhist(sys = sys, ctypes = ctypes, compfts = compfts, n0y0 = n0y0, beta = beta, tnowvec = tnowvec,
                  hor = hor, seqlen = seqlen, prior = prior, cu = cu, cp = cp, onecycle = onecycle, verbose = verbose)
  res2 <- subset(res, !is.na(res$gnow))
  cstars <- aggregate(gnow ~ tnow, res2, min)$gnow
  tstarsi <- aggregate(gnow ~ tnow, res2, which.min)
  tstarsi <- tstarsi$gnow + (seqlen - 1) * (0:(length(tnowvec) - 1))
  tstars <- res2$t[tstarsi]
  taustars <- tstars - tnowvec
  relstars <- res2$rel[tstarsi]
  cuint <- sapply(tnowvec, function(tnowi){
    trel <- subset(res, res$tnow == tnowi)
    f <- -diff(trel$rel)
    tauvec <- trel$t[-1] # absolute timescale!
    tau <- trel$t[which.min(trel$gnow)]
    sum(f[tauvec <= tau] / tauvec[tauvec <= tau])
  })
  #cat("cuint =", cuint)
  ctotal <- cp / tstars * res2$rel[tstarsi] + cu * cuint
  taus <- data.frame(tnow = tnowvec, taustar = taustars, tstar = tstars,
                     cstar = cstars, ctotal = ctotal, relstar = relstars)
  if (!tool)
    return(taus)
  else
    return(list(relg = res,
                taus = taus))
}

# function to process the output of taustarhist for use in the tool (figure plotting)
# taustarhistobj  object returned by taustarhist(..., tool = TRUE)
taustarhistprocessor <- function(taustarhistobj){
  taustarhistobj$relg$tnow <- as.factor(taustarhistobj$relg$tnow)
  names(taustarhistobj$relg)[c(3,4)] <- c("Reliability", "Costrate")
  #taustarhistobj$relg <- melt(taustarhistobj$relg, c("tnow", "t"))
  #taustarhistobj$relg <- subset(taustarhistobj$relg, !is.na(value))
  taustarhistobj$taus$tnowf <- factor(taustarhistobj$taus$tnow)
  taustarhistobj
}

# function to plot sysrel, gnow and taustar based on (unprocessed) taustarhist object
# taustarhistobj  object returned by taustarhist(..., tool = TRUE)
# gquantile       quantile to define the upper plotting bound for the g function plot
# taulabelvjust   postition adjustment for the numbers appearing on the tau plot
# taulabelround   number of digits for the numbers appearing on the tau plot
toolplot <- function(taustarhistobj, gquantile = 0.9, taulabelhjust = -1, taulabelround = 3){
  taustarhistobj <- taustarhistprocessor(taustarhistobj)
  relfig <- ggplot(taustarhistobj$relg, aes(x = t, y = Reliability)) +
    geom_line(aes(group = tnow, colour = tnow)) + # xlab(expression(t)) +
    coord_cartesian(ylim = c(0,1)) + scale_colour_discrete(drop = FALSE) +
    guides(colour=guide_legend(title=expression(t[now])))
  gfigylims <- c(quantile(taustarhistobj$relg$Costrate, 0, na.rm = TRUE),
                 quantile(taustarhistobj$relg$Costrate, gquantile, na.rm = TRUE))
  gfig <- ggplot(subset(taustarhistobj$relg, !is.na(Costrate)), aes(x = t, y = Costrate)) +
    geom_line(aes(group = tnow, colour = tnow)) + ylab("Cost rate") +
    coord_cartesian(ylim = gfigylims) + scale_colour_discrete(drop = FALSE) +
    guides(colour=guide_legend(title=expression(t[now])))
  taufigxlimsfactor <- diff(range(taustarhistobj$taus$tnow)) * 0.2
  taufigxlims <- range(taustarhistobj$taus$tnow) + c(-1 * taufigxlimsfactor, taufigxlimsfactor)
  taufigylims <- c(min(taustarhistobj$relg$t), max(taustarhistobj$relg$t))
  taufig <- ggplot(taustarhistobj$taus, aes(x = tnow, y = tstar, ymin = tnow, ymax = tstar, colour = tnowf)) +
    geom_linerange(size=1.5) + xlab("COMM") + #xlab(expression(tau['*']^(t[now]))) + 
    theme(axis.title.x = element_blank()) + coord_flip(xlim = taufigxlims, ylim = taufigylims) +
    geom_text(aes(label = round(tstar, taulabelround), hjust = taulabelhjust)) +
    scale_x_continuous(breaks = taustarhistobj$taus$tnow, minor_breaks = NULL) + scale_colour_discrete(drop = FALSE) +
    guides(colour=guide_legend(title=expression(t[now])))
  grid.arrange(relfig, gfig, taufig, nrow = 3, heights = c(1, 1, 0.5))
}

# functions to update prior to posterior parameters for Weibull likelihood with censored observations
# set fts = NULL when no failures occurred
# can be used to include test data and to update after an operational cycle ends
# fts   vector of observed failure times
# cts   vector of censoring times
# beta  Weibull shape parameter
nn <- function(n0, fts)
  n0 + length(fts)
yn <- function(n0, y0, fts, cts, beta)
  (n0*y0 + sum(fts^beta) + sum(cts^beta)) / (n0 + length(fts))
nnyn <- function(n0y0, fts, cts, beta)
  c(nn(n0y0[1], fts), yn(n0y0[1], n0y0[2], fts, cts, beta))
# update n0y0 list: 
nnynlist <- function(n0y0list, ftslist, ctslist, beta){
  K <- length(n0y0list)
  res <- as.list(1:K)
  for (k in 1:K){
    res[[k]] <- nnyn(n0y0list[[k]], ftslist[[k]], ctslist[[k]], beta[k])
  }
  res
}


# produces survival signature matrix for one component of type "name",
# for use in nonParBayesSystemInference()
oneCompSurvSign <- function(name){
  res <- data.frame(name=c(0,1), Probability=c(0,1))
  names(res)[1] <- name
  res
}



# function to print component prior predictive reliability functions for visualizing prior assumptions
# comphor  time horizon for the plot
# complen  number of points to evaluate for component reliability plots
brcomprelplotter <- function(comphor = 10, complen = 201){
  brfts0 <- list(C=NULL, H=NULL, M=NULL, P=NULL)
  onecomp <- graph.formula(s -- 1 -- t)
  Ccomp <- setCompTypes(onecomp, list(C = "1"))
  Hcomp <- setCompTypes(onecomp, list(H = "1"))
  Mcomp <- setCompTypes(onecomp, list(M = "1"))
  Pcomp <- setCompTypes(onecomp, list(P = "1"))
  Crel1 <- sysrelnowhor(computeSystemSurvivalSignature(Ccomp), br1n0y0[1], br1beta[1], brfts0[1],
                        tnow = 0, hor = comphor, seqlen = complen)
  Hrel1 <- sysrelnowhor(computeSystemSurvivalSignature(Hcomp), br1n0y0[2], br1beta[2], brfts0[2],
                        tnow = 0, hor = comphor, seqlen = complen)
  Mrel1 <- sysrelnowhor(computeSystemSurvivalSignature(Mcomp), br1n0y0[3], br1beta[3], brfts0[3],
                        tnow = 0, hor = comphor, seqlen = complen)
  Prel1 <- sysrelnowhor(computeSystemSurvivalSignature(Pcomp), br1n0y0[4], br1beta[4], brfts0[4],
                        tnow = 0, hor = comphor, seqlen = complen)
  compdf1 <- rbind(data.frame(Crel1, comp = "C"), data.frame(Hrel1, comp = "H"),
                   data.frame(Mrel1, comp = "M"), data.frame(Prel1, comp = "P"))
  compprior1fig1 <- ggplot(compdf1, aes(x = t, y = rel)) + geom_line(aes(group = comp)) + 
    ylab("Reliability") + scale_x_continuous(breaks=seq(0, 10, by=2), minor_breaks=0:10) +
    facet_wrap(~ comp, nrow = 2, scales = "free_y") + coord_cartesian(ylim = c(0, 1))
  print(compprior1fig1)
}



#