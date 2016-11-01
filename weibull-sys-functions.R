###############################################################################
#       Functions to calculate and plot system reliability functions          #
#       for systems consisting of Weibull components with fixed shape         #
###############################################################################

library(actuar) # provides dinvgamma(x, shape, scale) where shape = alpha and scale = beta
library(ReliabilityTheory)
library(ggplot2)
library(reshape2)

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
  ek <- unlist(lapply(fts, length))
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
#          time of this component; can be NULL if component does not fail
# n0y0     list of K prior parameter pairs c(n0,y0)
# beta     vector of K fixed weibull shape parameters
# tnowvec  vector of times for which to evaluate sysrelnow
# hor      how many time units from tnow into the future to calculate sysrelnow
# seqlen   at how many points to evaluate sysrelnow
sysrelnowhist <- function(sys, ctypes, compfts, n0y0, beta, tnowvec, hor, seqlen = 101, prior = FALSE){
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
  res <- data.frame()
  for (tnow in tnowvec){
    # index in siglist, ftslist, presentctypes for current (tnow) status
    tnowi <- max(which(ftschron <= tnow))
    prctps <- presentctypes[[tnowi]]
    cat("calculating sysrel for tnow =", tnow, "...\n")
    current <- data.frame(tnow = tnow,
                          sysrelnowhor(survsign = siglist[[tnowi]], n0y0 = n0y0[prctps],
                                       beta = beta[prctps], fts = ftslist[[tnowi]][prctps],
                                       tnow = tnow, hor = hor, seqlen = seqlen, prior = prior))
    res <- rbind(res, current)
  }
  return(res)
}

# now gtau and tau*
  

# calculates the system reliability function for the 2*2^K corner priors
# (all combinations of K lower and upper n's, for both lower and upper y)
# luckobjlist list of K of prior luckmodel objetcs
# survsign    data frame with the survival signature as output by computeSystemSurvivalSignature()
# beta        vector of K fixed weibull shape parameters
# fts         list of K vectors giving the observed failure times for the compents of type 1,...,K;
#             the list element should be NULL if no failure has been observed for type k, 
# tnow        time until the system is observed
# tvec        vector of timepoints for which to calculate P(T_sys > t), t > t_now
# nk          vector of length K giving the number of components of each type,
#             this is needed only when reduced survival signature is used
# prior       whether the prior system relibability should be calculated
fourKcornersSysrel <- function(luckobjlist, survsign, beta, fts, tnow, tvec, nk = NULL, prior = FALSE){
  K <- length(luckobjlist)
  ncomb <- 2^K
  combs <- as.list(rep(NA,K))
  for (k in 1:K)
    combs[[k]] <- c(1,2)
  # with all y0[1], jede spalte eine kombination, darunter dann tvec 
  lowercorners <- t(expand.grid(combs)) # transposed!
  NAmat <- matrix(NA, nrow = length(tvec), ncol = ncomb)
  lowercorners <- rbind(lowercorners, NAmat)
  uppercorners <- lowercorners
  # all n0 combinations with lower y's
  for (i in 1:ncomb){ # over columns!
    comblist <- as.list(rep(NA, K))
    for (k in 1:K)    # over rows!
      comblist[[k]] <- c(n0(luckobjlist[[k]])[lowercorners[k,i]], y0(luckobjlist[[k]])[1])
    rvec <- sapply(tvec, FUN = sysrel, n0y0 = comblist,
                   survsign = survsign, beta=beta, fts = fts, tnow = tnow, nk = nk, prior = prior)
    lowercorners[(1:length(rvec))+K,i] <- rvec
  }
  # all n0 combinations with upper y's
  for (i in 1:ncomb){ # over columns!
    comblist <- as.list(rep(NA, K))
    for (k in 1:K)    # over rows!
      comblist[[k]] <- c(n0(luckobjlist[[k]])[lowercorners[k,i]], y0(luckobjlist[[k]])[2])
    rvec <- sapply(tvec, FUN = sysrel, n0y0 = comblist,
                   survsign = survsign, beta=beta, fts = fts, tnow = tnow, nk = nk, prior = prior)
    uppercorners[(1:length(rvec))+K,i] <- rvec
  }
  list(lower = lowercorners, upper = uppercorners)
}

# plots one of the result data frames produced by fourKcornersSysrel()
# tvec    vector of t values, used for fourKcornersSysrel()
# rframe  one of the two data frames produced by fourKcornersSysrel()
# legend  if to print a legend
# add     just adds lines if TRUE
fourKcornersSysrelPlot <- function(tvec, rframe, legend = FALSE, add = FALSE,
                                   col = NULL, lty = NULL, ylim = c(0,1), ...){
  plotmat <- rframe[dimnames(rframe)[[1]] == "",]
  combmat <- rframe[dimnames(rframe)[[1]] != "",]
  K <- dim(combmat)[1]
  ncomb <- dim(combmat)[2]
  if(is.null(col)) col <- 1
  if(is.null(lty)) lty <- 1
  if(!add)
    matplot(tvec, plotmat, type = "l", col = col, lty = lty, ylim = ylim, ...)
  else
    matlines(tvec, plotmat, col = col, lty = lty, ...)
  if(legend){
    #if(is.null(col) & is.null(lty))
    #  stop("supply at least one of col and lty to get a legend")
    #if(is.null(col)) col <- 1
    #if(is.null(lty)) lty <- 1      
    legend("topright", legend=apply(combmat, 2, function(x) paste(x,collapse="")),
           lty = lty, col = col, lwd = 2)
  }
}

# calculates the lower and upper system reliability for fixed t,
# together with parameter pairs which achieve these bounds,
# using optim over the K n0-parameters
# this implements (21) (last equation in Sec 4.4)
#
# luckobjlist list of K of prior luckmodel objetcs
# survsign    data frame with the survival signature as output by computeSystemSurvivalSignature()
# beta        vector of K fixed weibull shape parameters
# fts         list of K vectors giving the observed failure times for the compents of type 1,...,K;
#             the list element should be NULL if no failure has been observed for type k
# tnow        time until the system is observed
# t           timepoint for which to calculate P(T_sys > t), t > t_now
# returnasvec if output should be a vector (lower bound, upper bound, l.b. parameters, u.b. parameters)
# nk          vector of length K giving the number of components of each type,
#             this is needed only when reduced survival signature is used
# prior       whether the prior system relibability should be calculated
sysrelLuck <- function(luckobjlist, survsign, beta, fts, tnow, t, returnasvec = FALSE, nk = NULL, prior = FALSE){
  K <- length(luckobjlist)
  optimfu <- function(n0vec, y0vec, K, survsign, beta, fts, tnow, t, nk, prior = prior){
    n0y0 <- as.list(rep(NA, K))
    for (k in 1:K)
      n0y0[[k]] <- c(n0vec[k], y0vec[k])
    #cat("optimfu called with n0y0vec:\n")
    #print(n0y0)
    sysrel(n0y0 = n0y0, survsign = survsign, beta = beta, fts = fts, tnow = tnow, t = t, nk = nk, prior = prior)
  }
  parl <- numeric(K)
  paru <- numeric(K)
  y0l <- y0u <- numeric(K)
  for (k in 1:K){
    parl[k] <- n0(luckobjlist[[k]])[1]
    paru[k] <- n0(luckobjlist[[k]])[2]
    y0l[k] <- y0(luckobjlist[[k]])[1]
    y0u[k] <- y0(luckobjlist[[k]])[2]
  }
  optl <- optim(par = parl, fn = optimfu, method = "L-BFGS-B", lower = parl, upper = paru,
                y0vec = y0l, K = K, survsign = survsign, beta = beta, fts = fts, tnow = tnow, t = t, nk = nk, prior = prior)
  optu <- optim(par = paru, fn = optimfu, method = "L-BFGS-B", lower = parl, upper = paru,
                y0vec = y0u, K = K, survsign = survsign, beta = beta, fts = fts, tnow = tnow, t = t, nk = nk, prior = prior,
                control = list(fnscale = -1))
  if(returnasvec)
    return(c(optl$value, optu$value, optl$par, optu$par))
  else
    return(list(rel = c(optl$value, optu$value), lpar = optl$par, upar = optu$par))  
}


# calculates the the lower and upper system reliability function for vector of timepoints tvec
# luckobjlist list of K of prior luckmodel objetcs
# survsign    data frame with the survival signature as output by computeSystemSurvivalSignature()
# beta        vector of K fixed weibull shape parameters
# fts         list of K vectors giving the observed failure times for the compents of type 1,...,K;
#             the list element should be NULL if no failure has been observed for type k, 
# tnow        time until the system is observed
# tvec        vector of timepoints for which to calculate P(T_sys > t), t > t_now
# nk          vector of length K giving the number of components of each type,
#             this is needed only when reduced survival signature is used
# prior       whether the prior system relibability should be calculated
sysrelPbox <- function(luckobjlist, survsign, beta, fts, tnow, tvec, nk = NULL, prior = FALSE, ...){
  res <- sapply(tvec, FUN = sysrelLuck, luckobjlist = luckobjlist, survsign = survsign,
                beta = beta, fts = fts, tnow = tnow, nk = nk, prior = prior, returnasvec = TRUE)
  res <- cbind(tvec, t(res))
  K <- length(luckobjlist)
  colnames(res) <- c("tvec", "lower", "upper", paste("ln", 1:K, sep=""), paste("un", 1:K, sep=""))
  return(res)
}

# plots the the lower and upper system reliability function for vector of timepoints tvec
# in default graphics as grey area via polygon()
#
# res  table of lower and upper system reliability values as output by sysrelPbox()
plotSysrelPbox <- function(res, add = FALSE, ylim = c(0,1), xlab = "t", ylab = expression(R[sys](t)),
                           polygonBorderCol = NA, polygonFillCol = "grey", ...){
  tvec <- res[,1]
  if(!add){
    plot(c(min(tvec), max(tvec)), c(0,0), type = "n", ylim = ylim, xlab = "", ylab = "", ...)
    mtext(xlab, 1, 2)
    mtext(ylab, 2, 2)
  }
  tvecrev <- numeric(length(tvec))
  lower <- res[,2]
  upper <- numeric(length(tvec))
  for (t in 1:length(tvec)){
    tvecrev[length(tvec)-t+1] <- tvec[t]
    upper[length(tvec)-t+1] <- res[t,3]
  }
  polygon(c(tvec, tvecrev), c(lower, upper), border = polygonBorderCol, col = polygonFillCol)
}

# -------------------------------------

# prior predictive density
prpr <- function(t, n0, y0, beta = 2)
  beta * t^(beta - 1) * (n0 + 1) * (n0 * y0)^(n0 + 1) * (n0 * y0 + t^beta)^(-n0 - 2)

# prior predictive reliability function, works only for t vectors starting at 0
prprRold <- function(t, n0, y0, beta = 2){
  dens <- pmax(prpr(t = t, n0 = n0, y0 = y0, beta = beta), 0)
  1-cumsum(dens)/(length(t)/max(t)) # divide by evaluation points per unit
}

# prior predictive reliability function with numerical integration
prprR <- function(t, n0, y0, beta = 2){
  1 - sapply(t, function(sgt) integrate(f = prpr, lower = 1e-5, upper = sgt, n0 = n0, y0 = y0, beta = beta)$value)
}

# prior predictive reliability function with numerical integration, with loop to see where the problem is
prprRloop <- function(t, n0, y0, beta = 2){
  res <- numeric(length(t))
  for (i in 1:length(t)){
    sgt <- t[i]
    cat(i, ": t = ", sgt, "\n")
    res[i] <- 1 - integrate(f = prpr, lower = 1e-5, upper = sgt, n0 = n0, y0 = y0, beta = beta)$value
  }
  res
}

# functions to update prior to posterior parameters for Weibull likelihood with censored observations
# set fts = NULL when no failures occurred
nn <- function(n0, fts)
  n0 + length(fts)
yn <- function(n0, y0, fts, beta, n = length(fts), tnow = max(fts))
  (n0*y0 + (n - length(fts))*tnow^beta + sum(fts^beta))/(n0 + length(fts))

# prior predictive reliability function, works only for t vectors starting at 0
# going over n0 grid, default is only lower and upper n0
prprRluck <- function(tvec, luckobj, beta = 2, n0gridlen = 2){
  n0grid <- seq(n0(luckobj)[1], n0(luckobj)[2], length.out = n0gridlen)
  lmat <- sapply(n0grid, function(n0) prprR(tvec, n0, y0(luckobj)[1], beta))
  umat <- sapply(n0grid, function(n0) prprR(tvec, n0, y0(luckobj)[2], beta))
  data.frame(tvec = tvec, lower = apply(lmat, 1, min), upper = apply(umat, 1, max))
}


# posterior predictive reliability function, works only for t vectors starting at 0
# going over n0 grid, default is only lower and upper n0
# n and tnow used only for update with censored observations
poprRluck <- function(tvec, luckobj, beta = 2, fts, n, tnow, n0gridlen = 2){
  n0grid <- seq(n0(luckobj)[1], n0(luckobj)[2], length.out = n0gridlen)
  lmat <- sapply(n0grid, function(n0) prprR(tvec, nn(n0, fts),
                                            yn(n0, y0(luckobj)[1], fts, beta, n, tnow), beta))
  umat <- sapply(n0grid, function(n0) prprR(tvec, nn(n0, fts),
                                            yn(n0, y0(luckobj)[2], fts, beta, n, tnow), beta))
  data.frame(tvec = tvec, lower = apply(lmat, 1, min), upper = apply(umat, 1, max))
}


# checking only the four corners of the posterior parameter set
poprRluck2 <- function(tvec, luckobj, beta = 2, fts, n, tnow){
  l1 <- prprR(tvec, nn(n0(luckobj)[1], fts), yn(n0(luckobj)[1], y0(luckobj)[1], fts, beta, n, tnow), beta)
  print("l1 ok")
  l2 <- prprR(tvec, nn(n0(luckobj)[2], fts), yn(n0(luckobj)[2], y0(luckobj)[1], fts, beta, n, tnow), beta)
  print("l2 ok")
  u1 <- prprR(tvec, nn(n0(luckobj)[1], fts), yn(n0(luckobj)[1], y0(luckobj)[2], fts, beta, n, tnow), beta)
  print("u1 ok")
  u2 <- prprR(tvec, nn(n0(luckobj)[2], fts), yn(n0(luckobj)[2], y0(luckobj)[2], fts, beta, n, tnow), beta)
  print("u2 ok")
  data.frame(tvec = tvec, lower = pmin(l1, l2), upper = pmax(u1, u2))
}


# produces survival signature matrix for one component of type "name",
# for use in nonParBayesSystemInference()
oneCompSurvSign <- function(name){
  res <- data.frame(name=c(0,1), Probability=c(0,1))
  names(res)[1] <- name
  res
}

#