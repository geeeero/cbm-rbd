# --------------------------------------------------------------------------- #
# CBM systems paper - code for example figs - simulation study
# age-based and corrective maintenance
# --------------------------------------------------------------------------- #

source("cbm-sys-other.R")

simsummary <- function(simlist){
  data.frame(id = 1:length(simlist),
             nfails = sapply(simlist, function(res) sum(res$tfunc != res$tend)),
             meantend = sapply(simlist, function(res) mean(res$tend)),
             meancostrate = sapply(simlist, function(res) weighted.mean(res$costrate, res$tend)))
}

simsummary2 <- function(simlist){
  nrep <- length(simlist)
  ncyc <- length(simlist[[1]]$tend)
  data.frame(repid = rep(1:nrep, each = ncyc),
             cycid = rep(1:ncyc, nrep),
             failed = as.vector(sapply(simlist, function(res) res$tfunc != res$tend)),
             tend = as.vector(sapply(simlist, function(res) res$tend)),
             costrate = as.vector(sapply(simlist, function(res) res$costrate)))
}

simlabeller <- as_labeller(c(nfails = "Number of system failures",
                             meantend = "Average cycle length",
                             meancostrate = "Average unit cost rate"))
simlabeller2 <- as_labeller(c(nfails = "Number of system failures",
                              meantend = "Average cycle length",
                              meancostrate = "Average unit cost rate",
                              "Continuous update" = "Continuous update",
                              "Cycle end update" = "Cycle end update",
                              "No update" = "No update"))

# sim 1: data fits to priors
# --------------------------

# simulate 20 5-cycle machines with priors as given
set.seed(2211)
br1sim5cycle20data1 <- list()
for (i in 1:20){
  br1sim5cycle20data1[[i]] <- brWeibullData(5, br1beta, br1mttf)
}

br1sim1C <- list()   # corrective
#br1sim1A <- list()   # age-based, update params at end of each cycle 
#br1sim1Apr <- list() # age-based, don't update params

for (i in 1:20){
  cat("Repetition", i, ": corrective\n")
  br1sim1C[[i]] <- simNcycleCorrective(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data1[[i]], tnowstep = 0.1)
}

br1sim1Csummary <- simsummary(br1sim1C)
#br1sim1Asummary <- simsummary(br1sim1A)
#br1sim1Aprsummary <- simsummary(br1sim1Apr)

br1sim1fig2C <- ggplot(melt(br1sim1Csummary, "id"), aes(x = id, y = value)) + 
  geom_line(aes(group = variable)) + geom_point(aes(group = variable)) + facet_wrap(~ variable, nrow = 1, scales = "free_y")
br1sim1fig2C


# sim 2: failures earlier than expected
# -------------------------------------

# simulate 20 5-cycle machines
set.seed(2211)
br1sim5cycle20data2 <- list()
for (i in 1:20){
  br1sim5cycle20data2[[i]] <- brWeibullData(5, br1beta, 0.5*br1mttf)
}

br1sim2C <- list()   # corrective
#br1sim2A <- list()   # age-based, update params at end of each cycle 
#br1sim2Apr <- list() # age-based, don't update params

for (i in 1:20){
  cat("Repetition", i, ": corrective\n")
  br1sim2C[[i]] <- simNcycleCorrective(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data2[[i]], tnowstep = 0.1)
}

br1sim2Csummary <- simsummary(br1sim2C)
#br1sim3Asummary <- simsummary(br1sim3A)
#br1sim3Aprsummary <- simsummary(br1sim3Apr)


# sim 3: failures later than expected
# -------------------------------------

# simulate 20 5-cycle machines
set.seed(2411)
br1sim5cycle20data3 <- list()
for (i in 1:20){
  br1sim5cycle20data3[[i]] <- brWeibullData(5, br1beta, 2*br1mttf)
}

br1sim3C <- list()   # corrective
#br1sim3A <- list()   # age-based, update params at end of each cycle 
#br1sim3Apr <- list() # age-based, don't update params

for (i in 1:20){
  cat("Repetition", i, ": corrective\n")
  br1sim3C[[i]] <- simNcycleCorrective(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data3[[i]], tnowstep = 0.1)
}

br1sim3Csummary <- simsummary(br1sim3C)
#br1sim3Asummary <- simsummary(br1sim3A)
#br1sim3Aprsummary <- simsummary(br1sim3Apr)

#