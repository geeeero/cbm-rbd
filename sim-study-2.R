# --------------------------------------------------------------------------- #
# CBM systems paper - code for example figs - simulation study
# using total unit cost rate now!
# choose lower thresh now???
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

br1sim1T <- list()     # our model
br1sim1prT <- list()   # do not update params during cycle, but at end of cycle 
br1sim1prprT <- list() # never update params

for (i in 1:20){
  cat("Repetition", i, ": full update\n")
  br1sim1T[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data1[[i]], n0y0 = br1n0y0,
                             beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.5, seqlen = 401, onecycle = FALSE)
  cat("Repetition", i, ": end of cycle update only\n")
  br1sim1prT[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data1[[i]], n0y0 = br1n0y0,
                               beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.5, seqlen = 401, prior = TRUE, onecycle = FALSE)
  cat("Repetition", i, ": no update\n")
  br1sim1prprT[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data1[[i]], n0y0 = br1n0y0,
                                 beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.5, seqlen = 401, prior = TRUE, cycleupdate = FALSE, onecycle = FALSE)
}

br1sim1Tsummary <- simsummary(br1sim1T)
br1sim1prTsummary <- simsummary(br1sim1prT)
br1sim1prprTsummary <- simsummary(br1sim1prprT)

br1sim1Tsummaryall <- rbind(data.frame(sim = "Continuous update", br1sim1Tsummary),
                            data.frame(sim = "Cycle end update", br1sim1prTsummary),
                            data.frame(sim = "No update", br1sim1prprTsummary))

br1sim1fig5T <- ggplot(melt(br1sim1Tsummaryall, c("id", "sim")), aes(x = id, y = value)) +
  geom_line(aes(group = sim)) + geom_point(aes(group = sim)) + 
  facet_grid(variable ~ sim, scales = "free_y", labeller = simlabeller2) +
  xlab("5-cycle repetition number") + theme(axis.title.y = element_blank())
pdf("br1sim1fig5T.pdf", width = 6, height = 6)
br1sim1fig5T
dev.off()


# sim 2: failures earlier than expected
# -------------------------------------

# simulate 20 5-cycle machines
set.seed(2211)
br1sim5cycle20data2 <- list()
for (i in 1:20){
  br1sim5cycle20data2[[i]] <- brWeibullData(5, br1beta, 0.5*br1mttf)
}

br1sim2T <- list()     # our model
br1sim2prT <- list()   # do not update params during cycle, but at end of cycle 
br1sim2prprT <- list() # never update params

for (i in 1:20){
  cat("Repetition", i, ": full update\n")
  br1sim2T[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data2[[i]], n0y0 = br1n0y0,
                             beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.5, seqlen = 401, onecycle = FALSE)
  cat("Repetition", i, ": end of cycle update only\n")
  br1sim2prT[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data2[[i]], n0y0 = br1n0y0,
                               beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.5, seqlen = 401, prior = TRUE, onecycle = FALSE)
  cat("Repetition", i, ": no update\n")
  br1sim2prprT[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data2[[i]], n0y0 = br1n0y0,
                                 beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.5, seqlen = 401, prior = TRUE, cycleupdate = FALSE, onecycle = FALSE)
}

br1sim2Tsummary <- simsummary(br1sim2T)
br1sim2prTsummary <- simsummary(br1sim2prT)
br1sim2prprTsummary <- simsummary(br1sim2prprT)

br1sim2Tsummaryall <- rbind(data.frame(sim = "Continuous update", br1sim2Tsummary),
                            data.frame(sim = "Cycle end update", br1sim2prTsummary),
                            data.frame(sim = "No update", br1sim2prprTsummary))

br1sim2fig5T <- ggplot(melt(br1sim2Tsummaryall, c("id", "sim")), aes(x = id, y = value)) +
  geom_line(aes(group = sim)) + geom_point(aes(group = sim)) + 
  facet_grid(variable ~ sim, scales = "free_y", labeller = simlabeller2) +
  xlab("5-cycle repetition number") + theme(axis.title.y = element_blank())
pdf("br1sim2fig5T.pdf", width = 6, height = 6)
br1sim2fig5T
dev.off()

br1sim2Tsummary2 <- simsummary2(br1sim2T)
br1sim2prTsummary2 <- simsummary2(br1sim2prT)
br1sim2prprTsummary2 <- simsummary2(br1sim2prprT)

br1sim2Tsummary2all <- rbind(data.frame(sim = "Continuous update", br1sim2Tsummary2),
                             data.frame(sim = "Cycle end update", br1sim2prTsummary2),
                             data.frame(sim = "No update", br1sim2prprTsummary2))

br1sim2fig2Ta <- ggplot(melt(br1sim2Tsummary2all, c("repid", "cycid", "sim")), aes(x = repid, y = value)) + 
  stat_boxplot(aes(group = repid)) +
  #geom_point(aes(group = variable, colour = cycid)) +
  facet_grid(variable ~ sim, scales = "free_y")
br1sim2fig2Ta

# sim 3: failures later than expected
# -------------------------------------

# simulate 20 5-cycle machines
set.seed(2411)
br1sim5cycle20data3 <- list()
for (i in 1:20){
  br1sim5cycle20data3[[i]] <- brWeibullData(5, br1beta, 2*br1mttf)
}

br1sim3T <- list()     # our model
br1sim3prT <- list()   # do not update params during cycle, but at end of cycle 
br1sim3prprT <- list() # never update params

for (i in 1:20){
  cat("Repetition", i, ": full update\n")
  br1sim3T[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data3[[i]], n0y0 = br1n0y0,
                             beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.5, seqlen = 401, onecycle = FALSE)
  cat("Repetition", i, ": end of cycle update only\n")
  br1sim3prT[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data3[[i]], n0y0 = br1n0y0,
                               beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.5, seqlen = 401, prior = TRUE, onecycle = FALSE)
  cat("Repetition", i, ": no update\n")
  br1sim3prprT[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data3[[i]], n0y0 = br1n0y0,
                                 beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.5, seqlen = 401, prior = TRUE, cycleupdate = FALSE, onecycle = FALSE)
}

br1sim3Tsummary <- simsummary(br1sim3T)
br1sim3prTsummary <- simsummary(br1sim3prT)
br1sim3prprTsummary <- simsummary(br1sim3prprT)

br1sim3Tsummaryall <- rbind(data.frame(sim = "Continuous update", br1sim3Tsummary),
                            data.frame(sim = "Cycle end update", br1sim3prTsummary),
                            data.frame(sim = "No update", br1sim3prprTsummary))

br1sim3fig5T <- ggplot(melt(br1sim3Tsummaryall, c("id", "sim")), aes(x = id, y = value)) +
  geom_line(aes(group = sim)) + geom_point(aes(group = sim)) + 
  facet_grid(variable ~ sim, scales = "free_y", labeller = simlabeller2) +
  xlab("5-cycle repetition number") + theme(axis.title.y = element_blank())
pdf("br1sim3fig5T.pdf", width = 6, height = 6)
br1sim3fig5T
dev.off()


# sim 4: failures much earlier than expected
# ------------------------------------------

# simulate 20 5-cycle machines
set.seed(2611)
br1sim5cycle20data4 <- list()
for (i in 1:20){
  br1sim5cycle20data4[[i]] <- brWeibullData(5, br1beta, 0.2*br1mttf)
}

br1sim4T <- list()     # our model
br1sim4prT <- list()   # do not update params during cycle, but at end of cycle 
br1sim4prprT <- list() # never update params

for (i in 1:20){
  cat("Repetition", i, ": full update\n")
  br1sim4T[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data4[[i]], n0y0 = br1n0y0,
                             beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.5, seqlen = 401, onecycle = FALSE)
  cat("Repetition", i, ": end of cycle update only\n")
  br1sim4prT[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data4[[i]], n0y0 = br1n0y0,
                               beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.5, seqlen = 401, prior = TRUE, onecycle = FALSE)
  cat("Repetition", i, ": no update\n")
  br1sim4prprT[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data4[[i]], n0y0 = br1n0y0,
                                 beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.5, seqlen = 401, prior = TRUE, cycleupdate = FALSE, onecycle = FALSE)
}

br1sim4Tsummary <- simsummary(br1sim4T)
br1sim4prTsummary <- simsummary(br1sim4prT)
br1sim4prprTsummary <- simsummary(br1sim4prprT)

br1sim4Tsummaryall <- rbind(data.frame(sim = "Continuous update", br1sim4Tsummary),
                            data.frame(sim = "Cycle end update", br1sim4prTsummary),
                            data.frame(sim = "No update", br1sim4prprTsummary))

#br1sim4Tmedians <- median(br1sim4Tsummary$meantend)

br1sim4fig5T <- ggplot(melt(br1sim4Tsummaryall, c("id", "sim")), aes(x = id, y = value)) +
  geom_line(aes(group = sim)) + geom_point(aes(group = sim)) + 
  facet_grid(variable ~ sim, scales = "free_y", labeller = simlabeller2) +
  xlab("5-cycle repetition number") + theme(axis.title.y = element_blank())
#br1sim4fig5T + geom_hline(data = br1sim4Tmedians, aes(yintercept = median(value), group = c(variable)))
pdf("br1sim4fig5T.pdf", width = 6, height = 6)
br1sim4fig5T
dev.off()

#