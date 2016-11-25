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
             failed = sapply(simlist, function(res) res$tfunc != res$tend),
             tend = sapply(simlist, function(res) res$tend),
             costrate = sapply(simlist, function(res) res$costrate))
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

br1sim2fig1T <- ggplot(melt(br1sim2Tsummary, "id"), aes(x = id, y = value)) + 
  geom_line(aes(colour = variable, group = variable)) + geom_point(aes(colour = variable, group = variable))
br1sim2fig1T

br1sim2fig2T <- ggplot(melt(br1sim2prTsummary, "id"), aes(x = id, y = value)) + 
  geom_line(aes(group = variable)) + geom_point(aes(group = variable)) + facet_wrap(~ variable, nrow = 1, scales = "free_y")
br1sim2fig2T

br1sim2fig3T <- ggplot(melt(br1sim2prprTsummary, "id"), aes(x = id, y = value)) + 
  geom_line(aes(colour = variable, group = variable)) + geom_point(aes(colour = variable, group = variable))
br1sim2fig3T

br1sim2prTsummary2 <- simsummary2(br1sim2prT)
br1sim2fig2Ta <- ggplot(melt(br1sim2prTsummary2, c("repid", "cycid")), aes(x = repid, y = value)) + 
  geom_point(aes(group = variable, colour = cycid)) +
  facet_wrap(~ variable, nrow = 1, scales = "free_y")
br1sim2fig2Ta



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


#