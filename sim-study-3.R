# --------------------------------------------------------------------------- #
# CBM systems paper - code for example figs - simulation study
# using total unit cost rate & tnowstep = 0.1
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
                              "No update" = "No update",
                              "Age-based with update" = "Age-based with update",
                              "Age-based without update" = "Age-based without update",
                              "Corrective" = "Corrective"))
simlabeller3 <- as_labeller(c(nfails = expression(n[failures]),
                              meantend = expression(bar(t)[sys]),
                              meancostrate = expression(bar(g)),
                              "Continuous update" = "CBMcpu",
                              "Cycle end update" = "CBMepu",
                              "No update" = "CBMnpu",
                              "Age-based with update" = "ABMepu",
                              "Age-based without update" = "ABMnpu",
                              "Corrective" = "CM"))


# sim 1: data fits to priors
# --------------------------

# simulate 20 5-cycle machines with priors as given
set.seed(2211)
br1sim5cycle20data1 <- list()
for (i in 1:20){
  br1sim5cycle20data1[[i]] <- brWeibullData(5, br1beta, br1mttf)
}

br1sim1Tt01 <- list()     # our model
br1sim1prTt01 <- list()   # do not update params during cycle, but at end of cycle 
br1sim1prprTt01 <- list() # never update params

for (i in 1:20){
  cat("Repetition", i, ": full update\n")
  br1sim1Tt01[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data1[[i]], n0y0 = br1n0y0,
                                beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.1, seqlen = 401, onecycle = FALSE)
  cat("Repetition", i, ": end of cycle update only\n")
  br1sim1prTt01[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data1[[i]], n0y0 = br1n0y0,
                                  beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.1, seqlen = 401, prior = TRUE, onecycle = FALSE)
  cat("Repetition", i, ": no update\n")
  br1sim1prprTt01[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data1[[i]], n0y0 = br1n0y0,
                                    beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.1, seqlen = 401, prior = TRUE, cycleupdate = FALSE, onecycle = FALSE)
}

br1sim1Tt01summary <- simsummary(br1sim1Tt01)
br1sim1prTt01summary <- simsummary(br1sim1prTt01)
br1sim1prprTt01summary <- simsummary(br1sim1prprTt01)
br1sim1Tt01summaryall <- rbind(data.frame(sim = "Continuous update", br1sim1Tt01summary),
                               data.frame(sim = "Cycle end update", br1sim1prTt01summary),
                               data.frame(sim = "No update", br1sim1prprTt01summary))
br1sim1fig5Tt01 <- ggplot(melt(br1sim1Tt01summaryall, c("id", "sim")), aes(x = id, y = value)) +
  geom_line(aes(group = sim)) + geom_point(aes(group = sim)) + 
  facet_grid(variable ~ sim, scales = "free_y", labeller = simlabeller2) +
  xlab("5-cycle repetition number") + theme(axis.title.y = element_blank())

# corrective and age based

br1sim1a <- list()   # age-based policy with parameter update
br1sim1apr <- list() # age-based policy without parameter update
br1sim1c <- list()   # corrective policy

for (i in 1:20){
  cat("Repetition", i, ": age-based with update\n")
  br1sim1a[[i]] <- simNcycleAgebased(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data1[[i]], n0y0 = br1n0y0,
                                     beta = br1beta, tnowstep = 0.1, hor = 4, seqlen = 401)
  cat("Repetition", i, ": age-based without update\n")
  br1sim1apr[[i]] <- simNcycleAgebased(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data1[[i]], n0y0 = br1n0y0,
                                       beta = br1beta, tnowstep = 0.1, hor = 4, seqlen = 401, cycleupdate = FALSE)
  cat("Repetition", i, ": corrective\n")
  br1sim1c[[i]] <- simNcycleCorrective(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data1[[i]], tnowstep = 0.1)
}

br1sim1asummary <- simsummary(br1sim1a)
br1sim1aprsummary <- simsummary(br1sim1apr)
br1sim1csummary <- simsummary(br1sim1c)
br1sim1acsummaryall <- rbind(data.frame(sim = "Age-based with update", br1sim1asummary),
                             data.frame(sim = "Age-based without update", br1sim1aprsummary),
                             data.frame(sim = "Corrective", br1sim1csummary))

br1sim1fig5ac <- ggplot(melt(br1sim1acsummaryall, c("id", "sim")), aes(x = id, y = value)) +
  geom_line(aes(group = sim)) + geom_point(aes(group = sim)) + 
  facet_grid(variable ~ sim, scales = "free_y", labeller = simlabeller2) +
  xlab("5-cycle repetition number") + theme(axis.title.y = element_blank())

br1sim1Tt01acsummaryall <- rbind(data.frame(sim = "CBM-cpu", br1sim1Tt01summary),
                                 data.frame(sim = "CBM-epu", br1sim1prTt01summary),
                                 data.frame(sim = "CBM-npu", br1sim1prprTt01summary),
                                 data.frame(sim = "ABM-epu", br1sim1asummary),
                                 data.frame(sim = "ABM-npu", br1sim1aprsummary),
                                 data.frame(sim = "CM", br1sim1csummary))
names(br1sim1Tt01acsummaryall)[3:5] <- c(expression(e[sys]), expression(bar(r)[sys]), expression(bar(g)))

br1sim1fig5Tt01ac <- ggplot(melt(br1sim1Tt01acsummaryall, c("id", "sim")), aes(x = id, y = value)) +
  geom_line(aes(group = sim), size = 0.3) + geom_point(aes(group = sim), size = 0.5) + 
  facet_grid(variable ~ sim, scales = "free_y", labeller = label_parsed) +
  xlab("5-cycle repetition number") + theme(axis.title.y = element_blank())
pdf("br1sim1fig5Tt01ac.pdf", width = 6, height = 3)
br1sim1fig5Tt01ac
dev.off()


# sim 2: failures earlier than expected
# -------------------------------------

# simulate 20 5-cycle machines
set.seed(2211)
br1sim5cycle20data2 <- list()
for (i in 1:20){
  br1sim5cycle20data2[[i]] <- brWeibullData(5, br1beta, 0.5*br1mttf)
}

br1sim2Tt01 <- list()     # our model
br1sim2prTt01 <- list()   # do not update params during cycle, but at end of cycle 
br1sim2prprTt01 <- list() # never update params

for (i in 1:20){
  cat("Repetition", i, ": full update\n")
  br1sim2Tt01[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data2[[i]], n0y0 = br1n0y0,
                                beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.1, seqlen = 401, onecycle = FALSE)
  cat("Repetition", i, ": end of cycle update only\n")
  br1sim2prTt01[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data2[[i]], n0y0 = br1n0y0,
                                  beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.1, seqlen = 401, prior = TRUE, onecycle = FALSE)
  cat("Repetition", i, ": no update\n")
  br1sim2prprTt01[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data2[[i]], n0y0 = br1n0y0,
                                    beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.1, seqlen = 401, prior = TRUE, cycleupdate = FALSE, onecycle = FALSE)
}

br1sim2Tt01summary <- simsummary(br1sim2Tt01)
br1sim2prTt01summary <- simsummary(br1sim2prTt01)
br1sim2prprTt01summary <- simsummary(br1sim2prprTt01)
br1sim2Tt01summaryall <- rbind(data.frame(sim = "Continuous update", br1sim2Tt01summary),
                               data.frame(sim = "Cycle end update", br1sim2prTt01summary),
                               data.frame(sim = "No update", br1sim2prprTt01summary))
br1sim2fig5Tt01 <- ggplot(melt(br1sim2Tt01summaryall, c("id", "sim")), aes(x = id, y = value)) +
  geom_line(aes(group = sim)) + geom_point(aes(group = sim)) + 
  facet_grid(variable ~ sim, scales = "free_y", labeller = simlabeller2) +
  xlab("5-cycle repetition number") + theme(axis.title.y = element_blank())

# corrective and age based

br1sim2a <- list()   # age-based policy with parameter update
br1sim2apr <- list() # age-based policy without parameter update
br1sim2c <- list()   # corrective policy

for (i in 1:20){
  cat("Repetition", i, ": age-based with update\n")
  br1sim2a[[i]] <- simNcycleAgebased(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data2[[i]], n0y0 = br1n0y0,
                                     beta = br1beta, tnowstep = 0.1, hor = 4, seqlen = 401)
  cat("Repetition", i, ": age-based without update\n")
  br1sim2apr[[i]] <- simNcycleAgebased(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data2[[i]], n0y0 = br1n0y0,
                                       beta = br1beta, tnowstep = 0.1, hor = 4, seqlen = 401, cycleupdate = FALSE)
  cat("Repetition", i, ": corrective\n")
  br1sim2c[[i]] <- simNcycleCorrective(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data2[[i]], tnowstep = 0.1)
}

br1sim2asummary <- simsummary(br1sim2a)
br1sim2aprsummary <- simsummary(br1sim2apr)
br1sim2csummary <- simsummary(br1sim2c)
br1sim2Tt01acsummaryall <- rbind(data.frame(sim = "CBM-cpu", br1sim2Tt01summary),
                                 data.frame(sim = "CBM-epu", br1sim2prTt01summary),
                                 data.frame(sim = "CBM-npu", br1sim2prprTt01summary),
                                 data.frame(sim = "ABM-epu", br1sim2asummary),
                                 data.frame(sim = "ABM-npu", br1sim2aprsummary),
                                 data.frame(sim = "CM", br1sim2csummary))
names(br1sim2Tt01acsummaryall)[3:5] <- c(expression(e[sys]), expression(bar(r)[sys]), expression(bar(g)))

br1sim2fig5Tt01ac <- ggplot(melt(br1sim2Tt01acsummaryall, c("id", "sim")), aes(x = id, y = value)) +
  geom_line(aes(group = sim), size = 0.3) + geom_point(aes(group = sim), size = 0.5) + 
  facet_grid(variable ~ sim, scales = "free_y", labeller = label_parsed) +
  xlab("5-cycle repetition number") + theme(axis.title.y = element_blank())
pdf("br1sim2fig5Tt01ac.pdf", width = 6, height = 3)
br1sim2fig5Tt01ac
dev.off()


# sim 3: failures later than expected
# -------------------------------------

# simulate 20 5-cycle machines
set.seed(2411)
br1sim5cycle20data3 <- list()
for (i in 1:20){
  br1sim5cycle20data3[[i]] <- brWeibullData(5, br1beta, 2*br1mttf)
}

br1sim3Tt01 <- list()     # our model
br1sim3prTt01 <- list()   # do not update params during cycle, but at end of cycle 
br1sim3prprTt01 <- list() # never update params

for (i in 1:20){
  cat("Repetition", i, ": full update\n")
  br1sim3Tt01[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data3[[i]], n0y0 = br1n0y0,
                                beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.1, seqlen = 401, onecycle = FALSE)
  cat("Repetition", i, ": end of cycle update only\n")
  br1sim3prTt01[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data3[[i]], n0y0 = br1n0y0,
                                  beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.1, seqlen = 401, prior = TRUE, onecycle = FALSE)
  cat("Repetition", i, ": no update\n")
  br1sim3prprTt01[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data3[[i]], n0y0 = br1n0y0,
                                    beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.1, seqlen = 401, prior = TRUE, cycleupdate = FALSE, onecycle = FALSE)
}

br1sim3Tt01summary <- simsummary(br1sim3Tt01)
br1sim3prTt01summary <- simsummary(br1sim3prTt01)
br1sim3prprTt01summary <- simsummary(br1sim3prprTt01)
br1sim3Tt01summaryall <- rbind(data.frame(sim = "Continuous update", br1sim3Tt01summary),
                               data.frame(sim = "Cycle end update", br1sim3prTt01summary),
                               data.frame(sim = "No update", br1sim3prprTt01summary))
br1sim3fig5Tt01 <- ggplot(melt(br1sim3Tt01summaryall, c("id", "sim")), aes(x = id, y = value)) +
  geom_line(aes(group = sim)) + geom_point(aes(group = sim)) + 
  facet_grid(variable ~ sim, scales = "free_y", labeller = simlabeller2) +
  xlab("5-cycle repetition number") + theme(axis.title.y = element_blank())

# corrective and age based

br1sim3a <- list()   # age-based policy with parameter update
br1sim3apr <- list() # age-based policy without parameter update
br1sim3c <- list()   # corrective policy

for (i in 1:20){
  cat("Repetition", i, ": age-based with update\n")
  br1sim3a[[i]] <- simNcycleAgebased(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data3[[i]], n0y0 = br1n0y0,
                                     beta = br1beta, tnowstep = 0.1, hor = 4, seqlen = 401)
  cat("Repetition", i, ": age-based without update\n")
  br1sim3apr[[i]] <- simNcycleAgebased(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data3[[i]], n0y0 = br1n0y0,
                                       beta = br1beta, tnowstep = 0.1, hor = 4, seqlen = 401, cycleupdate = FALSE)
  cat("Repetition", i, ": corrective\n")
  br1sim3c[[i]] <- simNcycleCorrective(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data3[[i]], tnowstep = 0.1)
}

br1sim3asummary <- simsummary(br1sim3a)
br1sim3aprsummary <- simsummary(br1sim3apr)
br1sim3csummary <- simsummary(br1sim3c)
br1sim3Tt01acsummaryall <- rbind(data.frame(sim = "CBM-cpu", br1sim3Tt01summary),
                                 data.frame(sim = "CBM-epu", br1sim3prTt01summary),
                                 data.frame(sim = "CBM-npu", br1sim3prprTt01summary),
                                 data.frame(sim = "ABM-epu", br1sim3asummary),
                                 data.frame(sim = "ABM-npu", br1sim3aprsummary),
                                 data.frame(sim = "CM", br1sim3csummary))
names(br1sim3Tt01acsummaryall)[3:5] <- c(expression(e[sys]), expression(bar(r)[sys]), expression(bar(g)))

br1sim3fig5Tt01ac <- ggplot(melt(br1sim3Tt01acsummaryall, c("id", "sim")), aes(x = id, y = value)) +
  geom_line(aes(group = sim), size = 0.3) + geom_point(aes(group = sim), size = 0.5) + 
  facet_grid(variable ~ sim, scales = "free_y", labeller = label_parsed) +
  xlab("5-cycle repetition number") + theme(axis.title.y = element_blank())
pdf("br1sim3fig5Tt01ac.pdf", width = 6, height = 3)
br1sim3fig5Tt01ac
dev.off()

mean(br1sim3Tt01summary$meancostrate); median(br1sim3Tt01summary$meancostrate)
mean(br1sim3prTt01summary$meancostrate); median(br1sim3prTt01summary$meancostrate)
mean(br1sim3prprTt01summary$meancostrate); median(br1sim3prprTt01summary$meancostrate)
mean(br1sim3asummary$meancostrate); median(br1sim3asummary$meancostrate)
mean(br1sim3aprsummary$meancostrate); median(br1sim3aprsummary$meancostrate)
mean(br1sim3csummary$meancostrate); median(br1sim3csummary$meancostrate)


# save stuff

ststr <- paste("br1sim", 1:3, sep="")
enstr <- c("Tt01", "prTt01", "prprTt01", "a", "apr", "c")
savelist <- as.vector(sapply(ststr, function(brstr) paste(brstr, enstr, sep="")))
save(list = savelist, file = "papersimobjects")
load("papersimobjects")

#