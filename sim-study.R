# --------------------------------------------------------------------------- #
# CBM systems paper - code for example figs - simulation study
# --------------------------------------------------------------------------- #

source("cbm-sys-other.R")

simsummary <- function(simlist){
  data.frame(id = 1:length(simlist),
             nfails = sapply(simlist, function(res) sum(res$tfunc != res$tend)),
             meantend = sapply(simlist, function(res) mean(res$tend)),
             meancostrate = sapply(simlist, function(res) weighted.mean(res$costrate, res$tend)))
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

# just one 5-cycle for the start and to illustrate
br1simN51 <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data1[[1]], n0y0 = br1n0y0, beta = br1beta,
                       tnowstep = 0.1, hor = 4, thresh = 0.5, seqlen = 401)
br1simN51$tend; br1simN51$tfunc; br1simN51$costrate
br1simN51fig1 <- ggplot(br1simN51$res, aes(x = tnow, y = taustar)) + 
  geom_line(aes(colour = cycle, group = cycle)) + geom_point(aes(colour = cycle, group = cycle)) +
  xlab(expression(t[now])) + ylab(expression(paste(tau["*"]^(t[now]), (t[now])))) +
  guides(colour = guide_legend(title="Cycle")) + scale_y_continuous(breaks=seq(0, 2, by=0.5), minor_breaks=seq(0, 2, by=0.25))
pdf("br1simN51fig1.pdf", width = 6, height = 4)
br1simN51fig1
dev.off()

br1sim1 <- list()     # our model
br1sim1pr <- list()   # do not update params during cycle, but at end of cycle 
br1sim1prpr <- list() # never update params

for (i in 1:20){
  cat("Repetition", i, ": full update\n")
  br1sim1[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data1[[i]], n0y0 = br1n0y0,
                            beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.5, seqlen = 401)
  cat("Repetition", i, ": end of cycle update only\n")
  br1sim1pr[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data1[[i]], n0y0 = br1n0y0,
                              beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.5, seqlen = 401, prior = TRUE)
  cat("Repetition", i, ": no update\n")
  br1sim1prpr[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data1[[i]], n0y0 = br1n0y0,
                                beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.5, seqlen = 401, prior = TRUE, cycleupdate = FALSE)
}

br1sim1summary <- simsummary(br1sim1)
br1sim1prsummary <- simsummary(br1sim1pr)
br1sim1prprsummary <- simsummary(br1sim1prpr)

br1sim1fig1 <- ggplot(melt(br1sim1summary, "id"), aes(x = id, y = value)) + 
  geom_line(aes(colour = variable, group = variable)) + geom_point(aes(colour = variable, group = variable))
br1sim1fig1

br1sim1fig2 <- ggplot(melt(br1sim1prsummary, "id"), aes(x = id, y = value)) + 
  geom_line(aes(group = variable)) + geom_point(aes(group = variable)) + facet_wrap(~ variable, nrow = 1, scales = "free_y")
br1sim1fig2

br1sim1fig3 <- ggplot(melt(br1sim1prprsummary, "id"), aes(x = id, y = value)) + 
  geom_line(aes(colour = variable, group = variable)) + geom_point(aes(colour = variable, group = variable))
br1sim1fig3

br1sim1summaryall <- rbind(data.frame(sim = "Continuous update", br1sim1summary),
                           data.frame(sim = "Cycle end update", br1sim1prsummary),
                           data.frame(sim = "No update", br1sim1prprsummary))

br1sim1fig4 <- ggplot(melt(br1sim1summaryall, c("id", "sim")), aes(x = id, y = value)) +
  geom_line(aes(linetype = sim, group = sim)) + geom_point(aes(shape = sim, group = sim)) +
  facet_wrap(~ variable, nrow = 1, scales = "free_y", labeller = simlabeller)
br1sim1fig4

br1sim1fig5 <- ggplot(melt(br1sim1summaryall, c("id", "sim")), aes(x = id, y = value)) +
  geom_line(aes(group = sim)) + geom_point(aes(group = sim)) +
  facet_grid(variable ~ sim, scales = "free_y", labeller = simlabeller2) +
  xlab("5-cycle repetition number") + theme(axis.title.y = element_blank())
pdf("br1sim1fig5.pdf", width = 6, height = 6)
br1sim1fig5
dev.off()

# sim 2: failures earlier than expected
# -------------------------------------

# simulate 20 5-cycle machines
set.seed(2211)
br1sim5cycle20data2 <- list()
for (i in 1:20){
  br1sim5cycle20data2[[i]] <- brWeibullData(5, br1beta, 0.5*br1mttf)
}

br1sim2 <- list()     # our model
br1sim2pr <- list()   # do not update params during cycle, but at end of cycle 
br1sim2prpr <- list() # never update params

for (i in 1:20){
  cat("Repetition", i, ": full update\n")
  br1sim2[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data2[[i]], n0y0 = br1n0y0,
                            beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.5, seqlen = 401)
  cat("Repetition", i, ": end of cycle update only\n")
  br1sim2pr[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data2[[i]], n0y0 = br1n0y0,
                              beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.5, seqlen = 401, prior = TRUE)
  cat("Repetition", i, ": no update\n")
  br1sim2prpr[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data2[[i]], n0y0 = br1n0y0,
                                beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.5, seqlen = 401, prior = TRUE, cycleupdate = FALSE)
}

br1sim2summary <- simsummary(br1sim2)
br1sim2prsummary <- simsummary(br1sim2pr)
br1sim2prprsummary <- simsummary(br1sim2prpr)

br1sim2fig1 <- ggplot(melt(br1sim2summary, "id"), aes(x = id, y = value)) + 
  geom_line(aes(colour = variable, group = variable)) + geom_point(aes(colour = variable, group = variable))
br1sim2fig1

br1sim2fig2 <- ggplot(melt(br1sim2prsummary, "id"), aes(x = id, y = value)) + 
  geom_line(aes(group = variable)) + geom_point(aes(group = variable)) + facet_wrap(~ variable, nrow = 1, scales = "free_y")
br1sim2fig2

br1sim2fig3 <- ggplot(melt(br1sim2prprsummary, "id"), aes(x = id, y = value)) + 
  geom_line(aes(colour = variable, group = variable)) + geom_point(aes(colour = variable, group = variable))
br1sim2fig3

br1sim2summaryall <- rbind(data.frame(sim = "Continuous update", br1sim2summary),
                           data.frame(sim = "Cycle end update", br1sim2prsummary),
                           data.frame(sim = "No update", br1sim2prprsummary))

br1sim2fig4 <- ggplot(melt(br1sim2summaryall, c("id", "sim")), aes(x = id, y = value)) +
  geom_line(aes(linetype = sim, group = sim)) + geom_point(aes(shape = sim, group = sim)) +
  facet_wrap(~ variable, nrow = 1, scales = "free_y", labeller = simlabeller) +
  bottomlegend
pdf("br1sim2fig4.pdf", width = 7, height = 3)
br1sim2fig4
dev.off()

br1sim2fig5 <- ggplot(melt(br1sim2summaryall, c("id", "sim")), aes(x = id, y = value)) +
  geom_line(aes(group = sim)) + geom_point(aes(group = sim)) + 
  facet_grid(variable ~ sim, scales = "free_y", labeller = simlabeller2) +
  xlab("5-cycle repetition number") + theme(axis.title.y = element_blank()) +
  geom_hline(yintercept = median(value))
pdf("br1sim2fig5.pdf", width = 6, height = 6)
br1sim2fig5
dev.off()

# sim 3: failures later than expected
# -----------------------------------

# simulate 20 5-cycle machines
set.seed(2411)
br1sim5cycle20data3 <- list()
for (i in 1:20){
  br1sim5cycle20data3[[i]] <- brWeibullData(5, br1beta, 2*br1mttf)
}

br1sim3 <- list()     # our model
br1sim3pr <- list()   # do not update params during cycle, but at end of cycle 
br1sim3prpr <- list() # never update params

for (i in 1:20){
  cat("Repetition", i, ": full update\n")
  br1sim3[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data3[[i]], n0y0 = br1n0y0,
                            beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.5, seqlen = 401)
  cat("Repetition", i, ": end of cycle update only\n")
  br1sim3pr[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data3[[i]], n0y0 = br1n0y0,
                              beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.5, seqlen = 401, prior = TRUE)
  cat("Repetition", i, ": no update\n")
  br1sim3prpr[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data3[[i]], n0y0 = br1n0y0,
                                beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.5, seqlen = 401, prior = TRUE, cycleupdate = FALSE)
}

br1sim3summary <- simsummary(br1sim3)
br1sim3prsummary <- simsummary(br1sim3pr)
br1sim3prprsummary <- simsummary(br1sim3prpr)

br1sim3fig5 <- ggplot(melt(br1sim3summaryall, c("id", "sim")), aes(x = id, y = value)) +
  geom_line(aes(group = sim)) + geom_point(aes(group = sim)) + 
  facet_grid(variable ~ sim, scales = "free_y", labeller = simlabeller2) +
  xlab("5-cycle repetition number") + theme(axis.title.y = element_blank()) +
  geom_hline(yintercept = median(value))
pdf("br1sim3fig5.pdf", width = 6, height = 6)
br1sim3fig5
dev.off()

#