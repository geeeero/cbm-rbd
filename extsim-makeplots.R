# --------------------------------------------------------------------------- #
# CBM systems paper - extended simulation study for revision - make plots
# --------------------------------------------------------------------------- #

# load packages
library(actuar)
library(ReliabilityTheory)
library(reshape2)
library(ggplot2)
library(gridExtra)

# code files (these contain also the plotting functions)
source("weibull-sys-functions.R") 
source("cbm-sim2.R")

# plot definitions and summary functions
source("plotdefs.R")
source("extsim-summaryfunctions.R")

# Case A: failure times as expected

## load simulation objects
# load()

## apply summary functions
br1sim1Tt01summary <- simsummary(br1sim1Tt01)
br1sim1prTt01summary <- simsummary(br1sim1prTt01)
br1sim1prprTt01summary <- simsummary(br1sim1prprTt01)

br1sim1asummary <- simsummary(br1sim1a)
br1sim1aprsummary <- simsummary(br1sim1apr)
br1sim1csummary <- simsummary(br1sim1c)

br1sim1Tt01acsummaryall <- rbind(data.frame(sim = "CBM-cpu", br1sim1Tt01summary),
                                 data.frame(sim = "CBM-epu", br1sim1prTt01summary),
                                 data.frame(sim = "CBM-npu", br1sim1prprTt01summary),
                                 data.frame(sim = "ABM-epu", br1sim1asummary),
                                 data.frame(sim = "ABM-npu", br1sim1aprsummary),
                                 data.frame(sim = "CM", br1sim1csummary))
names(br1sim1Tt01acsummaryall)[3:5] <- c(expression(e[sys]), expression(bar(r)[sys]), expression(bar(g)))

## the plot
br1sim1fig5Tt01ac <- ggplot(melt(br1sim1Tt01acsummaryall, c("id", "sim")), aes(x = id, y = value)) +
  geom_line(aes(group = sim), size = 0.3) + geom_point(aes(group = sim), size = 0.5) + 
  facet_grid(variable ~ sim, scales = "free_y", labeller = label_parsed) +
  xlab("5-cycle repetition number") + theme(axis.title.y = element_blank())
pdf("br1sim1fig5Tt01ac.pdf", width = 6, height = 3)
br1sim1fig5Tt01ac
dev.off()


# Case B: failure times earlier than expected

# Case C: failure times later than expected

# Case D: varying failure time behaviour

#