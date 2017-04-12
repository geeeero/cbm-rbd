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
load("extsim-a-objects.RData")

## apply summary functions & create joint data.frame
AsimCBMcpuSummary <- simsummary(AsimCBMcpu)
AsimCBMepuSummary <- simsummary(AsimCBMepu)
AsimCBMnpuSummary <- simsummary(AsimCBMnpu)
AsimABMepuSummary <- simsummary(AsimABMepu)
AsimABMnpuSummary <- simsummary(AsimABMnpu)
AsimCMSummary <- simsummary(AsimCM)
AsimRes <- rbind(data.frame(sim = "CBM-cpu", AsimCBMcpuSummary),
                 data.frame(sim = "CBM-epu", AsimCBMepuSummary),
                 data.frame(sim = "CBM-npu", AsimCBMnpuSummary),
                 data.frame(sim = "ABM-epu", AsimABMepuSummary),
                 data.frame(sim = "ABM-npu", AsimABMnpuSummary),
                 data.frame(sim = "CM", AsimCMSummary))
names(AsimRes)[3:5] <- c(expression(e[sys]), expression(bar(r)[sys]), expression(bar(g)))

## the plots
### with lines like before (change sizes?)
AsimPlotLines <- ggplot(melt(AsimRes, c("id", "sim")), aes(x = id, y = value)) +
  geom_line(aes(group = sim), size = 0.3) + geom_point(aes(group = sim), size = 0.5) + 
  facet_grid(variable ~ sim, scales = "free_y", labeller = label_parsed) +
  xlab("25-cycle repetition number") + theme(axis.title.y = element_blank())
pdf("AsimPlotLines.pdf", width = 6, height = 3)
AsimPlotLines
dev.off()
### one boxplot per policy
AsimPlotBoxplot <- ggplot(melt(AsimRes, c("id", "sim")), aes(x = sim, y = value, group = sim)) + stat_boxplot() +
  facet_grid(variable ~ ., scales = "free_y", labeller = label_parsed) + theme(axis.title = element_blank())
pdf("AsimPlotLines.pdf", width = 6, height = 3)
AsimPlotBoxplot
dev.off()
### other way to display the results?
#???

# boxplot test
#ggplot(melt(br1sim1Tt01acsummaryall, c("id", "sim")), aes(x = sim, y = value, group = sim)) + stat_boxplot() +
#  facet_grid(variable ~ ., scales = "free_y", labeller = label_parsed) + theme(axis.title = element_blank())

## mean and median of average unit cost
#mean(AsimCBMcpuSummary$meancostrate); median(AsimCBMcpuSummary$meancostrate) # there must be a more clever way
#mean(AsimRes$"bar(g)", by = sim) # something like this


# Case B: failure times earlier than expected

# Case C: failure times later than expected

# Case D: varying failure time behaviour

#