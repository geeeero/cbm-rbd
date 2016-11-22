# --------------------------------------------------------------------------- #
# CBM systems paper - code for example figs - simulation study
# --------------------------------------------------------------------------- #

source("cbm-sys-other.R")

# sim 1: data fits to priors
# -------------------

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

br1sim1 <- list()
for (i in 1:20){
  cat("Repetition", i, "\n")
  br1sim1[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data1[[i]], n0y0 = br1n0y0,
                            beta = br1beta, tnowstep = 0.1, hor = 4, thresh = 0.5, seqlen = 401)
}

br1sim1summary <- data.frame(id = 1:30,
                             meantfunc = sapply(br1sim5cycle20, function(res) mean(res$tfunc)),
                             meantend = sapply(br1sim5cycle20, function(res) mean(res$tend)),
                             meancostrate = sapply(br1sim5cycle20, function(res) weighted.mean(res$costrate, res$tend)))

br1sim1fig1 <- ggplot(melt(br1sim1summary, "id"), aes(x = id, y = value)) + 
  geom_line(aes(colour = variable, group = variable)) + geom_point(aes(colour = variable, group = variable))
br1sim1fig1


