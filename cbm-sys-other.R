# --------------------------------------------------------------------------- #
# CBM systems paper - code for example figs - other prior
# --------------------------------------------------------------------------- #

source("weibull-sys-functions.R")
source("cbm-sim.R")
source("plotdefs.R")

# full system for prior
br <- graph.formula(s -- M -- C1:C2:C3:C4, P1:P2:P3:P4 -- t,
                    C1 -- P1, C2 -- P2, C3 -- P3, C4 -- P4, s -- H -- P3:P4)
brctypes <- list("C"=c("C1", "C2", "C3", "C4"), "H"=c("H"), "M"=c("M"), "P"=c("P1", "P2", "P3", "P4"))
br <- setCompTypes(br, brctypes)
brsysign <- computeSystemSurvivalSignature(br)

# order is C H M P
br1beta <- c(2, 1, 2.5, 1.5)
br1mttf <- c(3, 10, 8, 5)
br1n0 <- c(1, 1, 1, 1) # weak prior
br1n0y0 <- list(C=c(br1n0[1], failuretolambda(br1mttf[1], br1beta[1])),
                H=c(br1n0[2], failuretolambda(br1mttf[2], br1beta[2])),
                M=c(br1n0[3], failuretolambda(br1mttf[3], br1beta[3])),
                P=c(br1n0[4], failuretolambda(br1mttf[4], br1beta[4])))
brfts0 <- list(C=NULL, H=NULL, M=NULL, P=NULL)

# component prior predictive reliability functions to visualize prior assumptions
comphor <- 10; complen <- 201
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
# in different facets
compprior1fig1 <- ggplot(compdf1, aes(x = t, y = rel)) + geom_line(aes(group = comp)) + 
  ylab("Reliability") + scale_x_continuous(breaks=seq(0, 10, by=2), minor_breaks=0:10) +
  facet_wrap(~ comp, nrow = 2, scales = "free_y") + coord_cartesian(ylim = c(0, 1))
pdf("compprior1fig1.pdf", width = 6, height = 4)
compprior1fig1
dev.off()
# as different linetypes
compprior1fig2 <- ggplot(compdf1, aes(x = t, y = rel)) + geom_line(aes(linetype = comp)) + 
  ylab("Reliability") + scale_x_continuous(breaks=seq(0, 10, by=2), minor_breaks=0:10) +
  guides(linetype = guide_legend(title=NULL)) + coord_cartesian(ylim = c(0, 1))
pdf("compprior1fig2.pdf", width = 6, height = 4)
compprior1fig2
dev.off()

# -------------------------------------------------------------------------------------------------

# exemplary failure history
brcompfts <- list(C1 = NA, C2 = 6, C3 = 7, C4 = NA, H = 8,  M = NA, P1 = NA, P2 = 3, P3 = 4, P4 = NA)
# failure history plot
brfstdf <- compfts2df(brcompfts)
plotfts(brcompfts)
#plotfts(brcompfts,10) 
#plotfts(brcompfts,6) 

# gnow & sysrelnow histories
br1ghist1 <- gnowhist(br, brctypes, brcompfts, br1n0y0, br1beta, c(0,2,4,6,8), hor = 10, seqlen = 1001)
br1ghist1$tnow <- as.factor(br1ghist1$tnow)
br1ghist1a <- br1ghist1
names(br1ghist1a)[c(3,4)] <- c("Reliability", "Unit cost rate")
br1ghist1a <- melt(br1ghist1a, c("tnow", "t"))
br1ghist1a <- subset(br1ghist1a, !is.na(value))
ghist1fig2r <- ggplot(subset(br1ghist1a, variable == "Reliability"), aes(x = t, y = value)) +
  geom_line(aes(colour = tnow)) + coord_cartesian(xlim = c(0, 18), ylim = c(0, 1)) +
  xlab(expression(t)) + ylab(expression(paste(R[sys]^(t[now]), (t)))) + guides(colour = 'none') + 
  scale_x_continuous(breaks = seq(0, 18, by = 2), minor_breaks = 0:18) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), minor_breaks = seq(0, 1, by = 0.1)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ghist1fig2g <- ggplot(subset(br1ghist1a, variable == "Unit cost rate"), aes(x = t, y = value)) +
  geom_line(aes(colour = tnow)) + coord_cartesian(xlim = c(0, 18), ylim = c(0, 2.2)) +
  xlab(expression(t)) + ylab(expression(paste(g[sys]^(t[now]), (t)))) +
  guides(colour=guide_legend(title=expression(t[now]))) +
  theme(legend.position = 'bottom', legend.direction = 'horizontal') +
  scale_x_continuous(breaks = seq(0, 18, by = 2), minor_breaks = 0:18) +
  scale_y_continuous(breaks = seq(0, 3, by = 0.5), minor_breaks = seq(0, 3, by = 0.25))
pdf("ghist1fig2.pdf", width = 6, height = 4)
grid.arrange(ghist1fig2r, ghist1fig2g, nrow = 2, heights = c(1, 1.5))
dev.off()

# -------------------------------------------------------------------------------------------------

br1taus1fine <- taustarhist(br, brctypes, brcompfts, br1n0y0, br1beta, seq(0,8,by=0.1), hor=4, seqlen=401)
tauhist1fig2 <- ggplot(melt(br1taus1fine, "tnow"), aes(x = tnow, y = value)) + xlab(expression(t[now])) +
  geom_line(aes(colour = variable, group = variable)) + ylab("") + guides(colour = guide_legend(title=NULL)) +
  geom_point(aes(colour = variable, group = variable), size = 0.15) +
  facet_wrap(~ variable, nrow = 2, scales = "free_y")
pdf("tauhist1fig2.pdf", width = 6, height = 4)
tauhist1fig2
dev.off()
br1taus1fine[40:41,]
br1taus1finepr <- taustarhist(br, brctypes, brcompfts, br1n0y0, br1beta, seq(0,8,by=0.1), hor=4, seqlen=401, prior=T)
br1taus1prpo <- rbind(data.frame(melt(br1taus1fine, "tnow"), variable2 = "update"),
                      data.frame(melt(br1taus1finepr, "tnow"), variable2 = "no update"))
tauhist1fig3 <- ggplot(br1taus1prpo, aes(x = tnow, y = value)) + xlab(expression(t[now])) +
  geom_line(aes(colour = variable2)) + ylab("") + guides(colour = guide_legend(title=NULL)) +
  geom_point(aes(colour = variable2), size = 0.15) +
  facet_wrap(~ variable, nrow = 2, scales = "free_y")
pdf("tauhist1fig3.pdf", width = 6, height = 4)
tauhist1fig3
dev.off()

# -------------------------------------------------------------------------

# earlier failures
brcompfts2 <- list(C1 = NA, C2 = 1, C3 = 2, C4 = NA, H = 8,  M = NA, P1 = NA, P2 = 0.5, P3 = 1.5, P4 = NA)
#plotfts(brcompfts2, 8)
br1taus2fine <- taustarhist(br, brctypes, brcompfts2, br1n0y0, br1beta, seq(0, 2, by=0.05), hor=4, seqlen=401)
br1taus2finepr <- taustarhist(br, brctypes, brcompfts2, br1n0y0, br1beta, seq(0, 2, by=0.05), hor=4, seqlen=401, prior=T)
br1taus2prpo <- rbind(data.frame(melt(br1taus2fine, "tnow"), variable2 = "update"),
                      data.frame(melt(br1taus2finepr, "tnow"), variable2 = "no update"))
tauhist2fig3 <- ggplot(br1taus2prpo, aes(x = tnow, y = value)) + xlab(expression(t[now])) +
  geom_line(aes(colour = variable2)) + ylab("") + guides(colour = guide_legend(title=NULL)) +
  geom_point(aes(colour = variable2), size = 0.15) +
  facet_wrap(~ variable, nrow = 2, scales = "free_y")
pdf("tauhist2fig3.pdf", width = 6, height = 4)
tauhist2fig3
dev.off()

# later failures
brcompfts3 <- list(C1 = NA, C2 = 10, C3 = 12, C4 = NA, H = 14,  M = NA, P1 = NA, P2 = 8, P3 = 9, P4 = NA)
#plotfts(brcompfts3, 14)
br1taus3fine <- taustarhist(br, brctypes, brcompfts3, br1n0y0, br1beta, seq(0, 14, by=0.2), hor=4, seqlen=401)
br1taus3finepr <- taustarhist(br, brctypes, brcompfts3, br1n0y0, br1beta, seq(0, 14, by=0.2), hor=4, seqlen=401, prior=T)
br1taus3prpo <- rbind(data.frame(melt(br1taus3fine, "tnow"), variable2 = "update"),
                      data.frame(melt(br1taus3finepr, "tnow"), variable2 = "no update"))
tauhist3fig3 <- ggplot(br1taus3prpo, aes(x = tnow, y = value)) + xlab(expression(t[now])) +
  geom_line(aes(colour = variable2)) + ylab("") + guides(colour = guide_legend(title=NULL)) +
  geom_point(aes(colour = variable2), size = 0.15) +
  facet_wrap(~ variable, nrow = 2, scales = "free_y")
pdf("tauhist3fig3.pdf", width = 6, height = 4)
tauhist3fig3
dev.off()

# -------------------------------------------------------------------------

# small simulation study
# simulate 20 5-cycle machines with priors as given
set.seed(2211)
br1sim5cycle20data <- list()
for (i in 1:20){
  br1sim5cycle20data[[i]] <- brWeibullData(5, br1beta, br1mttf)
}

br1simN51 <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data[[1]], n0y0 = br1n0y0, beta = br1beta,
                      tnowstep = 0.1, hor = 4, tprep = 0.5, trepa = 0, seqlen = 401)
br1simN51$downtime; br1simN51$tend; br1simN51$tfunc; br1simN51$costrate
br1simN51fig1 <- ggplot(br1simN51$res, aes(x = tnow, y = taustar)) + 
  geom_line(aes(colour = cycle, group = cycle)) + geom_point(aes(colour = cycle, group = cycle)) +
  xlab(expression(t[now])) + ylab(expression(paste(tau["*"]^(t[now]), (t[now])))) +
  guides(colour = guide_legend(title="Cycle")) + scale_y_continuous(breaks=seq(0, 2, by=0.5), minor_breaks=seq(0, 2, by=0.25))
pdf("br1simN51fig1.pdf", width = 5, height = 3)
br1simN51fig1
dev.off()

br1sim5cycle20 <- list()
for (i in 1:20){
  br1sim5cycle20[[i]] <- simNcycle(sys = br, ctypes = brctypes, compfts = br1sim5cycle20data[[i]], n0y0 = br1n0y0,
                                   beta = br1beta, tnowstep = 0.1, hor = 4, tprep = 0.5, trepa = 0, seqlen = 401)
}

br1sim5cycle20summary <- data.frame(id = 1:20,
                                    downtime = sapply(br1sim5cycle20, function(res) sum(res$downtime)),
                                    meantend = sapply(br1sim5cycle20, function(res) mean(res$tend)),
                                    meancostrate = sapply(br1sim5cycle20, function(res) weighted.mean(res$costrate, res$tend)))

br1sim5cycle20fig1 <- ggplot(melt(br1sim5cycle20summary, "id"), aes(x = id, y = value)) + 
  geom_line(aes(colour = variable, group = variable)) + geom_point(aes(colour = variable, group = variable))
br1sim5cycle20fig1




#