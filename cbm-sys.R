# --------------------------------------------------------------------------- #
# CBM systems paper - code for example figs
# --------------------------------------------------------------------------- #

source("weibull-sys-functions.R")
source("plotdefs.R")

# full system for prior
br <- graph.formula(s -- M -- C1:C2:C3:C4, P1:P2:P3:P4 -- t,
                    C1 -- P1, C2 -- P2, C3 -- P3, C4 -- P4, s -- H -- P3:P4)
brctypes <- list("C"=c("C1", "C2", "C3", "C4"), "H"=c("H"), "M"=c("M"), "P"=c("P1", "P2", "P3", "P4"))
br <- setCompTypes(br, brctypes)
brsysign <- computeSystemSurvivalSignature(br)

# tests # order is C H M P
brbeta <- c(2, 1.2, 2.5, 1.5)
brmttf <- c(8, 2, 5, 3)
brn0 <- c(1, 1, 2, 1)
brn0y0 <- list(C=c(brn0[1], failuretolambda(brmttf[1], brbeta[1])),
               H=c(brn0[2], failuretolambda(brmttf[2], brbeta[2])),
               M=c(brn0[3], failuretolambda(brmttf[3], brbeta[3])),
               P=c(brn0[4], failuretolambda(brmttf[4], brbeta[4])))
# system (prior) predictive at tnow = 0
brfts0 <- list(C=NULL, H=NULL, M=NULL, P=NULL)
sysrelnow(brsysign, brn0y0, brbeta, brfts0, tnow=0, t=1)
sysrelnowvec(brsysign, brn0y0, brbeta, brfts0, tnow=0, tvec=1:20)
sysrelnowhor(brsysign, brn0y0, brbeta, brfts0, tnow=0, hor=10, seqlen=5)
brtnow0 <- sysrelnowhor(brsysign, brn0y0, brbeta, brfts0, tnow=0, hor=10)
qplot(brtnow0$t, brtnow0$rel)

# posterior, assuming C2, C3, P2, P3 have failed
brfts8 <- list(c(6, 7), NULL, NULL, c(3, 4))
br8 <- induced_subgraph(br, vids=c("s", "M", "C1", "P1", "C4", "P4", "H", "t"))
br8sysign <- computeSystemSurvivalSignature(br8)
sysrelnow(br8sysign, brn0y0, brbeta, brfts8, tnow=8, t=8)
sysrelnowvec(br8sysign, brn0y0, brbeta, brfts8, tnow=8, tvec=8:20)
sysrelnowhor(br8sysign, brn0y0, brbeta, brfts8, tnow=8, hor=10, seqlen=5)
brtnow8 <- sysrelnowhor(br8sysign, brn0y0, brbeta, brfts8, tnow=8, hor=10)
qplot(brtnow8$t, brtnow8$rel)

# posterior, assuming C2, C3, P2, P3 and H have failed -> take out elements corresponding to H
brfts8H <- list(c(6, 7), 8, NULL, c(3, 4))
br8H <- induced_subgraph(br, vids=c("s", "M", "C1", "P1", "C4", "P4", "t"))
br8Hsysign <- computeSystemSurvivalSignature(br8H)
sysrelnow(br8Hsysign, brn0y0[-2], brbeta[-2], brfts8H[-2], tnow=8, t=8)
sysrelnowvec(br8Hsysign, brn0y0[-2], brbeta[-2], brfts8H[-2], tnow=8, tvec=8:20)
sysrelnowhor(br8Hsysign, brn0y0[-2], brbeta[-2], brfts8H[-2], tnow=8, hor=10, seqlen=5)
brtnow8H <- sysrelnowhor(br8Hsysign, brn0y0[-2], brbeta[-2], brfts8[-2], tnow=8, hor=10)
qplot(brtnow8H$t, brtnow8H$rel)

# component prior predictive reliability functions to visualize prior assumptions
comphor <- 10; complen <- 201
onecomp <- graph.formula(s -- 1 -- t)
Ccomp <- setCompTypes(onecomp, list(C = "1"))
Hcomp <- setCompTypes(onecomp, list(H = "1"))
Mcomp <- setCompTypes(onecomp, list(M = "1"))
Pcomp <- setCompTypes(onecomp, list(P = "1"))
Crel <- sysrelnowhor(computeSystemSurvivalSignature(Ccomp), brn0y0[1], brbeta[1], brfts0[1],
                     tnow = 0, hor = comphor, seqlen = complen)
Hrel <- sysrelnowhor(computeSystemSurvivalSignature(Hcomp), brn0y0[2], brbeta[2], brfts0[2],
                     tnow = 0, hor = comphor, seqlen = complen)
Mrel <- sysrelnowhor(computeSystemSurvivalSignature(Mcomp), brn0y0[3], brbeta[3], brfts0[3],
                     tnow = 0, hor = comphor, seqlen = complen)
Prel <- sysrelnowhor(computeSystemSurvivalSignature(Pcomp), brn0y0[4], brbeta[4], brfts0[4],
                     tnow = 0, hor = comphor, seqlen = complen)
compdf <- rbind(data.frame(Crel, comp = "C"), data.frame(Hrel, comp = "H"),
                data.frame(Mrel, comp = "M"), data.frame(Prel, comp = "P"))
# in different facets
comppriorfig1 <- ggplot(compdf, aes(x = t, y = rel)) + geom_line(aes(group = comp)) + 
  ylab("Reliability") + scale_x_continuous(breaks=seq(0, 10, by=2), minor_breaks=0:10) +
  facet_wrap(~ comp, nrow = 2, scales = "free_y") + coord_cartesian(ylim = c(0, 1))
pdf("comppriorfig1.pdf", width = 6, height = 4)
comppriorfig1
dev.off()
# as different linetypes
comppriorfig2 <- ggplot(compdf, aes(x = t, y = rel)) + geom_line(aes(linetype = comp)) + 
  ylab("Reliability") + scale_x_continuous(breaks=seq(0, 10, by=2), minor_breaks=0:10) +
  guides(linetype = guide_legend(title=NULL)) + coord_cartesian(ylim = c(0, 1))
pdf("comppriorfig2.pdf", width = 6, height = 4)
comppriorfig2
dev.off()

# -------------------------------------------------------------------------------------------------

# to plot sysrelnow history
#brcompfts <- list(C1 = NA, C2 = 6, C3 = 7, C4 = NA, H = NA, M = NA, P1 = NA, P2 = 3, P3 = 4, P4 = NA)
brcompfts <- list(C1 = NA, C2 = 6, C3 = 7, C4 = NA, H = 8,  M = NA, P1 = NA, P2 = 3, P3 = 4, P4 = NA)
# failure history plot
# turn compfts into a data.frame
brfstdf <- compfts2df(brcompfts)
# plot directly
plotfts(brcompfts)
#plotfts(brcompfts,10) 
#plotfts(brcompfts,6) 

brhist1 <- sysrelnowhist(br, brctypes, brcompfts, brn0y0, brbeta, c(0,8), hor=10)
#brhist1 <- sysrelnowhist(br, brctypes, brcompfts, brn0y0, brbeta, c(0:8), hor=10)
brhist1$tnow <- as.factor(brhist1$tnow)
histfig1 <- ggplot(brhist1, aes(x = t, y = rel)) + ylim(0, 1) +
  geom_line(aes(group = tnow, colour = tnow)) + xlab(expression(t)) + ylab(expression(R(t)))
histfig1

# sysrelnowhist does the right thing
#plot(brhist1$t, brhist1$rel, type="l")
#lines(brtnow0$t, brtnow0$rel, col=2)
#lines(brtnow8$t, brtnow8$rel, col=2)

# what is the difference when parameters do not get updated?
brhist1prio <- sysrelnowhist(br, brctypes, brcompfts, brn0y0, brbeta, c(0,8), hor=10, prior = TRUE)
#brhist1prio <- sysrelnowhist(br, brctypes, brcompfts, brn0y0, brbeta, c(0:8), hor=10, prior = TRUE)
#brhist1prio$tnow <- as.factor(brhist1prio$tnow)
histfig2 <- ggplot(rbind(data.frame(updated = "yes", brhist1),
                         data.frame(updated = "no", brhist1prio)), aes(x = t, y = rel)) + ylim(0, 1) +
  geom_line(aes(colour = tnow, linetype = updated)) + xlab(expression(t)) + ylab(expression(R(t)))
histfig2

# unit cost rate gnow
brgnow0 <- gnowhor(brsysign, brn0y0, brbeta, brfts0, tnow=0, hor=10, seqlen=201)
qplot(brgnow0$tau, brgnow0$gnow, geom="line") + coord_cartesian(ylim=c(0,1))
brgnow8 <- gnowhor(br8sysign, brn0y0, brbeta, brfts8, tnow=8, hor=10, seqlen=201)
qplot(brgnow8$tau, brgnow8$gnow, geom="line") + coord_cartesian(ylim=c(0,1))

# one component plot
oc <- graph.formula(s -- M -- t)
occtypes <- list("M"=c("M"))
oc <- setCompTypes(oc, occtypes)
ocsysign <- computeSystemSurvivalSignature(oc)
ocbeta <- 2.5
ocn0y0 <- list(M = c(1, failuretolambda(10, ocbeta)))
ocgnow0 <- gnowhor(ocsysign, ocn0y0, ocbeta, list(M = NULL), tnow=0, hor=20, seqlen=201)
ggplot(melt(ocgnow0, id="tau"), aes(x=tau, y=value)) + geom_line(aes(colour = variable))  + coord_cartesian(ylim=c(0,1))
ocgnow5 <- gnowhor(ocsysign, ocn0y0, ocbeta, list(M = NULL), tnow=5, hor=20, seqlen=201)
ggplot(melt(ocgnow5, id="tau"), aes(x=tau, y=value)) + geom_line(aes(colour = variable))  + coord_cartesian(ylim=c(0,1))

# gnow & sysrelnow histories
brghist1 <- gnowhist(br, brctypes, brcompfts, brn0y0, brbeta, c(0,2,4,6,8), hor = 10, seqlen = 1001)
#brghist2 <- gnowhist(br, brctypes, brcompfts, brn0y0, brbeta, c(0,8), hor=4, seqlen = 401)
brghist1$tnow <- as.factor(brghist1$tnow)

brghist1a <- brghist1
names(brghist1a)[c(3,4)] <- c("Reliability", "Unit cost rate")
brghist1a <- melt(brghist1a, c("tnow", "t"))
brghist1a <- subset(brghist1a, !is.na(value))
ghistfig1 <- ggplot(brghist1a, aes(x = t, y = value)) + coord_cartesian(ylim = c(0, 2.9)) +
  geom_line(aes(colour = tnow, linetype = variable)) + xlab(expression(t)) +
  #geom_line(aes(colour = tnow)) + xlab(expression(t)) + facet_wrap(~variable, nrow=2, scales = "free_y") +
  ylab(expression(paste(g^(t[now]), (t), " and ", R[sys]^(t[now]), (t)))) +
  #theme(legend.position = 'bottom', legend.direction = 'horizontal') +
  guides(colour=guide_legend(title=expression(t[now])), linetype=guide_legend(title=NULL)) +
  scale_x_continuous(breaks=seq(0, 18, by=2), minor_breaks=0:18)
pdf("ghistfig1.pdf", width = 6, height = 3)
ghistfig1
dev.off()

ghistfig2r <- ggplot(subset(brghist1a, variable == "Reliability"), aes(x = t, y = value)) +
  geom_line(aes(colour = tnow)) + coord_cartesian(xlim = c(0, 18), ylim = c(0, 1)) +
  xlab(expression(t)) + ylab(expression(paste(R[sys]^(t[now]), (t)))) + guides(colour = 'none') + 
  scale_x_continuous(breaks = seq(0, 18, by = 2), minor_breaks = 0:18) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), minor_breaks = seq(0, 1, by = 0.1)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
ghistfig2g <- ggplot(subset(brghist1a, variable == "Unit cost rate"), aes(x = t, y = value)) +
  geom_line(aes(colour = tnow)) + coord_cartesian(xlim = c(0, 18), ylim = c(0, 3)) +
  xlab(expression(t)) + ylab(expression(paste(g[sys]^(t[now]), (t)))) +
  guides(colour=guide_legend(title=expression(t[now]))) +
  theme(legend.position = 'bottom', legend.direction = 'horizontal') +
  scale_x_continuous(breaks = seq(0, 18, by = 2), minor_breaks = 0:18) +
  scale_y_continuous(breaks = seq(0, 3, by = 0.5), minor_breaks = seq(0, 3, by = 0.25))
pdf("ghistfig2.pdf", width = 6, height = 4)
grid.arrange(ghistfig2r, ghistfig2g, nrow = 2, heights = c(1, 1.5))
dev.off()

#brtaus1 <- taustarhist(br, brctypes, brcompfts, brn0y0, brbeta, c(0,2,4,6,8), hor=10)
#brtaus1 <- taustarhist(br, brctypes, brcompfts, brn0y0, brbeta, c(0,8), hor=2, seqlen=201)
#brtaus1 <- taustarhist(br, brctypes, brcompfts, brn0y0, brbeta, seq(0, 8, by=0.25), hor=5, seqlen=501)
#ggplot(melt(brtaus1, "tnow"), aes(x = tnow, y = value)) + xlab(expression(t[now])) + ylab("") +
#  geom_line(aes(colour = variable, group = variable)) + guides(colour = guide_legend(title=NULL)) +
#  geom_point(aes(colour = variable, group = variable), size = 0.15) + facet_wrap(~ variable, nrow = 2, scales = "free_y")

brtaus1fine <- taustarhist(br, brctypes, brcompfts, brn0y0, brbeta, seq(0,8,by=0.1), hor=4, seqlen=401)
tauhistfig1 <- ggplot(melt(brtaus1fine, "tnow"), aes(x = tnow, y = value)) + xlab(expression(t[now])) +
  geom_line(aes(colour = variable, group = variable)) + ylab(element_blank()) +
  geom_point(aes(colour = variable, group = variable), size = 0.5) + guides(colour = guide_legend(title=NULL)) +
  scale_x_continuous(breaks=seq(0, 8, by=2), minor_breaks=0:8)
pdf("tauhistfig1.pdf", width = 6, height = 3)
tauhistfig1
dev.off()

tauhistfig2 <- ggplot(melt(brtaus1fine, "tnow"), aes(x = tnow, y = value)) + xlab(expression(t[now])) +
  geom_line(aes(colour = variable, group = variable)) + ylab("") + guides(colour = guide_legend(title=NULL)) +
  geom_point(aes(colour = variable, group = variable), size = 0.15) +
#  facet_grid(variable ~ ., scales = "free_y")
  facet_wrap(~ variable, nrow = 2, scales = "free_y")
#  scale_x_continuous(breaks=seq(0, 8, by=2), minor_breaks=0:8)
pdf("tauhistfig2.pdf", width = 6, height = 4)
tauhistfig2
dev.off()

# compare with prior (just dynamic, not adaptive, i.e. no updating of weibull parameters)
brtaus1finepr <- taustarhist(br, brctypes, brcompfts, brn0y0, brbeta, seq(0,8,by=0.1), hor=4, seqlen=401, prior=T)
brtaus1prpo <- rbind(data.frame(melt(brtaus1fine, "tnow"), variable2 = "update"),
                     data.frame(melt(brtaus1finepr, "tnow"), variable2 = "no update"))

tauhistfig3 <- ggplot(brtaus1prpo, aes(x = tnow, y = value)) + xlab(expression(t[now])) +
  geom_line(aes(colour = variable2)) + ylab("") + guides(colour = guide_legend(title=NULL)) +
  geom_point(aes(colour = variable2), size = 0.15) +
  #  facet_grid(variable ~ ., scales = "free_y")
  facet_wrap(~ variable, nrow = 2, scales = "free_y")
pdf("tauhistfig3.pdf", width = 6, height = 4)
tauhistfig3
dev.off()

# -------------------------------------------------------------------------

# earlier failures
brcompfts2 <- list(C1 = NA, C2 = 1, C3 = 2, C4 = NA, H = 8,  M = NA, P1 = NA, P2 = 0.5, P3 = 1.5, P4 = NA)
plotfts(brcompfts2, 8)

# gnow & sysrelnow histories
brghist2 <- gnowhist(br, brctypes, brcompfts2, brn0y0, brbeta, seq(0, 2, by = 0.5), hor = 10, seqlen = 1001)
brghist2$tnow <- as.factor(brghist2$tnow)
brghist2a <- brghist2
names(brghist2a)[c(3,4)] <- c("Reliability", "Unit cost rate")
brghist2a <- melt(brghist2a, c("tnow", "t"))
brghist2a <- subset(brghist2a, !is.na(value))
ghist2fig1 <- ggplot(brghist2a, aes(x = t, y = value)) + coord_cartesian(ylim = c(0, 1.25)) +
  geom_line(aes(colour = tnow, linetype = variable)) + xlab(expression(t)) +
  #geom_line(aes(colour = tnow)) + xlab(expression(t)) + facet_wrap(~variable, nrow=2, scales = "free_y") +
  ylab(expression(paste(g^(t[now]), (t), " and ", R[sys]^(t[now]), (t)))) +
  #theme(legend.position = 'bottom', legend.direction = 'horizontal') +
  guides(colour=guide_legend(title=expression(t[now])), linetype=guide_legend(title=NULL)) +
  scale_x_continuous(breaks=seq(0, 18, by=2), minor_breaks=0:18)
pdf("ghist2fig1.pdf", width = 6, height = 3)
ghist2fig1
dev.off()

brtaus2fine <- taustarhist(br, brctypes, brcompfts2, brn0y0, brbeta, seq(0, 2, by=0.05), hor=4, seqlen=401)
tauhist2fig2 <- ggplot(melt(brtaus2fine, "tnow"), aes(x = tnow, y = value)) + xlab(expression(t[now])) +
  geom_line(aes(colour = variable, group = variable)) + ylab("") + guides(colour = guide_legend(title=NULL)) +
  geom_point(aes(colour = variable, group = variable), size = 0.15) +
  #  facet_grid(variable ~ ., scales = "free_y")
  facet_wrap(~ variable, nrow = 2, scales = "free_y")
#  scale_x_continuous(breaks=seq(0, 8, by=2), minor_breaks=0:8)
pdf("tauhist2fig2.pdf", width = 6, height = 4)
tauhist2fig2
dev.off()

brtaus2finepr <- taustarhist(br, brctypes, brcompfts2, brn0y0, brbeta, seq(0, 2, by=0.05), hor=4, seqlen=401, prior=T)
brtaus2prpo <- rbind(data.frame(melt(brtaus2fine, "tnow"), variable2 = "update"),
                     data.frame(melt(brtaus2finepr, "tnow"), variable2 = "no update"))
tauhist2fig3 <- ggplot(brtaus2prpo, aes(x = tnow, y = value)) + xlab(expression(t[now])) +
  geom_line(aes(colour = variable2)) + ylab("") + guides(colour = guide_legend(title=NULL)) +
  geom_point(aes(colour = variable2), size = 0.15) +
  #  facet_grid(variable ~ ., scales = "free_y")
  facet_wrap(~ variable, nrow = 2, scales = "free_y")
pdf("tauhist2fig3.pdf", width = 6, height = 4)
tauhist2fig3
dev.off()


# -------------------------------------------------------------------------

# later failures
brcompfts3 <- list(C1 = NA, C2 = 10, C3 = 12, C4 = NA, H = 14,  M = NA, P1 = NA, P2 = 8, P3 = 9, P4 = NA)
plotfts(brcompfts3, 14)

# gnow & sysrelnow histories
brghist3 <- gnowhist(br, brctypes, brcompfts3, brn0y0, brbeta, seq(0, 14, by = 2), hor = 10, seqlen = 1001)
brghist3$tnow <- as.factor(brghist3$tnow)
brghist3a <- brghist3
names(brghist3a)[c(3,4)] <- c("Reliability", "Unit cost rate")
brghist3a <- melt(brghist3a, c("tnow", "t"))
brghist3a <- subset(brghist3a, !is.na(value))
ghist3fig1 <- ggplot(brghist3a, aes(x = t, y = value)) + coord_cartesian(ylim = c(0, 2.5)) +
  geom_line(aes(colour = tnow, linetype = variable)) + xlab(expression(t)) +
  #geom_line(aes(colour = tnow)) + xlab(expression(t)) + facet_wrap(~variable, nrow=2, scales = "free_y") +
  ylab(expression(paste(g^(t[now]), (t), " and ", R[sys]^(t[now]), (t)))) +
  #theme(legend.position = 'bottom', legend.direction = 'horizontal') +
  guides(colour=guide_legend(title=expression(t[now])), linetype=guide_legend(title=NULL)) +
  scale_x_continuous(breaks=seq(0, 24, by=2), minor_breaks=0:24)
pdf("ghist3fig1.pdf", width = 6, height = 3)
ghist3fig1
dev.off()

brtaus3fine <- taustarhist(br, brctypes, brcompfts3, brn0y0, brbeta, seq(0, 14, by=0.2), hor=4, seqlen=401)
tauhist3fig2 <- ggplot(melt(brtaus3fine, "tnow"), aes(x = tnow, y = value)) + xlab(expression(t[now])) +
  geom_line(aes(colour = variable, group = variable)) + ylab("") + guides(colour = guide_legend(title=NULL)) +
  geom_point(aes(colour = variable, group = variable), size = 0.15) +
  #  facet_grid(variable ~ ., scales = "free_y")
  facet_wrap(~ variable, nrow = 2, scales = "free_y")
#  scale_x_continuous(breaks=seq(0, 8, by=2), minor_breaks=0:8)
pdf("tauhist3fig2.pdf", width = 6, height = 4)
tauhist3fig2
dev.off()

brtaus3finepr <- taustarhist(br, brctypes, brcompfts3, brn0y0, brbeta, seq(0, 14, by=0.2), hor=4, seqlen=401, prior=T)
brtaus3prpo <- rbind(data.frame(melt(brtaus3fine, "tnow"), variable2 = "update"),
                     data.frame(melt(brtaus3finepr, "tnow"), variable2 = "no update"))
tauhist3fig3 <- ggplot(brtaus3prpo, aes(x = tnow, y = value)) + xlab(expression(t[now])) +
  geom_line(aes(colour = variable2)) + ylab("") + guides(colour = guide_legend(title=NULL)) +
  geom_point(aes(colour = variable2), size = 0.15) +
  #  facet_grid(variable ~ ., scales = "free_y")
  facet_wrap(~ variable, nrow = 2, scales = "free_y")
pdf("tauhist3fig3.pdf", width = 6, height = 4)
tauhist3fig3
dev.off()


# -------------------------------------------------------------------------

# simulate one operational cycle
# new compfts with failure times for all components such that sim1cycle definitively terminates
brcompftssim <- list(C1 = 10, C2 = 6, C3 = 7, C4 = 10, H = 8,  M = 12, P1 = 10, P2 = 3, P3 = 4, P4 = 10)

brsim1 <- sim1cycle(sys = br, ctypes = brctypes, compfts = brcompftssim, n0y0 = brn0y0, beta = brbeta,
                    tnowstep = 0.1, hor = 4, tprep = 0.5, trepa = 0, seqlen = 401, prior = FALSE)
brsim1pr <- sim1cycle(sys = br, ctypes = brctypes, compfts = brcompftssim, n0y0 = brn0y0, beta = brbeta,
                      tnowstep = 0.1, hor = 4, tprep = 0.5, trepa = 0, seqlen = 401, prior = TRUE)

# quick sexample where system fails at time 3.7
brcompftssim0 <- list(C1 = 10, C2 = 6, C3 = 7, C4 = 10, H = 3.6,  M = 3.7, P1 = 10, P2 = 3, P3 = 4, P4 = 10)
brsim0 <- sim1cycle(sys = br, ctypes = brctypes, compfts = brcompftssim0, n0y0 = brn0y0, beta = brbeta,
                    tnowstep = 0.25, hor = 4, tprep = 1, trepa = 0, seqlen = 401, prior = FALSE)

# simulate several operational cycles
brcompftssim2 <- list(C1 = c(10, 10), C2 = c(6, 6), C3 = c(7, 7), C4 = c(10, 10), H = c(8, 8),
                      M = c(12, 12), P1 = c(10, 10), P2 = c(3, 3), P3 = c(4, 4), P4 = c(10, 10))
brsim2pr <- simNcycle(sys = br, ctypes = brctypes, compfts = brcompftssim2, n0y0 = brn0y0, beta = brbeta,
                      tnowstep = 0.5, hor = 4, tprep = 0.5, trepa = 0, seqlen = 401, prior = TRUE)
brsim2pr <- simNcycle(sys = br, ctypes = brctypes, compfts = brcompftssim2, n0y0 = brn0y0, beta = brbeta,
                      tnowstep = 0.5, hor = 4, tprep = 0.5, trepa = 0, seqlen = 401, prior = TRUE, cycleupdate = FALSE)

# simulate Weibull failure times according to prior parameter choices (note different parametrization!)
set.seed(1328)
ncycles <- 5
C1sim1 <- rweibull(ncycles, shape = brbeta[1], scale = (failuretolambda(brmttf[1], brbeta[1]))^(1/brbeta[1]))
C2sim1 <- rweibull(ncycles, shape = brbeta[1], scale = (failuretolambda(brmttf[1], brbeta[1]))^(1/brbeta[1]))
C3sim1 <- rweibull(ncycles, shape = brbeta[1], scale = (failuretolambda(brmttf[1], brbeta[1]))^(1/brbeta[1]))
C4sim1 <- rweibull(ncycles, shape = brbeta[1], scale = (failuretolambda(brmttf[1], brbeta[1]))^(1/brbeta[1]))
Hsim1  <- rweibull(ncycles, shape = brbeta[2], scale = (failuretolambda(brmttf[2], brbeta[2]))^(1/brbeta[2]))
Msim1  <- rweibull(ncycles, shape = brbeta[3], scale = (failuretolambda(brmttf[3], brbeta[3]))^(1/brbeta[3]))
P1sim1 <- rweibull(ncycles, shape = brbeta[4], scale = (failuretolambda(brmttf[4], brbeta[4]))^(1/brbeta[4]))
P2sim1 <- rweibull(ncycles, shape = brbeta[4], scale = (failuretolambda(brmttf[4], brbeta[4]))^(1/brbeta[4]))
P3sim1 <- rweibull(ncycles, shape = brbeta[4], scale = (failuretolambda(brmttf[4], brbeta[4]))^(1/brbeta[4]))
P4sim1 <- rweibull(ncycles, shape = brbeta[4], scale = (failuretolambda(brmttf[4], brbeta[4]))^(1/brbeta[4]))
brcompftssim1 <- list(C1 = C1sim1, C2 = C2sim1, C3 = C3sim1, C4 = C4sim1, H = Hsim1,
                      M = Msim1, P1 = P1sim1, P2 = P2sim1, P3 = P3sim1, P4 = P4sim1)
brsimN5 <- simNcycle(sys = br, ctypes = brctypes, compfts = brcompftssim1, n0y0 = brn0y0, beta = brbeta,
                     tnowstep = 0.1, hor = 4, tprep = 0.5, trepa = 0, seqlen = 401)
sapply(brsimN5, function(res) res$tend)
sapply(brsimN5, function(res) res$downtime)
sapply(brsimN5, function(res) res$costrate)
# overall costrate
weighted.mean(sapply(brsimN5, function(res) res$costrate), sapply(brsimN5, function(res) res$tend))
brsimN5taudf <- rbind(data.frame(cycle = 1, brsimN5[[1]]$res),
                      data.frame(cycle = 2, brsimN5[[2]]$res),
                      data.frame(cycle = 3, brsimN5[[3]]$res),
                      data.frame(cycle = 4, brsimN5[[4]]$res),
                      data.frame(cycle = 5, brsimN5[[5]]$res))
brsimN5taudf$cycle <- as.factor(brsimN5taudf$cycle)
brsimN5fig1 <- ggplot(brsimN5taudf, aes(x = tnow, y = taustar)) + geom_line(aes(colour = cycle, group = cycle)) +
  xlab(expression(t[now])) + ylab(expression(paste(tau["*"]^(t[now]), (t[now])))) +
  guides(colour = guide_legend(title="Cycle"))
pdf("brsimN5fig1.pdf", width = 5, height = 3)
brsimN5fig1
dev.off()

# -------------------------------------------------------------------------




#