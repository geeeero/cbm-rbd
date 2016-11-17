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
# prior
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






# -------------------------------------------------------------------------


brprtvec <- seq(0, 20, length.out=101)

# prior system reliability function at tnow = 0
brprio <- sysrelPbox(luckobjlist = brpr, survsign = brsysign, nk = brprnk, beta = brbeta,
                     fts = brfts0, tnow = 0, tvec = brprtvec, prior = TRUE)

# component reliability functions via sysrelPbox, posterior does not work
brMprrel <- sysrelPbox(luckobjlist = list(brpr[[1]]), survsign = oneCompSurvSign("M"), nk = 1,
                       beta = brbeta[1], fts = list(NULL), tnow = 0, tvec = brprtvec, prior = TRUE)
brHprrel <- sysrelPbox(luckobjlist = list(brpr[[2]]), survsign = oneCompSurvSign("H"), nk = 1,
                       beta = brbeta[2], fts = list(NULL), tnow = 0, tvec = brprtvec, prior = TRUE)
brCprrel <- sysrelPbox(luckobjlist = list(brpr[[3]]), survsign = oneCompSurvSign("C"), nk = 1,
                       beta = brbeta[3], fts = list(NULL), tnow = 0, tvec = brprtvec, prior = TRUE)
brPprrel <- sysrelPbox(luckobjlist = list(brpr[[4]]), survsign = oneCompSurvSign("P"), nk = 1,
                       beta = brbeta[4], fts = list(NULL), tnow = 0, tvec = brprtvec, prior = TRUE)
brprdf <- rbind(data.frame(brMprrel, Part = "M", Item = "Prior"), data.frame(brHprrel, Part = "H", Item = "Prior"),
                data.frame(brCprrel, Part = "C", Item = "Prior"), data.frame(brPprrel, Part = "P", Item = "Prior"))
brprdf$Item <- ordered(brprdf$Item, levels=c("Prior", "Posterior"))
brprdf$Part <- ordered(brprdf$Part, levels=c("M", "H", "C", "P"))
fig3 <- ggplot(brprdf, aes(x = tvec)) + theme_bw() + ijarcol + ijarfill + ylim(0, 1) +  
  geom_ribbon(aes(ymin = lower, ymax = upper, group = Item, colour = Item, fill = Item), alpha = 0.5) +
  geom_line(aes(y = lower, group = Item, colour = Item)) + 
  geom_line(aes(y = upper, group = Item, colour = Item)) + 
  facet_wrap(~Part, nrow=2) + xlab(expression(t)) + ylab(expression(R(t))) + nolegend 
fig3

# component reliability functions using the prior predictive
prprtvec <- seq(0, 20, length.out=301)
prprC <- prprRluck(prprtvec, brpr[[1]], beta = brbeta[1], n0gridlen = 2000)
prprH <- prprRluck(prprtvec, brpr[[2]], beta = brbeta[2], n0gridlen = 2000)
prprM <- prprRluck(prprtvec, brpr[[3]], beta = brbeta[3], n0gridlen = 2000)
prprP <- prprRluck(prprtvec, brpr[[4]], beta = brbeta[4], n0gridlen = 2000)
prprdf <- rbind(data.frame(prprM, Part = "M", Item = "Prior"), data.frame(prprH, Part = "H", Item = "Prior"),
                data.frame(prprC, Part = "C", Item = "Prior"), data.frame(prprP, Part = "P", Item = "Prior"))
prprdf$Item <- ordered(prprdf$Item, levels=c("Prior", "Posterior"))
prprdf$Part <- ordered(prprdf$Part, levels=c("M", "H", "C", "P"))
prprdf$lower[prprdf$lower < 0] <- 0
prprdf$upper[prprdf$upper > 1] <- 1
fig4 <- ggplot(prprdf, aes(x = tvec)) + theme_bw() + ijarcol + ijarfill + ylim(0, 1) +  
  geom_ribbon(aes(ymin = lower, ymax = upper, group = Item, colour = Item, fill = Item), alpha = 0.5) +
  geom_line(aes(y = lower, group = Item, colour = Item)) + 
  geom_line(aes(y = upper, group = Item, colour = Item)) + 
  facet_wrap(~Part, nrow=2) + xlab(expression(t)) + ylab(expression(R(t))) + rightlegend #nolegend 
fig4 # warnings refer to the four instances which have been rounded down to 1, bug in ggplot?

#setEPS()
#postscript("fig4.eps",width=8,height=6)
pdf("fig4.pdf", width=8, height=6)
fig4
dev.off()


# data 1: fitting failure times ---------------------------------------------------------------- #
brfts1 <- list(c(6, 7), NULL, NULL, c(3, 4))
tnow1 <- 8

# components posterior # CHMP
# n and tnow are used for the updating!
# prprtvec <- seq(0, 20, length.out=501)
#prprtvec[1] <- 1e-5
popr1C <- poprRluck(prprtvec, brpr[[1]], beta = brbeta[1], fts = brfts1[[1]], n = brprnk[1], tnow = tnow1, n0gridlen = 2000)
popr1H <- poprRluck(prprtvec, brpr[[2]], beta = brbeta[2], fts = brfts1[[2]], n = brprnk[2], tnow = tnow1, n0gridlen = 2000)
popr1M <- poprRluck(prprtvec, brpr[[3]], beta = brbeta[3], fts = brfts1[[3]], n = brprnk[3], tnow = tnow1, n0gridlen = 2000)
popr1P <- poprRluck(prprtvec, brpr[[4]], beta = brbeta[4], fts = brfts1[[4]], n = brprnk[4], tnow = tnow1, n0gridlen = 2000)
popr1df <- rbind(data.frame(prprC, Part = "C", Item = "Prior"), data.frame(popr1C, Part = "C", Item = "Posterior"),
                 data.frame(prprH, Part = "H", Item = "Prior"), data.frame(popr1H, Part = "H", Item = "Posterior"),
                 data.frame(prprM, Part = "M", Item = "Prior"), data.frame(popr1M, Part = "M", Item = "Posterior"),
                 data.frame(prprP, Part = "P", Item = "Prior"), data.frame(popr1P, Part = "P", Item = "Posterior"))
popr1df$Item <- ordered(popr1df$Item, levels=c("Prior", "Posterior"))
popr1df$Part <- ordered(popr1df$Part, levels=c("M", "H", "C", "P"))
popr1df$lower[popr1df$lower < 0] <- 0
popr1df$upper[popr1df$upper > 1] <- 1
fig5a <- ggplot(popr1df, aes(x = tvec)) + theme_bw() + ijarcol + ijarfill + ylim(0, 1) +  
  geom_ribbon(aes(ymin = lower, ymax = upper, group = Item, colour = Item, fill = Item), alpha = 0.5) +
  geom_line(aes(y = lower, group = Item, colour = Item)) + 
  geom_line(aes(y = upper, group = Item, colour = Item)) + 
  facet_wrap(~Part, nrow=2) + xlab(expression(t)) + ylab(expression(R(t))) + rightlegend 
fig5a # warnings refer to the instances which have been rounded down to 1, bug in ggplot?

pdf("fig5a.pdf", width=8, height=6)
fig5a
dev.off()

# posterior system reliability function at tnow = 8
brpos1tvec <- seq(tnow1, 20, length.out=61)
brpost1 <- sysrelPbox(luckobjlist = brpr, survsign = brpossysign, nk = brprnk, beta = brbeta,
                      fts = brfts1, tnow = tnow1, tvec = brpos1tvec, prior = FALSE)
br1df <- rbind(data.frame(brprio, Panel = "Prior", Item = "Prior"),
               data.frame(brprio, Panel = "Prior and Posterior 1", Item = "Prior"),
               data.frame(brpost1, Panel = "Prior and Posterior 1", Item = "Posterior"))
br1df$Item <- ordered(br1df$Item, levels=c("Prior", "Posterior"))
fig5 <- ggplot(br1df, aes(x = tvec)) + theme_bw() + ijarcol + ijarfill + ylim(0, 1) +  
  geom_ribbon(aes(ymin = lower, ymax = upper, group = Item, colour = Item, fill = Item), alpha = 0.5) +
  geom_line(aes(y = lower, group = Item, colour = Item)) + 
  geom_line(aes(y = upper, group = Item, colour = Item)) + facet_wrap(~Panel, nrow=1) + 
  xlab(expression(t)) + ylab(expression(R(t))) + rightlegend
fig5

#setEPS()
#postscript("fig5.eps",width=8,height=6)
pdf("fig5.pdf", width=8, height=3)
fig5
dev.off()


# data 2: early failure times ---------------------------------------------------------------- #
brfts2 <- list(c(1, 2), NULL, NULL, c(0.25, 0.5))
tnow2 <- 2

# components posterior # CHMP
popr2C <- poprRluck(prprtvec, brpr[[1]], beta = brbeta[1], fts = brfts2[[1]], n = brprnk[1], tnow = tnow2, n0gridlen = 2000)
popr2H <- poprRluck(prprtvec, brpr[[2]], beta = brbeta[2], fts = brfts2[[2]], n = brprnk[2], tnow = tnow2, n0gridlen = 2000)
popr2M <- poprRluck(prprtvec, brpr[[3]], beta = brbeta[3], fts = brfts2[[3]], n = brprnk[3], tnow = tnow2, n0gridlen = 2000)
popr2P <- poprRluck(prprtvec, brpr[[4]], beta = brbeta[4], fts = brfts2[[4]], n = brprnk[4], tnow = tnow2, n0gridlen = 2000)
popr2df <- rbind(data.frame(prprC, Part = "C", Item = "Prior"), data.frame(popr2C, Part = "C", Item = "Posterior"),
                 data.frame(prprH, Part = "H", Item = "Prior"), data.frame(popr2H, Part = "H", Item = "Posterior"),
                 data.frame(prprM, Part = "M", Item = "Prior"), data.frame(popr2M, Part = "M", Item = "Posterior"),
                 data.frame(prprP, Part = "P", Item = "Prior"), data.frame(popr2P, Part = "P", Item = "Posterior"))
popr2df$Item <- ordered(popr2df$Item, levels=c("Prior", "Posterior"))
popr2df$Part <- ordered(popr2df$Part, levels=c("M", "H", "C", "P"))
popr2df$lower[popr2df$lower < 0] <- 0
popr2df$upper[popr2df$upper > 1] <- 1
fig6a <- ggplot(popr2df, aes(x = tvec)) + theme_bw() + ijarcol + ijarfill + ylim(0, 1) +  
  geom_ribbon(aes(ymin = lower, ymax = upper, group = Item, colour = Item, fill = Item), alpha = 0.5) +
  geom_line(aes(y = lower, group = Item, colour = Item)) + 
  geom_line(aes(y = upper, group = Item, colour = Item)) + 
  facet_wrap(~Part, nrow=2) + xlab(expression(t)) + ylab(expression(R(t))) + rightlegend 
fig6a # warnings refer to the instances which have been rounded down to 1, bug in ggplot?
pdf("fig6a.pdf", width=8, height=6)
fig6a
dev.off()

# posterior system reliability function at tnow = 2
brpos2tvec <- seq(tnow2, 20, length.out=91)
brpost2 <- sysrelPbox(luckobjlist = brpr, survsign = brpossysign, nk = brprnk, beta = brbeta,
                      fts = brfts2, tnow = tnow2, tvec = brpos2tvec, prior = FALSE)
br2df <- rbind(data.frame(brprio, Panel = "Prior", Item = "Prior"),
               data.frame(brprio, Panel = "Prior and Posterior 2", Item = "Prior"),
               data.frame(brpost2, Panel = "Prior and Posterior 2", Item = "Posterior"))
br2df$Item <- ordered(br2df$Item, levels=c("Prior", "Posterior"))
fig6 <- ggplot(br2df, aes(x = tvec)) + theme_bw() + ijarcol + ijarfill + ylim(0, 1) +  
  geom_ribbon(aes(ymin = lower, ymax = upper, group = Item, colour = Item, fill = Item), alpha = 0.5) +
  geom_line(aes(y = lower, group = Item, colour = Item)) + 
  geom_line(aes(y = upper, group = Item, colour = Item)) + facet_wrap(~Panel, nrow=1) + 
  xlab(expression(t)) + ylab(expression(R(t))) + rightlegend
fig6

#setEPS()
#postscript("fig6.eps",width=8,height=3)
pdf("fig6.pdf", width=8, height=3)
fig6
dev.off()


# data 3: late failure times ---------------------------------------------------------------- #
brfts3 <- list(c(11, 12), NULL, NULL, c(8, 9))
tnow3 <- 12

# components posterior # CHMP
popr3C <- poprRluck(prprtvec, brpr[[1]], beta = brbeta[1], fts = brfts3[[1]], n = brprnk[1], tnow = tnow3, n0gridlen = 2000)
popr3H <- poprRluck(prprtvec, brpr[[2]], beta = brbeta[2], fts = brfts3[[2]], n = brprnk[2], tnow = tnow3, n0gridlen = 2000)
popr3M <- poprRluck(prprtvec, brpr[[3]], beta = brbeta[3], fts = brfts3[[3]], n = brprnk[3], tnow = tnow3, n0gridlen = 2000)
popr3P <- poprRluck(prprtvec, brpr[[4]], beta = brbeta[4], fts = brfts3[[4]], n = brprnk[4], tnow = tnow3, n0gridlen = 2000)
popr3df <- rbind(data.frame(prprC, Part = "C", Item = "Prior"), data.frame(popr3C, Part = "C", Item = "Posterior"),
                 data.frame(prprH, Part = "H", Item = "Prior"), data.frame(popr3H, Part = "H", Item = "Posterior"),
                 data.frame(prprM, Part = "M", Item = "Prior"), data.frame(popr3M, Part = "M", Item = "Posterior"),
                 data.frame(prprP, Part = "P", Item = "Prior"), data.frame(popr3P, Part = "P", Item = "Posterior"))
popr3df$Item <- ordered(popr3df$Item, levels=c("Prior", "Posterior"))
popr3df$Part <- ordered(popr3df$Part, levels=c("M", "H", "C", "P"))
popr3df$lower[popr3df$lower < 0] <- 0
popr3df$upper[popr3df$upper > 1] <- 1
fig7a <- ggplot(popr3df, aes(x = tvec)) + theme_bw() + ijarcol + ijarfill + ylim(0, 1) +  
  geom_ribbon(aes(ymin = lower, ymax = upper, group = Item, colour = Item, fill = Item), alpha = 0.5) +
  geom_line(aes(y = lower, group = Item, colour = Item)) + 
  geom_line(aes(y = upper, group = Item, colour = Item)) + 
  facet_wrap(~Part, nrow=2) + xlab(expression(t)) + ylab(expression(R(t))) + rightlegend 
fig7a # warnings refer to the instances which have been rounded down to 1, bug in ggplot?
pdf("fig7a.pdf", width=8, height=6)
fig7a
dev.off()

# posterior system reliability function at tnow = 12
brpos3tvec <- seq(tnow3, 20, length.out=41)
brpost3 <- sysrelPbox(luckobjlist = brpr, survsign = brpossysign, nk = brprnk, beta = brbeta,
                      fts = brfts3, tnow = tnow3, tvec = brpos3tvec, prior = FALSE)
br3df <- rbind(data.frame(brprio, Panel = "Prior", Item = "Prior"),
               data.frame(brprio, Panel = "Prior and Posterior 3", Item = "Prior"),
               data.frame(brpost3, Panel = "Prior and Posterior 3", Item = "Posterior"))
br3df$Item <- ordered(br3df$Item, levels=c("Prior", "Posterior"))
fig7 <- ggplot(br3df, aes(x = tvec)) + theme_bw() + ijarcol + ijarfill + ylim(0, 1) +  
  geom_ribbon(aes(ymin = lower, ymax = upper, group = Item, colour = Item, fill = Item), alpha = 0.5) +
  geom_line(aes(y = lower, group = Item, colour = Item)) + 
  geom_line(aes(y = upper, group = Item, colour = Item)) + facet_wrap(~Panel, nrow=1) + 
  xlab(expression(t)) + ylab(expression(R(t))) + rightlegend
fig7

#setEPS()
#postscript("fig7.eps",width=8,height=3)
pdf("fig7.pdf", width=8, height=3)
fig7
dev.off()

# figure with prior and the three posterior system reliabilties in four panels
br4df <- rbind(data.frame(brprio,  Panel = "Prior", Item = "Prior"),
               data.frame(brprio,  Panel = "Failures as expected", Item = "Prior"),
               data.frame(brpost1, Panel = "Failures as expected", Item = "Posterior"),
               data.frame(brprio,  Panel = "Surprisingly early failures", Item = "Prior"),
               data.frame(brpost2, Panel = "Surprisingly early failures", Item = "Posterior"),
               data.frame(brprio,  Panel = "Surprisingly late failures", Item = "Prior"),
               data.frame(brpost3, Panel = "Surprisingly late failures", Item = "Posterior"))
br4df$Item <- ordered(br4df$Item, levels=c("Prior", "Posterior"))
br3df$Panel <- ordered(br3df$Panel, levels=c("Prior", "Failures as expected", "Surprisingly early failures",
                                             "Surprisingly late failures"))
fig8 <- ggplot(br4df, aes(x = tvec)) + theme_bw() + ijarcol + ijarfill + ylim(0, 1) +  
  geom_ribbon(aes(ymin = lower, ymax = upper, group = Item, colour = Item, fill = Item), alpha = 0.5) +
  geom_line(aes(y = lower, group = Item, colour = Item)) + 
  geom_line(aes(y = upper, group = Item, colour = Item)) + facet_wrap(~Panel, nrow=2) + 
  xlab(expression(t)) + ylab(expression(R[sys](t))) + rightlegend
fig8

#setEPS()
#postscript("fig8.eps",width=8,height=6)
pdf("fig8.pdf", width=8, height=6)
fig8
dev.off()

# plot of sysrel prior / posterior imprecision
brpost1e <- brpost1
brpost1e[,1] <- brpost1e[,1] - tnow1
brpost2e <- brpost2
brpost2e[,1] <- brpost2e[,1] - tnow2
brpost3e <- brpost3
brpost3e[,1] <- brpost3e[,1] - tnow3
br5df <- rbind(data.frame(br4df, Panel2 = "Elapsed time"),
               data.frame(brprio,  Panel = "Prior", Item = "Prior", Panel2 = "Prospective time"),
               #data.frame(brprio,  Panel = "Failures as expected", Item = "Prior", Panel2 = "Prospective time"),
               data.frame(brpost1e, Panel = "Failures as expected", Item = "Posterior", Panel2 = "Prospective time"),
               #data.frame(brprio,  Panel = "Surprisingly early failures", Item = "Prior", Panel2 = "Prospective time"),
               data.frame(brpost2e, Panel = "Surprisingly early failures", Item = "Posterior", Panel2 = "Prospective time"),
               #data.frame(brprio,  Panel = "Surprisingly late failures", Item = "Prior", Panel2 = "Prospective time"),
               data.frame(brpost3e, Panel = "Surprisingly late failures", Item = "Posterior", Panel2 = "Prospective time"))
#br4dfe <- br4df
#br4dfe$tvec <- 
#br5df <- rbind(data.frame(br4df, Panel2 = "Elapsed time"), data.frame(br4dfe, Panel2 = "Elapsed time"))
               
fig9 <- ggplot(br5df, aes(x = tvec)) + theme_bw() + ylim(0, 1) + ijarcol + rightlegend + #ijarfill + 
  geom_line(aes(y = upper-lower, group = interaction(Item, Panel), colour = Item, lty = Panel)) +
  facet_wrap(~Panel2, nrow=1) + xlab(expression(t)) + ylab(bquote(bar(R)[sys](t) - ~underline(R)[sys](t)))
fig9

#setEPS()
#postscript("fig9.eps",width=8,height=6)
pdf("fig9.pdf", width=8, height=3)
fig9
dev.off()


#