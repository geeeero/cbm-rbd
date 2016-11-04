# --------------------------------------------------------------------------- #
# cbm systems paper - plots for campi slotevent talk
# --------------------------------------------------------------------------- #

source("weibull-sys-functions.R")
source("plotdefs.R")

# component reliability plots
sysA <- graph.formula(s -- A -- t)
betaA <- 2
compA <- list(survsign = computeSystemSurvivalSignature(setCompTypes(sysA, list("A" = "A"))),
              n0y0 = list(A = c(1, failuretolambda(5, betaA))),
              beta = betaA)
compArel0 <-sysrelnowhor(compA$survsign, compA$n0y0, compA$beta, list(A = NULL), tnow = 0, hor = 10, seqlen=101)
compArel5 <-sysrelnowhor(compA$survsign, compA$n0y0, compA$beta, list(A = 5), tnow = 5, hor = 10, seqlen=101)
pdf("talk-campi/compArel0.pdf", width = minirelsize, height = minirelsize)
ggplot(compArel0, aes(x = t, y = rel)) + geom_line() + minirel
dev.off()
pdf("talk-campi/compArel5.pdf", width = minirelsize, height = minirelsize)
ggplot(compArel5, aes(x = t, y = rel)) + geom_line() + minirel
dev.off()

sysB <- graph.formula(s -- B -- t)
betaB <- 1.2
compB <- list(survsign = computeSystemSurvivalSignature(setCompTypes(sysB, list("B" = "B"))),
              n0y0 = list(B = c(5, failuretolambda(8, betaB))),
              beta = betaB)
compBrel0 <-sysrelnowhor(compB$survsign, compB$n0y0, compB$beta, list(B = NULL), tnow = 0, hor = 10, seqlen=101)
compBrel5 <-sysrelnowhor(compB$survsign, compB$n0y0, compB$beta, list(B = NULL), tnow = 5, hor = 10, seqlen=101)
pdf("talk-campi/compBrel0.pdf", width = minirelsize, height = minirelsize)
ggplot(compBrel0, aes(x = t, y = rel)) + geom_line() + minirel
dev.off()
pdf("talk-campi/compBrel5.pdf", width = minirelsize, height = minirelsize)
ggplot(compBrel5, aes(x = t, y = rel)) + geom_line() + minirel
dev.off()


# AB system reliability plots
sysAB <- graph.formula(s -- A1:A2:B -- t)
sysABtypes <- list("A"=c("A1", "A2"), "B"=c("B"))
sysAB <- setCompTypes(sysAB, sysABtypes)
sysABsign <- computeSystemSurvivalSignature(sysAB)
sysABfts <- list(A1 = NULL, A2 = 5, B = NULL)
sysABn0y0 <- list(A = compA$n0y0$A, B = compB$n0y0$B)
sysABrelhist05 <- sysrelnowhist(sysAB, sysABtypes, sysABfts, sysABn0y0, c(betaA, betaB), c(0,5), hor=15)
pdf("talk-campi/sysABrel0.pdf", width = 1.5*minirelsize, height = minirelsize)
ggplot(subset(sysABrelhist05, tnow == 0), aes(x = t, y = rel)) + geom_line() + minirel
dev.off()
pdf("talk-campi/sysABrel5.pdf", width = 1.5*minirelsize, height = minirelsize)
ggplot(subset(sysABrelhist05, tnow == 5), aes(x = t, y = rel)) + geom_line() + minirel
dev.off()

sysABghist <- gnowhist(sysAB, sysABtypes, sysABfts, sysABn0y0, c(betaA, betaB), seq(0, 10, by=2.5), hor = 15, seqlen = 301)
sysABghist$tnow <- as.factor(sysABghist$tnow)

pdf("talk-campi/sysABrelhist5.pdf", width = 5.3, height = 1.6)
ggplot(sysABghist, aes(x = t, y = rel)) + coord_cartesian(xlim = c(0,25), ylim = c(0,1)) + xlab(expression(t)) +
  geom_line(aes(group = tnow, colour = tnow)) + #ylab(expression(paste(R[sys]^(t[now]), (t),))) +
  theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  #guides(colour=guide_legend(title=expression(t[now]))) + 
  guides(colour = 'none') +
  scale_y_continuous(breaks=c(0, 0.5, 1))
dev.off()

pdf("talk-campi/sysABghist5.pdf", width = 6, height = 1.7)
ggplot(sysABghist, aes(x = t, y = gnow)) + coord_cartesian(xlim = c(0,25), ylim = c(0,0.4)) + xlab(expression(t)) +
  geom_line(aes(group = tnow, colour = tnow)) + #ylab(expression(paste(g^(t[now]), (t),))) +
  theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(axis.title.y = element_blank()) + guides(colour=guide_legend(title=expression(t[now]))) +
  scale_y_continuous(breaks=seq(0.1, 0.4, 0.1))
dev.off()

sysABtauhist <- taustarhist(sysAB, sysABtypes, sysABfts, sysABn0y0, c(betaA, betaB), seq(0, 10, by=2.5), hor = 5, seqlen = 501)
sysABtauhist$tnowf <- as.factor(sysABtauhist$tnow)

pdf("talk-campi/sysABtau5.pdf", width = 4.5, height = 1)
ggplot(sysABtauhist, aes(x = tnow, y = tstar, ymin = tnow, ymax = tstar, colour = tnowf)) + geom_linerange(size=1.5) +
  scale_x_continuous(breaks=seq(0, 10, 2.5), minor_breaks = NULL) + coord_flip(ylim = c(0,25)) + ylab("t") +
  #guides(colour=guide_legend(title=expression(t[now]))) +
  guides(colour = 'none') +
  theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
dev.off()

# more complex system?


#